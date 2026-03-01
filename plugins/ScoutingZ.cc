// ScoutingZ.cc
// Standard candles in Run3 Scouting: Z -> mumu, J/psi -> mumu (all/prompt/nonprompt)
//
// Refactor (2026-02):
// - W removed completely (MET left for future, but no code paths).
// - Tag&Probe kept formally strict: Tag does NOT require global; global is measured on probe.
// - Tag tightened on tracker: >=4 tracker layers, >=1 pixel layer, chi2/ndof cut.
// - Probe remains inclusive (pt>=3 by default).
// - Z: main Tag pT = 15 GeV (configurable) + Tag pT scan (3/10/15/25 by default) in separate folders.
// - Vertexed categories added (require common vtxIndx between the two muons) in parallel folders.
// - passMuonSystem_* removed entirely (and SA pass/fail histograms removed).
// - muonVertexCollection is optional and handled safely; needed only for J/psi prompt/nonprompt classification.
// - Track branch "TrkSta" renamed by meaning: "TrkMuMatch" (track matched to any scouting muon) and "TrkMuGlobal".
//
// Build: EDM plugin.

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"

#include <map>
#include <unordered_set>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <utility>

namespace {
  constexpr double MUON_MASS  = 0.1056583755;
  constexpr double MUON_BOUND = 0.9; // keep as-is: BB if both |eta|<0.9; EE if both >=0.9; else BE
  constexpr double MZ = 91.1876;
  constexpr double MJPSI = 3.0969;

  inline bool isCentral(double eta) { return std::abs(eta) < MUON_BOUND; }

  inline float relIso(const Run3ScoutingMuon& mu) {
    const float absIso = mu.trackIso() + mu.ecalIso() + mu.hcalIso();
    return absIso / std::max(1.f, mu.pt());
  }

  enum class Region { BB=0, BE=1, EE=2, N=3 };

  inline Region pairRegion(bool cen1, bool cen2) {
    if (cen1 && cen2) return Region::BB;
    if (!cen1 && !cen2) return Region::EE;
    return Region::BE;
  }

  // ---- J/psi class ----
  enum class JpsiClass { Unknown=0, Prompt=1, NonPrompt=2 };

  // ---- SFINAE helpers to get chi2/ndof if the method names exist ----
  template <typename T, typename = void> struct has_tk_chi2_ndof : std::false_type {};
  template <typename T>
  struct has_tk_chi2_ndof<T, std::void_t<
      decltype(std::declval<T>().tk_chi2()),
      decltype(std::declval<T>().tk_ndof())
    >> : std::true_type {};

  inline double safeChi2OverNdof(double chi2, double ndof) {
    return chi2 / std::max(1.0, ndof);
  }

  inline double muTrkChi2OverNdof(const Run3ScoutingMuon& mu) {
    // In Run3ScoutingMuon the accessors are trk_chi2 / trk_ndof
    return safeChi2OverNdof(mu.trk_chi2(), mu.trk_ndof());
  }

  inline double trkChi2OverNdof(const Run3ScoutingTrack& trk) {
    if constexpr (has_tk_chi2_ndof<Run3ScoutingTrack>::value) {
      return safeChi2OverNdof(trk.tk_chi2(), trk.tk_ndof());
    } else {
      // most common case in scouting: only tk_chi2() exists -> treat as chi2 (no ndof)
      return trk.tk_chi2();
    }
  }

  inline bool haveCommonVtxIndx(const Run3ScoutingMuon& a, const Run3ScoutingMuon& b) {
    const auto& ia = a.vtxIndx();
    const auto& ib = b.vtxIndx();
    if (ia.empty() || ib.empty()) return false;
    std::unordered_set<int> sa(ia.begin(), ia.end());
    for (int x : ib) {
      if (sa.count(x)) return true;
    }
    return false;
  }

  inline const char* regionName(Region r) {
    switch (r) {
      case Region::BB: return "BB";
      case Region::BE: return "BE";
      case Region::EE: return "EE";
      default: return "UNK";
    }
  }

  inline int regionIndex(Region r) {
    return (r==Region::BB ? 0 : (r==Region::BE ? 1 : 2));
  }

  inline std::string ptLabel(double pt) {
    // folder-friendly label: pt3, pt10, pt15, pt25
    const int ipt = static_cast<int>(std::lround(pt));
    return "pt" + std::to_string(ipt);
  }
}

class ScoutingZ : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingZ(const edm::ParameterSet&);
  ~ScoutingZ() override = default;

  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

private:
  // ---- Inputs ----
  edm::EDGetTokenT<std::vector<Run3ScoutingMuon>>    muonToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> pvToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingTrack>>  trackToken_;

  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> muonVtxToken_; // displacedVtx for prompt/nonprompt
  bool hasMuonVtxToken_{false};

  bool doZ_{true};
  bool doJpsi_{true};

  // ---- Tag/probe base ----
  double tagEtaMax_{2.4};
  double tagTrkChi2NdofMax_{3.0};
  int    tagMinTrkLayers_{4};
  int    tagMinPixLayers_{1};

  double probePtMin_{3.0};
  double probeEtaMax_{2.4};

  // ---- Global measured component (NOT required in tag) ----
  double gloNormChi2Max_{999.0};  // only used if isGlobalMuon==true

  // ---- Z tag pT ----
  double zTagPtMain_{15.0};
  std::vector<double> zTagPtScan_{3.0, 10.0, 15.0, 25.0};

  // ---- Track probe selection ----
  double trkPtMin_{3.0};
  double trkEtaMax_{2.4};
  double trkChi2NdofMax_{3.0};
  int    trkMinLayers_{4};
  double dr2Match_{0.09}; // deltaR^2 match track <-> muon

  // ---- Mass windows / binning ----
  // Z working window: 60-120
  int    zMassBins_{60};
  double zMassMin_{60.0};
  double zMassMax_{120.0};
  // Z control mass wide: 40-160 (1D only)
  int    zMassBinsCtrl_{120};
  double zMassMinCtrl_{40.0};
  double zMassMaxCtrl_{160.0};

  // J/psi window
  int    jMassBins_{50};
  double jMassMin_{2.6};
  double jMassMax_{3.6};

  // Prompt split
  double lxySigPromptMax_{3.0};

  // PV selection (optional control)
  double vtxNdofMin_{4.0};
  double vtxAbsZMax_{24.0};
  double vtxRhoMax_{2.0};

  // J/psi candidate selection controls
  unsigned jpsiTopN_{2};          // keep up to N disjoint pairs per event
  bool jpsiUseDisjointPairs_{true};

  // ---- Services ----
  edm::Service<TFileService> fs_;

  // ---- Buffers (per LS) ----
  using Buf = std::map<unsigned, std::vector<double>>;

  struct PerRegionBuf {
    Buf bb, be, ee;
  };

  struct PerRegionH2 {
    TH2D* bb{nullptr}; TH2D* be{nullptr}; TH2D* ee{nullptr};
  };

  static void fillPerRegion(PerRegionBuf& r, unsigned ls, Region reg, double val) {
    if (reg == Region::BB) r.bb[ls].push_back(val);
    else if (reg == Region::EE) r.ee[ls].push_back(val);
    else r.be[ls].push_back(val);
  }

  static inline void flushMapToH2(TH2D* h, const Buf& m) {
    if (!h) return;
    for (const auto& [ls, vals] : m) {
      for (double v : vals) h->Fill(static_cast<double>(ls), v);
    }
  }

  static inline void flushPerRegion(const PerRegionH2& h2, const PerRegionBuf& buf) {
    flushMapToH2(h2.bb, buf.bb);
    flushMapToH2(h2.be, buf.be);
    flushMapToH2(h2.ee, buf.ee);
  }

  // -----------------------
  // Z packs: inclusive + vertexed
  // -----------------------
  struct ZBufPack {
    PerRegionBuf twoTag, oneTag;
    PerRegionBuf id_pass, id_fail;
    PerRegionBuf glo_pass, glo_fail;

    // track branch: tag mu + track probe
    PerRegionBuf trkMu_pass, trkMu_fail;       // track matched to any scouting muon
    PerRegionBuf trkGlo_pass, trkGlo_fail;     // matched muon is global (measured)

    // Vertexed (common vtxIndx between tag and probe muons) - muon-muon only
    PerRegionBuf v_twoTag, v_oneTag;
    PerRegionBuf v_id_pass, v_id_fail;
    PerRegionBuf v_glo_pass, v_glo_fail;
  };

  struct ZH2Pack {
    PerRegionH2 h2tag, h1tag;
    PerRegionH2 hidp, hidf;
    PerRegionH2 hglop, hglof;

    PerRegionH2 htrkMu_p, htrkMu_f;
    PerRegionH2 htrkGlo_p, htrkGlo_f;

    PerRegionH2 vh2tag, vh1tag;
    PerRegionH2 vhidp, vhidf;
    PerRegionH2 vhglop, vhglof;
  };

  ZBufPack zMainBuf_;
  std::vector<ZBufPack> zScanBuf_;   // size = zTagPtScan_.size()

  ZH2Pack zMainH2_;
  std::vector<ZH2Pack> zScanH2_;

  // -----------------------
  // J/psi packs: 3 classes, inclusive + vertexed
  // -----------------------
  struct JpsiCats {
    PerRegionBuf twoTag, oneTag;
    PerRegionBuf id_pass, id_fail;
    PerRegionBuf glo_pass, glo_fail;

    PerRegionBuf trkMu_pass, trkMu_fail;
    PerRegionBuf trkGlo_pass, trkGlo_fail;

    PerRegionBuf v_twoTag, v_oneTag;
    PerRegionBuf v_id_pass, v_id_fail;
    PerRegionBuf v_glo_pass, v_glo_fail;
  };

  struct JH2Pack {
    PerRegionH2 h2tag, h1tag;
    PerRegionH2 hidp, hidf;
    PerRegionH2 hglop, hglof;

    PerRegionH2 htrkMu_p, htrkMu_f;
    PerRegionH2 htrkGlo_p, htrkGlo_f;

    PerRegionH2 vh2tag, vh1tag;
    PerRegionH2 vhidp, vhidf;
    PerRegionH2 vhglop, vhglof;
  };

  JpsiCats j_all_, j_prompt_, j_nonprompt_;
  JH2Pack  jAllH2_, jPrH2_, jNpH2_;

  // PV per-LS (control)
  std::map<unsigned, std::vector<unsigned>> z_npv_, j_npv_;

  // ---- Control histos (1D only) ----
  // Z control
  TH1I* h_z_cutflow_{nullptr};
  TH1D* h_z_mass_ctrl_{nullptr};

  // keep tag control per region (pair-region), but remove muon-system histos
  TH1D* h_z_tag_pt_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_eta_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_relIso_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_trkChi2Ndof_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_nTrkLay_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_nPixLay_[3]{nullptr,nullptr,nullptr};
  TH1I* h_z_tag_isGlobal_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_normChi2_[3]{nullptr,nullptr,nullptr}; // filled only if isGlobal

  // J/psi control
  TH1I* h_j_cutflow_{nullptr};
  TH1D* h_j_lxy_{nullptr};
  TH1D* h_j_lxySig_{nullptr};

  // ---- Helpers: selections ----
  bool passPV_(const Run3ScoutingVertex& v) const {
    if (!v.isValidVtx()) return false;
    if (v.ndof() < vtxNdofMin_) return false;
    if (std::abs(v.z()) > vtxAbsZMax_) return false;
    const double rho = std::hypot(v.x(), v.y());
    if (rho > vtxRhoMax_) return false;
    return true;
  }

  unsigned countGoodPVs_(const std::vector<Run3ScoutingVertex>& pvs) const {
    unsigned n=0;
    for (const auto& v : pvs) if (passPV_(v)) ++n;
    return n;
  }

  bool passProbeBase_(const Run3ScoutingMuon& mu) const {
    if (mu.pt() < probePtMin_) return false;
    if (std::abs(mu.eta()) > probeEtaMax_) return false;
    return true;
  }

  bool passTrackerQuality_(const Run3ScoutingMuon& mu) const {
    if (mu.nTrackerLayersWithMeasurement() < tagMinTrkLayers_) return false;
    if (mu.nPixelLayersWithMeasurement()   < tagMinPixLayers_) return false;

    const double chi2ndof = muTrkChi2OverNdof(mu);
    if (chi2ndof >= 0.0 && chi2ndof > tagTrkChi2NdofMax_) return false;
    return true;
  }

  bool passTagWithPt_(const Run3ScoutingMuon& mu, double ptMin) const {
    if (mu.pt() < ptMin) return false;
    if (std::abs(mu.eta()) > tagEtaMax_) return false;
    if (!passTrackerQuality_(mu)) return false;
    return true;
  }

  bool passGlobalMeasured_(const Run3ScoutingMuon& mu) const {
    if (!mu.isGlobalMuon()) return false;
    if (mu.normalizedChi2() > gloNormChi2Max_) return false;
    return true;
  }

  bool passTrackProbe_(const Run3ScoutingTrack& trk) const {
    if (trk.tk_pt() < trkPtMin_) return false;
    if (std::abs(trk.tk_eta()) > trkEtaMax_) return false;
    if (trk.tk_nTrackerLayersWithMeasurement() < trkMinLayers_) return false;

    const double chi2ndof = trkChi2OverNdof(trk);
    if (chi2ndof > trkChi2NdofMax_) return false;
    return true;
  }

  // ---- J/psi prompt/nonprompt classifier ----
  JpsiClass classifyJpsi_(const Run3ScoutingMuon& mu1,
                          const Run3ScoutingMuon& mu2,
                          const std::vector<Run3ScoutingVertex>& muVtx,
                          double avgX, double avgY) const
  {
    const auto& idx1 = mu1.vtxIndx();
    const auto& idx2 = mu2.vtxIndx();
    if (idx1.empty() || idx2.empty()) return JpsiClass::Unknown;

    std::unordered_set<int> s1(idx1.begin(), idx1.end());
    for (int vidx : idx2) {
      if (!s1.count(vidx)) continue;
      if (vidx < 0 || static_cast<size_t>(vidx) >= muVtx.size()) continue;
      const auto& v = muVtx[static_cast<size_t>(vidx)];

      const double dx = v.x() - avgX;
      const double dy = v.y() - avgY;
      const double r2 = dx*dx + dy*dy;
      if (r2 <= 0.0) return JpsiClass::Unknown;

      const double Lxy = std::sqrt(r2);
      const double sx2 = double(v.xError()) * double(v.xError());
      const double sy2 = double(v.yError()) * double(v.yError());
      const double num = dx*dx*sx2 + dy*dy*sy2;
      const double LxyErr = (num > 0.0) ? std::sqrt(num) / Lxy : -1.0;
      const double LxySig = (LxyErr > 0.0) ? (Lxy / LxyErr) : -1.0;

      if (h_j_lxy_)    h_j_lxy_->Fill(Lxy);
      if (h_j_lxySig_) h_j_lxySig_->Fill(LxySig);

      if (LxyErr <= 0.0) return JpsiClass::Unknown;
      return (LxySig < lxySigPromptMax_) ? JpsiClass::Prompt : JpsiClass::NonPrompt;
    }
    return JpsiClass::Unknown;
  }

  // ---- Booking helpers (endJob) ----
  static void scanLS(const Buf& b, unsigned& minLS, unsigned& maxLS) {
    for (const auto& kv : b) { minLS = std::min(minLS, kv.first); maxLS = std::max(maxLS, kv.first); }
  }
  static void scanLS(const PerRegionBuf& r, unsigned& minLS, unsigned& maxLS) {
    scanLS(r.bb, minLS, maxLS); scanLS(r.be, minLS, maxLS); scanLS(r.ee, minLS, maxLS);
  }

  TH2D* makeLSMass2D(TFileDirectory& dir,
                    const std::string& name, const std::string& title,
                    unsigned minLS, unsigned maxLS, int ybins, double ymin, double ymax) const {
    const int nLumiBins = (maxLS >= minLS) ? int(maxLS - minLS + 1) : 1;
    const double lmin = double(minLS) - 0.5;
    const double lmax = double(maxLS) + 0.5;
    auto* h = dir.make<TH2D>(name.c_str(), title.c_str(),
                            nLumiBins, lmin, lmax, ybins, ymin, ymax);
    h->GetXaxis()->SetTitle("luminosity section");
    h->GetYaxis()->SetTitle("mass [GeV]");
    return h;
  }

  void bookPerRegionH2(PerRegionH2& out,
                      TFileDirectory& bbDir, TFileDirectory& beDir, TFileDirectory& eeDir,
                      const std::string& baseName, const std::string& baseTitle,
                      unsigned minLS, unsigned maxLS, int ybins, double ymin, double ymax) const {
    out.bb = makeLSMass2D(bbDir, baseName, baseTitle + " BB", minLS, maxLS, ybins, ymin, ymax);
    out.be = makeLSMass2D(beDir, baseName, baseTitle + " BE", minLS, maxLS, ybins, ymin, ymax);
    out.ee = makeLSMass2D(eeDir, baseName, baseTitle + " EE", minLS, maxLS, ybins, ymin, ymax);
  }

  void bookZPack(ZH2Pack& out, TFileDirectory& baseDir,
                 unsigned minLS, unsigned maxLS,
                 int ybins, double ymin, double ymax,
                 bool withTrackBranch,
                 bool withVertexed) const
  {
    auto bb = baseDir.mkdir("BB");
    auto be = baseDir.mkdir("BE");
    auto ee = baseDir.mkdir("EE");

    bookPerRegionH2(out.h2tag, bb, be, ee, "h_mass_2TAG",     "Z 2TAG",       minLS, maxLS, ybins, ymin, ymax);
    bookPerRegionH2(out.h1tag, bb, be, ee, "h_mass_1TAG",     "Z 1TAG",       minLS, maxLS, ybins, ymin, ymax);

    bookPerRegionH2(out.hidp,  bb, be, ee, "h_mass_ID_pass",  "Z ID pass",    minLS, maxLS, ybins, ymin, ymax);
    bookPerRegionH2(out.hidf,  bb, be, ee, "h_mass_ID_fail",  "Z ID fail",    minLS, maxLS, ybins, ymin, ymax);

    bookPerRegionH2(out.hglop, bb, be, ee, "h_mass_Glo_pass", "Z Global pass",minLS, maxLS, ybins, ymin, ymax);
    bookPerRegionH2(out.hglof, bb, be, ee, "h_mass_Glo_fail", "Z Global fail",minLS, maxLS, ybins, ymin, ymax);

    if (withTrackBranch) {
      bookPerRegionH2(out.htrkMu_p, bb, be, ee, "h_mass_TrkMuMatch_pass", "Z TrkMuMatch pass", minLS, maxLS, ybins, ymin, ymax);
      bookPerRegionH2(out.htrkMu_f, bb, be, ee, "h_mass_TrkMuMatch_fail", "Z TrkMuMatch fail", minLS, maxLS, ybins, ymin, ymax);
      bookPerRegionH2(out.htrkGlo_p,bb, be, ee, "h_mass_TrkMuGlobal_pass","Z TrkMuGlobal pass",minLS, maxLS, ybins, ymin, ymax);
      bookPerRegionH2(out.htrkGlo_f,bb, be, ee, "h_mass_TrkMuGlobal_fail","Z TrkMuGlobal fail",minLS, maxLS, ybins, ymin, ymax);
    }

    if (withVertexed) {
      auto vtxDir = baseDir.mkdir("Vtx");
      auto vbb = vtxDir.mkdir("BB");
      auto vbe = vtxDir.mkdir("BE");
      auto vee = vtxDir.mkdir("EE");

      bookPerRegionH2(out.vh2tag, vbb, vbe, vee, "h_mass_2TAG",     "Z(Vtx) 2TAG",       minLS, maxLS, ybins, ymin, ymax);
      bookPerRegionH2(out.vh1tag, vbb, vbe, vee, "h_mass_1TAG",     "Z(Vtx) 1TAG",       minLS, maxLS, ybins, ymin, ymax);

      bookPerRegionH2(out.vhidp,  vbb, vbe, vee, "h_mass_ID_pass",  "Z(Vtx) ID pass",    minLS, maxLS, ybins, ymin, ymax);
      bookPerRegionH2(out.vhidf,  vbb, vbe, vee, "h_mass_ID_fail",  "Z(Vtx) ID fail",    minLS, maxLS, ybins, ymin, ymax);

      bookPerRegionH2(out.vhglop, vbb, vbe, vee, "h_mass_Glo_pass", "Z(Vtx) Global pass",minLS, maxLS, ybins, ymin, ymax);
      bookPerRegionH2(out.vhglof, vbb, vbe, vee, "h_mass_Glo_fail", "Z(Vtx) Global fail",minLS, maxLS, ybins, ymin, ymax);
    }
  }

  void bookJPack(JH2Pack& out, TFileDirectory& baseDir,
                 const std::string& titlePrefix,
                 unsigned minLS, unsigned maxLS,
                 int ybins, double ymin, double ymax,
                 bool withTrackBranch,
                 bool withVertexed) const
  {
    auto bb = baseDir.mkdir("BB");
    auto be = baseDir.mkdir("BE");
    auto ee = baseDir.mkdir("EE");

    bookPerRegionH2(out.h2tag, bb, be, ee, "h_mass_2TAG",     titlePrefix + " 2TAG",        minLS, maxLS, ybins, ymin, ymax);
    bookPerRegionH2(out.h1tag, bb, be, ee, "h_mass_1TAG",     titlePrefix + " 1TAG",        minLS, maxLS, ybins, ymin, ymax);

    bookPerRegionH2(out.hidp,  bb, be, ee, "h_mass_ID_pass",  titlePrefix + " ID pass",     minLS, maxLS, ybins, ymin, ymax);
    bookPerRegionH2(out.hidf,  bb, be, ee, "h_mass_ID_fail",  titlePrefix + " ID fail",     minLS, maxLS, ybins, ymin, ymax);

    bookPerRegionH2(out.hglop, bb, be, ee, "h_mass_Glo_pass", titlePrefix + " Global pass", minLS, maxLS, ybins, ymin, ymax);
    bookPerRegionH2(out.hglof, bb, be, ee, "h_mass_Glo_fail", titlePrefix + " Global fail", minLS, maxLS, ybins, ymin, ymax);

    if (withTrackBranch) {
      bookPerRegionH2(out.htrkMu_p, bb, be, ee, "h_mass_TrkMuMatch_pass", titlePrefix + " TrkMuMatch pass", minLS, maxLS, ybins, ymin, ymax);
      bookPerRegionH2(out.htrkMu_f, bb, be, ee, "h_mass_TrkMuMatch_fail", titlePrefix + " TrkMuMatch fail", minLS, maxLS, ybins, ymin, ymax);
      bookPerRegionH2(out.htrkGlo_p,bb, be, ee, "h_mass_TrkMuGlobal_pass",titlePrefix + " TrkMuGlobal pass",minLS, maxLS, ybins, ymin, ymax);
      bookPerRegionH2(out.htrkGlo_f,bb, be, ee, "h_mass_TrkMuGlobal_fail",titlePrefix + " TrkMuGlobal fail",minLS, maxLS, ybins, ymin, ymax);
    }

    if (withVertexed) {
      auto vtxDir = baseDir.mkdir("Vtx");
      auto vbb = vtxDir.mkdir("BB");
      auto vbe = vtxDir.mkdir("BE");
      auto vee = vtxDir.mkdir("EE");

      bookPerRegionH2(out.vh2tag, vbb, vbe, vee, "h_mass_2TAG",     titlePrefix + "(Vtx) 2TAG",        minLS, maxLS, ybins, ymin, ymax);
      bookPerRegionH2(out.vh1tag, vbb, vbe, vee, "h_mass_1TAG",     titlePrefix + "(Vtx) 1TAG",        minLS, maxLS, ybins, ymin, ymax);

      bookPerRegionH2(out.vhidp,  vbb, vbe, vee, "h_mass_ID_pass",  titlePrefix + "(Vtx) ID pass",     minLS, maxLS, ybins, ymin, ymax);
      bookPerRegionH2(out.vhidf,  vbb, vbe, vee, "h_mass_ID_fail",  titlePrefix + "(Vtx) ID fail",     minLS, maxLS, ybins, ymin, ymax);

      bookPerRegionH2(out.vhglop, vbb, vbe, vee, "h_mass_Glo_pass", titlePrefix + "(Vtx) Global pass", minLS, maxLS, ybins, ymin, ymax);
      bookPerRegionH2(out.vhglof, vbb, vbe, vee, "h_mass_Glo_fail", titlePrefix + "(Vtx) Global fail", minLS, maxLS, ybins, ymin, ymax);
    }
  }

  // Fill helpers
  void fillZPack_(ZBufPack& b, unsigned ls, Region reg, double mass,
                 const Run3ScoutingMuon& tag, const Run3ScoutingMuon& probe,
                 double tagPtMin, bool commonVtx) const
  {
    if (!passTagWithPt_(tag, tagPtMin)) return;

    const bool probeTight = passTagWithPt_(probe, tagPtMin);
    if (probeTight) fillPerRegion(b.twoTag, ls, reg, mass);
    else           fillPerRegion(b.oneTag, ls, reg, mass);

    const bool idok = passTrackerQuality_(probe);
    if (idok) fillPerRegion(b.id_pass, ls, reg, mass);
    else      fillPerRegion(b.id_fail, ls, reg, mass);

    const bool glok = passGlobalMeasured_(probe);
    if (glok) fillPerRegion(b.glo_pass, ls, reg, mass);
    else      fillPerRegion(b.glo_fail, ls, reg, mass);

    if (commonVtx) {
      if (probeTight) fillPerRegion(b.v_twoTag, ls, reg, mass);
      else           fillPerRegion(b.v_oneTag, ls, reg, mass);

      if (idok) fillPerRegion(b.v_id_pass, ls, reg, mass);
      else      fillPerRegion(b.v_id_fail, ls, reg, mass);

      if (glok) fillPerRegion(b.v_glo_pass, ls, reg, mass);
      else      fillPerRegion(b.v_glo_fail, ls, reg, mass);
    }
  }

  void fillJPack_(JpsiCats& b, unsigned ls, Region reg, double mass,
                 const Run3ScoutingMuon& tag, const Run3ScoutingMuon& probe,
                 bool commonVtx) const
  {
    // J/psi tag pT is not scanned here: tag definition is passTagWithPt_(..., pt>=probePtMin_).
    // You can later add a J/psi tag-scan similarly if needed.
    const double tagPtMin = probePtMin_;
    if (!passTagWithPt_(tag, tagPtMin)) return;

    const bool probeTight = passTagWithPt_(probe, tagPtMin);
    if (probeTight) fillPerRegion(b.twoTag, ls, reg, mass);
    else           fillPerRegion(b.oneTag, ls, reg, mass);

    const bool idok = passTrackerQuality_(probe);
    if (idok) fillPerRegion(b.id_pass, ls, reg, mass);
    else      fillPerRegion(b.id_fail, ls, reg, mass);

    const bool glok = passGlobalMeasured_(probe);
    if (glok) fillPerRegion(b.glo_pass, ls, reg, mass);
    else      fillPerRegion(b.glo_fail, ls, reg, mass);

    if (commonVtx) {
      if (probeTight) fillPerRegion(b.v_twoTag, ls, reg, mass);
      else           fillPerRegion(b.v_oneTag, ls, reg, mass);

      if (idok) fillPerRegion(b.v_id_pass, ls, reg, mass);
      else      fillPerRegion(b.v_id_fail, ls, reg, mass);

      if (glok) fillPerRegion(b.v_glo_pass, ls, reg, mass);
      else      fillPerRegion(b.v_glo_fail, ls, reg, mass);
    }
  }
};

ScoutingZ::ScoutingZ(const edm::ParameterSet& cfg) {
  usesResource("TFileService");

  muonToken_  = consumes<std::vector<Run3ScoutingMuon>>(   cfg.getParameter<edm::InputTag>("muonCollection") );
  pvToken_    = consumes<std::vector<Run3ScoutingVertex>>( cfg.getParameter<edm::InputTag>("primaryVertexCollection") );
  trackToken_ = consumes<std::vector<Run3ScoutingTrack>>(  cfg.getParameter<edm::InputTag>("trackCollection") );

  if (cfg.existsAs<edm::InputTag>("muonVertexCollection")) {
    muonVtxToken_ = consumes<std::vector<Run3ScoutingVertex>>( cfg.getParameter<edm::InputTag>("muonVertexCollection") );
    hasMuonVtxToken_ = true;
  }

  doZ_    = cfg.getUntrackedParameter<bool>("doZ", true);
  doJpsi_ = cfg.getUntrackedParameter<bool>("doJpsi", true);

  // Tag/probe
  tagEtaMax_ = cfg.getUntrackedParameter<double>("tagEtaMax", 2.4);
  tagTrkChi2NdofMax_ = cfg.getUntrackedParameter<double>("tagTrkChi2NdofMax", 3.0);
  tagMinTrkLayers_ = cfg.getUntrackedParameter<int>("tagMinTrackerLayers", 4);
  tagMinPixLayers_ = cfg.getUntrackedParameter<int>("tagMinPixelLayers", 1);

  probePtMin_ = cfg.getUntrackedParameter<double>("probePtMin", 3.0);
  probeEtaMax_ = cfg.getUntrackedParameter<double>("probeEtaMax", 2.4);

  gloNormChi2Max_ = cfg.getUntrackedParameter<double>("gloNormChi2Max", 999.0);

  // Z tag pT main + scan
  zTagPtMain_ = cfg.getUntrackedParameter<double>("zTagPtMain", 15.0);
  zTagPtScan_ = cfg.getUntrackedParameter<std::vector<double>>("zTagPtScan", zTagPtScan_);
  if (zTagPtScan_.empty()) zTagPtScan_ = {3.0, 10.0, 15.0, 25.0};

  // Track probe
  trkPtMin_ = cfg.getUntrackedParameter<double>("trkPtMin", 3.0);
  trkEtaMax_ = cfg.getUntrackedParameter<double>("trkEtaMax", 2.4);
  trkChi2NdofMax_ = cfg.getUntrackedParameter<double>("trkChi2NdofMax", 3.0);
  trkMinLayers_ = cfg.getUntrackedParameter<int>("trkMinLayers", 4);
  dr2Match_ = cfg.getUntrackedParameter<double>("dr2Match", 0.09);

  // Mass windows
  zMassBins_ = cfg.getUntrackedParameter<int>("zMassBins", 60);
  zMassMin_  = cfg.getUntrackedParameter<double>("zMassMin", 60.0);
  zMassMax_  = cfg.getUntrackedParameter<double>("zMassMax", 120.0);

  zMassBinsCtrl_ = cfg.getUntrackedParameter<int>("zMassBinsCtrl", 120);
  zMassMinCtrl_  = cfg.getUntrackedParameter<double>("zMassMinCtrl", 40.0);
  zMassMaxCtrl_  = cfg.getUntrackedParameter<double>("zMassMaxCtrl", 160.0);

  jMassBins_ = cfg.getUntrackedParameter<int>("jpsiMassBins", 50);
  jMassMin_  = cfg.getUntrackedParameter<double>("jpsiMassMin", 2.6);
  jMassMax_  = cfg.getUntrackedParameter<double>("jpsiMassMax", 3.6);

  lxySigPromptMax_ = cfg.getUntrackedParameter<double>("lxySigPromptMax", 3.0);

  // PV
  vtxNdofMin_ = cfg.getUntrackedParameter<double>("vtxNdofMin", 4.0);
  vtxAbsZMax_ = cfg.getUntrackedParameter<double>("vtxAbsZMax", 24.0);
  vtxRhoMax_  = cfg.getUntrackedParameter<double>("vtxRhoMax", 2.0);

  // J/psi pairing controls
  jpsiTopN_ = cfg.getUntrackedParameter<unsigned>("jpsiTopN", 2);
  jpsiUseDisjointPairs_ = cfg.getUntrackedParameter<bool>("jpsiUseDisjointPairs", true);

  // allocate scan packs
  zScanBuf_.assign(zTagPtScan_.size(), ZBufPack{});
  zScanH2_.assign(zTagPtScan_.size(), ZH2Pack{});
}

void ScoutingZ::beginJob() {
  if (!fs_.isAvailable()) return;

  // Z/Control
  auto zdir = fs_->mkdir("Z");
  auto zc   = zdir.mkdir("Control");

  h_z_cutflow_ = zc.make<TH1I>("cutflow", "Z cutflow", 10, 0.5, 10.5);
  h_z_cutflow_->GetXaxis()->SetBinLabel(1, "events");
  h_z_cutflow_->GetXaxis()->SetBinLabel(2, ">=2 mu");
  h_z_cutflow_->GetXaxis()->SetBinLabel(3, ">=1 tagQ (quality)");
  h_z_cutflow_->GetXaxis()->SetBinLabel(4, "OS pair");
  h_z_cutflow_->GetXaxis()->SetBinLabel(5, "mass ctrl (40-160)");
  h_z_cutflow_->GetXaxis()->SetBinLabel(6, "mass work (Z win)");

  h_z_mass_ctrl_ = zc.make<TH1D>("mass_ctrl", "Z mass control; m_{#mu#mu} [GeV]", zMassBinsCtrl_, zMassMinCtrl_, zMassMaxCtrl_);

  const char* rname[3] = {"BB","BE","EE"};
  for (int ir=0; ir<3; ++ir) {
    h_z_tag_pt_[ir]  = zc.make<TH1D>(std::string("tag_pt_"+std::string(rname[ir])).c_str(),
                                     ("tag pT "+std::string(rname[ir])+"; p_{T} [GeV]").c_str(), 200, 0, 200);
    h_z_tag_eta_[ir] = zc.make<TH1D>(std::string("tag_eta_"+std::string(rname[ir])).c_str(),
                                     ("tag eta "+std::string(rname[ir])+"; #eta").c_str(), 60, -3, 3);
    h_z_tag_relIso_[ir] = zc.make<TH1D>(std::string("tag_relIso_"+std::string(rname[ir])).c_str(),
                                        ("tag relIso "+std::string(rname[ir])+"; relIso").c_str(), 100, 0, 2);
    h_z_tag_trkChi2Ndof_[ir] = zc.make<TH1D>(std::string("tag_trkChi2Ndof_"+std::string(rname[ir])).c_str(),
                                             ("tag trk #chi^{2}/ndof "+std::string(rname[ir])+"; #chi^{2}/ndof").c_str(), 120, 0, 12);
    h_z_tag_nTrkLay_[ir] = zc.make<TH1D>(std::string("tag_nTrkLayers_"+std::string(rname[ir])).c_str(),
                                         ("tag nTrackerLayers "+std::string(rname[ir])+"; n").c_str(), 30, -0.5, 29.5);
    h_z_tag_nPixLay_[ir] = zc.make<TH1D>(std::string("tag_nPixLayers_"+std::string(rname[ir])).c_str(),
                                         ("tag nPixelLayers "+std::string(rname[ir])+"; n").c_str(), 15, -0.5, 14.5);
    h_z_tag_isGlobal_[ir] = zc.make<TH1I>(std::string("tag_isGlobal_"+std::string(rname[ir])).c_str(),
                                          ("tag isGlobal "+std::string(rname[ir])+"; flag").c_str(), 2, -0.5, 1.5);
    h_z_tag_normChi2_[ir] = zc.make<TH1D>(std::string("tag_normChi2_"+std::string(rname[ir])).c_str(),
                                          ("tag normalizedChi2 (global only) "+std::string(rname[ir])+"; #chi^{2}/ndof").c_str(), 200, 0, 200);
  }

  // Jpsi/Control
  auto jdir = fs_->mkdir("Jpsi");
  auto jc   = jdir.mkdir("Control");
  h_j_cutflow_ = jc.make<TH1I>("cutflow", "J/psi cutflow", 8, 0.5, 8.5);
  h_j_cutflow_->GetXaxis()->SetBinLabel(1, "events");
  h_j_cutflow_->GetXaxis()->SetBinLabel(2, ">=2 mu");
  h_j_cutflow_->GetXaxis()->SetBinLabel(3, "OS pair");
  h_j_cutflow_->GetXaxis()->SetBinLabel(4, "mass win");

  h_j_lxy_    = jc.make<TH1D>("lxy", "L_{xy}; L_{xy}", 1000, 0.0, 1.0);
  h_j_lxySig_ = jc.make<TH1D>("lxySig", "L_{xy}/#sigma; L_{xy}/#sigma", 500, 0.0, 50.0);
}

void ScoutingZ::analyze(const edm::Event& ev, const edm::EventSetup&) {
  const unsigned ls = ev.luminosityBlock();

  edm::Handle<std::vector<Run3ScoutingMuon>> hMu;
  edm::Handle<std::vector<Run3ScoutingVertex>> hPV;
  edm::Handle<std::vector<Run3ScoutingTrack>> hTrk;

  ev.getByToken(muonToken_, hMu);
  ev.getByToken(pvToken_, hPV);
  ev.getByToken(trackToken_, hTrk);

  edm::Handle<std::vector<Run3ScoutingVertex>> hMuVtx;
  bool haveMuVtx = false;
  if (hasMuonVtxToken_) {
    ev.getByToken(muonVtxToken_, hMuVtx);
    haveMuVtx = hMuVtx.isValid() && !hMuVtx->empty();
  }

  //if (!hMu.isValid() || !hPV.isValid() || !hTrk.isValid()) return;
  if (hMu->empty() || hPV->empty() || hTrk->empty()) return;

  const unsigned npv = countGoodPVs_(*hPV);
  if (doZ_)    z_npv_[ls].push_back(npv);
  if (doJpsi_) j_npv_[ls].push_back(npv);

  // avg PV position for Lxy
  double avgX = 0.0, avgY = 0.0;
  if (!hPV->empty()) {
    avgX = std::accumulate(hPV->begin(), hPV->end(), 0.0, [](double a, const auto& v){return a+v.x();}) / hPV->size();
    avgY = std::accumulate(hPV->begin(), hPV->end(), 0.0, [](double a, const auto& v){return a+v.y();}) / hPV->size();
  }

  const size_t nMu = hMu->size();
  if (doZ_ && h_z_cutflow_) h_z_cutflow_->Fill(1);
  if (doJpsi_ && h_j_cutflow_) h_j_cutflow_->Fill(1);

  if (nMu < 2) return;
  if (doZ_ && h_z_cutflow_) h_z_cutflow_->Fill(2);
  if (doJpsi_ && h_j_cutflow_) h_j_cutflow_->Fill(2);

  // =========================
  // Z: muon-muon T&P, tag = tracker-quality + (pt threshold per category)
  // Probe base: inclusive pt>=probePtMin_
  // =========================
  if (doZ_) {
    TLorentzVector v1, v2;

    // collect tag-quality indices (no pT cut here, only quality+eta)
    std::vector<size_t> tagsQ;
    tagsQ.reserve(nMu);
    for (size_t i=0; i<nMu; ++i) {
      const auto& mu = hMu->at(i);
      if (std::abs(mu.eta()) > tagEtaMax_) continue;
      if (!passTrackerQuality_(mu)) continue;
      tagsQ.push_back(i);
    }
    if (!tagsQ.empty() && h_z_cutflow_) h_z_cutflow_->Fill(3);

    for (size_t it=0; it<tagsQ.size(); ++it) {
      const auto& tag = hMu->at(tagsQ[it]);
      v1.SetPtEtaPhiM(tag.pt(), tag.eta(), tag.phi(), MUON_MASS);
      const bool cenTag = isCentral(tag.eta());

      for (size_t j=0; j<nMu; ++j) {
        if (j == tagsQ[it]) continue;
        const auto& probe = hMu->at(j);
        if (!passProbeBase_(probe)) continue;
        if (tag.charge() == probe.charge()) continue;

        v2.SetPtEtaPhiM(probe.pt(), probe.eta(), probe.phi(), MUON_MASS);
        const double mass = (v1+v2).M();

        if (mass >= zMassMinCtrl_ && mass <= zMassMaxCtrl_) {
          if (h_z_cutflow_) h_z_cutflow_->Fill(5);
          if (h_z_mass_ctrl_) h_z_mass_ctrl_->Fill(mass);
        }

        if (mass < zMassMin_ || mass > zMassMax_) continue;
        if (h_z_cutflow_) h_z_cutflow_->Fill(6);
        if (h_z_cutflow_) h_z_cutflow_->Fill(4);

        const bool cenProbe = isCentral(probe.eta());
        const Region reg = pairRegion(cenTag, cenProbe);
        const int ir = regionIndex(reg);

        // control fill for tag (per pair region); note: tag-quality only, no pt requirement here
        if (h_z_tag_pt_[ir])  h_z_tag_pt_[ir]->Fill(tag.pt());
        if (h_z_tag_eta_[ir]) h_z_tag_eta_[ir]->Fill(tag.eta());
        if (h_z_tag_relIso_[ir]) h_z_tag_relIso_[ir]->Fill(relIso(tag));
        const double tchi = muTrkChi2OverNdof(tag);
        if (h_z_tag_trkChi2Ndof_[ir]) h_z_tag_trkChi2Ndof_[ir]->Fill(tchi >= 0.0 ? tchi : -0.1);
        if (h_z_tag_nTrkLay_[ir]) h_z_tag_nTrkLay_[ir]->Fill(tag.nTrackerLayersWithMeasurement());
        if (h_z_tag_nPixLay_[ir]) h_z_tag_nPixLay_[ir]->Fill(tag.nPixelLayersWithMeasurement());
        if (h_z_tag_isGlobal_[ir]) h_z_tag_isGlobal_[ir]->Fill(tag.isGlobalMuon() ? 1 : 0);
        if (tag.isGlobalMuon() && h_z_tag_normChi2_[ir]) h_z_tag_normChi2_[ir]->Fill(tag.normalizedChi2());

        const bool commonVtx = haveCommonVtxIndx(tag, probe);

        // main pack (default tagPt = 15 GeV)
        fillZPack_(zMainBuf_, ls, reg, mass, tag, probe, zTagPtMain_, commonVtx);

        // scan packs
        for (size_t k=0; k<zTagPtScan_.size(); ++k) {
          fillZPack_(zScanBuf_[k], ls, reg, mass, tag, probe, zTagPtScan_[k], commonVtx);
        }
      }
    }

    // =========================
    // Z: muon-track branch (tag with pt threshold, track probe), mass in Z window
    // Fill ONLY inclusive (no vertexed) for now (track has no vtxIndx in this context).
    // =========================
    auto fillZTrack = [&](ZBufPack& b, double tagPtMin) {
      for (size_t it=0; it<tagsQ.size(); ++it) {
        const auto& tag = hMu->at(tagsQ[it]);
        if (!passTagWithPt_(tag, tagPtMin)) continue;

        TLorentzVector vtag, vtrk;
        vtag.SetPtEtaPhiM(tag.pt(), tag.eta(), tag.phi(), MUON_MASS);
        const bool cenTag = isCentral(tag.eta());

        for (const auto& trk : *hTrk) {
          if (!passTrackProbe_(trk)) continue;
          if (trk.tk_charge() == tag.charge()) continue;

          vtrk.SetPtEtaPhiM(trk.tk_pt(), trk.tk_eta(), trk.tk_phi(), MUON_MASS);
          const double mass = (vtag+vtrk).M();
          if (mass < zMassMin_ || mass > zMassMax_) continue;

          const bool cenTrk = isCentral(trk.tk_eta());
          const Region reg = pairRegion(cenTag, cenTrk);

          // find matched muon
          int match = -1;
          for (size_t j=0; j<nMu; ++j) {
            const auto& mu = hMu->at(j);
            if (reco::deltaR2(mu.eta(), mu.phi(), trk.tk_eta(), trk.tk_phi()) < dr2Match_) { match = (int)j; break; }
          }

          if (match >= 0) {
            fillPerRegion(b.trkMu_pass, ls, reg, mass);
            const bool glok = passGlobalMeasured_(hMu->at((size_t)match));
            if (glok) fillPerRegion(b.trkGlo_pass, ls, reg, mass);
            else      fillPerRegion(b.trkGlo_fail, ls, reg, mass);
          } else {
            fillPerRegion(b.trkMu_fail, ls, reg, mass);
          }
        }
      }
    };

    // main + scan (keep consistent)
    fillZTrack(zMainBuf_, zTagPtMain_);
    for (size_t k=0; k<zTagPtScan_.size(); ++k) fillZTrack(zScanBuf_[k], zTagPtScan_[k]);
  }

  // =========================
  // J/psi: all/prompt/nonprompt, per-muon T&P with limited multi-cand selection
  // =========================
  if (doJpsi_) {
    TLorentzVector v1, v2;

    struct PairCand {
      size_t i, j;
      double mass;
      bool commonVtx;
      double score; // lower is better
    };

    std::vector<PairCand> cands;
    cands.reserve(nMu*nMu/2);

    for (size_t i=0; i<nMu; ++i) {
      const auto& mu1 = hMu->at(i);
      if (!passProbeBase_(mu1)) continue;
      v1.SetPtEtaPhiM(mu1.pt(), mu1.eta(), mu1.phi(), MUON_MASS);

      for (size_t j=i+1; j<nMu; ++j) {
        const auto& mu2 = hMu->at(j);
        if (!passProbeBase_(mu2)) continue;
        if (mu1.charge() == mu2.charge()) continue;

        v2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), MUON_MASS);
        const double mass = (v1+v2).M();
        if (mass < jMassMin_ || mass > jMassMax_) continue;

        if (h_j_cutflow_) h_j_cutflow_->Fill(4);
        if (h_j_cutflow_) h_j_cutflow_->Fill(3);

        const bool commonVtx = haveCommonVtxIndx(mu1, mu2);
        const double dm = std::abs(mass - MJPSI);

        // ranking: prefer common-vtx, then closeness to J/psi
        const double score = (commonVtx ? 0.0 : 1000.0) + dm;

        cands.push_back(PairCand{i, j, mass, commonVtx, score});
      }
    }

    if (!cands.empty()) {
      std::sort(cands.begin(), cands.end(), [](const PairCand& a, const PairCand& b){
        return a.score < b.score;
      });

      std::vector<PairCand> chosen;
      chosen.reserve(std::min<size_t>(jpsiTopN_, cands.size()));

      if (jpsiUseDisjointPairs_) {
        std::vector<char> used(nMu, 0);
        for (const auto& c : cands) {
          if (chosen.size() >= jpsiTopN_) break;
          if (used[c.i] || used[c.j]) continue;
          used[c.i] = 1; used[c.j] = 1;
          chosen.push_back(c);
        }
      } else {
        for (size_t k=0; k<std::min<size_t>(jpsiTopN_, cands.size()); ++k) chosen.push_back(cands[k]);
      }

      const bool canClassify = haveMuVtx && hMuVtx.isValid() && !hMuVtx->empty() && !hPV->empty();

      for (const auto& c : chosen) {
        const auto& mu1 = hMu->at(c.i);
        const auto& mu2 = hMu->at(c.j);

        JpsiClass jc = JpsiClass::Unknown;
        if (canClassify) jc = classifyJpsi_(mu1, mu2, *hMuVtx, avgX, avgY);

        const bool cen1 = isCentral(mu1.eta());
        const bool cen2 = isCentral(mu2.eta());
        const Region reg12 = pairRegion(cen1, cen2);
        const Region reg21 = pairRegion(cen2, cen1);

        // per-muon T&P: both orientations
        auto doFill = [&](JpsiCats& cat, const Run3ScoutingMuon& tag, const Run3ScoutingMuon& probe, Region reg) {
          fillJPack_(cat, ls, reg, c.mass, tag, probe, c.commonVtx);
        };

        doFill(j_all_, mu1, mu2, reg12);
        doFill(j_all_, mu2, mu1, reg21);

        if (jc == JpsiClass::Prompt) {
          doFill(j_prompt_, mu1, mu2, reg12);
          doFill(j_prompt_, mu2, mu1, reg21);
        } else if (jc == JpsiClass::NonPrompt) {
          doFill(j_nonprompt_, mu1, mu2, reg12);
          doFill(j_nonprompt_, mu2, mu1, reg21);
        }
      }
    }

    // J/psi track branch: tag (pt>=probePtMin_) + track probe, match to muon
    // Fill only "All" for now (prompt/nonprompt requires a dimuon vertex).
    {
      // build tag-quality list once (no pt), then require pt>=probePtMin_ in branch
      std::vector<size_t> tagsQ;
      tagsQ.reserve(nMu);
      for (size_t i=0; i<nMu; ++i) {
        const auto& mu = hMu->at(i);
        if (std::abs(mu.eta()) > tagEtaMax_) continue;
        if (!passTrackerQuality_(mu)) continue;
        tagsQ.push_back(i);
      }

      for (size_t it=0; it<tagsQ.size(); ++it) {
        const auto& tag = hMu->at(tagsQ[it]);
        if (!passTagWithPt_(tag, probePtMin_)) continue;

        TLorentzVector vtag, vtrk;
        vtag.SetPtEtaPhiM(tag.pt(), tag.eta(), tag.phi(), MUON_MASS);
        const bool cenTag = isCentral(tag.eta());

        for (const auto& trk : *hTrk) {
          if (!passTrackProbe_(trk)) continue;
          if (trk.tk_charge() == tag.charge()) continue;

          vtrk.SetPtEtaPhiM(trk.tk_pt(), trk.tk_eta(), trk.tk_phi(), MUON_MASS);
          const double mass = (vtag+vtrk).M();
          if (mass < jMassMin_ || mass > jMassMax_) continue;

          const bool cenTrk = isCentral(trk.tk_eta());
          const Region reg = pairRegion(cenTag, cenTrk);

          int match = -1;
          for (size_t j=0; j<nMu; ++j) {
            const auto& mu = hMu->at(j);
            if (reco::deltaR2(mu.eta(), mu.phi(), trk.tk_eta(), trk.tk_phi()) < dr2Match_) { match = (int)j; break; }
          }

          if (match >= 0) {
            fillPerRegion(j_all_.trkMu_pass, ls, reg, mass);
            const bool glok = passGlobalMeasured_(hMu->at((size_t)match));
            if (glok) fillPerRegion(j_all_.trkGlo_pass, ls, reg, mass);
            else      fillPerRegion(j_all_.trkGlo_fail, ls, reg, mass);
          } else {
            fillPerRegion(j_all_.trkMu_fail, ls, reg, mass);
          }
        }
      }
    }
  }
}

void ScoutingZ::endJob() {
  if (!fs_.isAvailable()) return;

  // Determine LS span from all buffers
  unsigned minLS = std::numeric_limits<unsigned>::max();
  unsigned maxLS = 0;

  auto scanAll = [&](const PerRegionBuf& r){ scanLS(r, minLS, maxLS); };
  auto scanZPack = [&](const ZBufPack& b){
    scanAll(b.twoTag); scanAll(b.oneTag);
    scanAll(b.id_pass); scanAll(b.id_fail);
    scanAll(b.glo_pass); scanAll(b.glo_fail);
    scanAll(b.trkMu_pass); scanAll(b.trkMu_fail);
    scanAll(b.trkGlo_pass); scanAll(b.trkGlo_fail);

    scanAll(b.v_twoTag); scanAll(b.v_oneTag);
    scanAll(b.v_id_pass); scanAll(b.v_id_fail);
    scanAll(b.v_glo_pass); scanAll(b.v_glo_fail);
  };
  auto scanJPack = [&](const JpsiCats& b){
    scanAll(b.twoTag); scanAll(b.oneTag);
    scanAll(b.id_pass); scanAll(b.id_fail);
    scanAll(b.glo_pass); scanAll(b.glo_fail);
    scanAll(b.trkMu_pass); scanAll(b.trkMu_fail);
    scanAll(b.trkGlo_pass); scanAll(b.trkGlo_fail);

    scanAll(b.v_twoTag); scanAll(b.v_oneTag);
    scanAll(b.v_id_pass); scanAll(b.v_id_fail);
    scanAll(b.v_glo_pass); scanAll(b.v_glo_fail);
  };

  if (doZ_) {
    scanZPack(zMainBuf_);
    for (const auto& b : zScanBuf_) scanZPack(b);
  }
  if (doJpsi_) {
    scanJPack(j_all_);
    scanJPack(j_prompt_);
    scanJPack(j_nonprompt_);
  }

  if (minLS == std::numeric_limits<unsigned>::max()) return;

  // -----------------------
  // Z booking
  // -----------------------
  if (doZ_) {
    auto zdir = fs_->mkdir("Z");

    // Main (keeps old-ish layout: Z/BB,BE,EE plus Z/Vtx)
    {
      auto mainDir = zdir.mkdir("Main"); // keep explicit to avoid ambiguity with scan
      bookZPack(zMainH2_, mainDir, minLS, maxLS, zMassBins_, zMassMin_, zMassMax_,
                /*withTrackBranch=*/true, /*withVertexed=*/true);

      flushPerRegion(zMainH2_.h2tag, zMainBuf_.twoTag);
      flushPerRegion(zMainH2_.h1tag, zMainBuf_.oneTag);
      flushPerRegion(zMainH2_.hidp,  zMainBuf_.id_pass);
      flushPerRegion(zMainH2_.hidf,  zMainBuf_.id_fail);
      flushPerRegion(zMainH2_.hglop, zMainBuf_.glo_pass);
      flushPerRegion(zMainH2_.hglof, zMainBuf_.glo_fail);

      flushPerRegion(zMainH2_.htrkMu_p, zMainBuf_.trkMu_pass);
      flushPerRegion(zMainH2_.htrkMu_f, zMainBuf_.trkMu_fail);
      flushPerRegion(zMainH2_.htrkGlo_p,zMainBuf_.trkGlo_pass);
      flushPerRegion(zMainH2_.htrkGlo_f,zMainBuf_.trkGlo_fail);

      flushPerRegion(zMainH2_.vh2tag, zMainBuf_.v_twoTag);
      flushPerRegion(zMainH2_.vh1tag, zMainBuf_.v_oneTag);
      flushPerRegion(zMainH2_.vhidp,  zMainBuf_.v_id_pass);
      flushPerRegion(zMainH2_.vhidf,  zMainBuf_.v_id_fail);
      flushPerRegion(zMainH2_.vhglop, zMainBuf_.v_glo_pass);
      flushPerRegion(zMainH2_.vhglof, zMainBuf_.v_glo_fail);
    }

    // Tag pT scan: Z/TagPtScan/ptX/...
    {
      auto scanDir = zdir.mkdir("TagPtScan");
      for (size_t k=0; k<zTagPtScan_.size(); ++k) {
        auto ptDir = scanDir.mkdir(ptLabel(zTagPtScan_[k]));
        bookZPack(zScanH2_[k], ptDir, minLS, maxLS, zMassBins_, zMassMin_, zMassMax_,
                  /*withTrackBranch=*/true, /*withVertexed=*/true);

        const auto& b = zScanBuf_[k];

        flushPerRegion(zScanH2_[k].h2tag, b.twoTag);
        flushPerRegion(zScanH2_[k].h1tag, b.oneTag);
        flushPerRegion(zScanH2_[k].hidp,  b.id_pass);
        flushPerRegion(zScanH2_[k].hidf,  b.id_fail);
        flushPerRegion(zScanH2_[k].hglop, b.glo_pass);
        flushPerRegion(zScanH2_[k].hglof, b.glo_fail);

        flushPerRegion(zScanH2_[k].htrkMu_p, b.trkMu_pass);
        flushPerRegion(zScanH2_[k].htrkMu_f, b.trkMu_fail);
        flushPerRegion(zScanH2_[k].htrkGlo_p,b.trkGlo_pass);
        flushPerRegion(zScanH2_[k].htrkGlo_f,b.trkGlo_fail);

        flushPerRegion(zScanH2_[k].vh2tag, b.v_twoTag);
        flushPerRegion(zScanH2_[k].vh1tag, b.v_oneTag);
        flushPerRegion(zScanH2_[k].vhidp,  b.v_id_pass);
        flushPerRegion(zScanH2_[k].vhidf,  b.v_id_fail);
        flushPerRegion(zScanH2_[k].vhglop, b.v_glo_pass);
        flushPerRegion(zScanH2_[k].vhglof, b.v_glo_fail);
      }
    }
  }

  // -----------------------
  // J/psi booking
  // -----------------------
  if (doJpsi_) {
    auto jdir = fs_->mkdir("Jpsi");

    auto bookClass = [&](const std::string& cls, JH2Pack& out, const JpsiCats& buf, const std::string& title) {
      auto clsDir = jdir.mkdir(cls);
      bookJPack(out, clsDir, title, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_,
                /*withTrackBranch=*/true, /*withVertexed=*/true);

      flushPerRegion(out.h2tag, buf.twoTag);
      flushPerRegion(out.h1tag, buf.oneTag);
      flushPerRegion(out.hidp,  buf.id_pass);
      flushPerRegion(out.hidf,  buf.id_fail);
      flushPerRegion(out.hglop, buf.glo_pass);
      flushPerRegion(out.hglof, buf.glo_fail);

      flushPerRegion(out.htrkMu_p, buf.trkMu_pass);
      flushPerRegion(out.htrkMu_f, buf.trkMu_fail);
      flushPerRegion(out.htrkGlo_p,buf.trkGlo_pass);
      flushPerRegion(out.htrkGlo_f,buf.trkGlo_fail);

      flushPerRegion(out.vh2tag, buf.v_twoTag);
      flushPerRegion(out.vh1tag, buf.v_oneTag);
      flushPerRegion(out.vhidp,  buf.v_id_pass);
      flushPerRegion(out.vhidf,  buf.v_id_fail);
      flushPerRegion(out.vhglop, buf.v_glo_pass);
      flushPerRegion(out.vhglof, buf.v_glo_fail);
    };

    bookClass("All",       jAllH2_, j_all_,       "J/psi All");
    bookClass("Prompt",    jPrH2_,  j_prompt_,    "J/psi Prompt");
    bookClass("NonPrompt", jNpH2_,  j_nonprompt_, "J/psi NonPrompt");
  }
}

DEFINE_FWK_MODULE(ScoutingZ);