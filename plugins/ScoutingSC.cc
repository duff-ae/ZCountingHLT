// ScoutingSC.cc
// Standard candles in Run3 Scouting: Z -> mumu, J/psi -> mumu (all/prompt/nonprompt), W -> mu nu (monitor)
// Per-LS mass-vs-LS histograms for T&P pass/fail categories (DQM-like), plus minimal per-channel Control.
//
// Key design choices:
// - Single Tag WP: "TagTight" (pt>=3, |eta|<2.4 by default), no enforced global via normalizedChi2.
// - Everything else measured via Tag&Probe pass/fail categories.
// - Separate folders: Z/, Jpsi/(All|Prompt|NonPrompt)/, W/ and Z/Control, Jpsi/Control, W/Control.
// - Dynamic LS range: buffer per LS and book 2D at endJob with observed [minLS,maxLS].
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
#include "DataFormats/Math/interface/deltaPhi.h"

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
  template <typename T, typename = void> struct has_trk_chi2_ndof : std::false_type {};
  template <typename T>
  struct has_trk_chi2_ndof<T, std::void_t<
      decltype(std::declval<T>().trk_chi2()),
      decltype(std::declval<T>().trk_ndof())
    >> : std::true_type {};

  template <typename T, typename = void> struct has_trkChi2_trkNdof : std::false_type {};
  template <typename T>
  struct has_trkChi2_trkNdof<T, std::void_t<
      decltype(std::declval<T>().trkChi2()),
      decltype(std::declval<T>().trkNdof())
    >> : std::true_type {};

  template <typename T, typename = void> struct has_trackChi2_trackNdof : std::false_type {};
  template <typename T>
  struct has_trackChi2_trackNdof<T, std::void_t<
      decltype(std::declval<T>().trackChi2()),
      decltype(std::declval<T>().trackNdof())
    >> : std::true_type {};

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
}

class ScoutingSC : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingSC(const edm::ParameterSet&);
  ~ScoutingSC() override = default;

  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

private:
  // ---- Inputs ----
  edm::EDGetTokenT<std::vector<Run3ScoutingMuon>>    muonToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> pvToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingTrack>>  trackToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> muonVtxToken_; // displacedVtx for prompt/nonprompt

  edm::EDGetTokenT<double> metPtToken_;
  edm::EDGetTokenT<double> metPhiToken_;

  bool doZ_{true};
  bool doJpsi_{true};
  bool doW_{true};

  // ---- W: configurable WPs and iso sidebands ----
  std::vector<std::string> wWpNames_{"Loose","Tight"};
  std::vector<double> wWpMuonPtMin_{22.0, 27.0};
  std::vector<double> wWpMetMin_{20.0, 30.0};
  std::vector<double> wWpDphiMin_{1.5, 2.0};
  std::vector<double> wWpMtMin_{40.0, 50.0};
  double wIsoPassMax_{0.15};     // iso-pass upper bound
  double wIsoFailMax_{0.50};     // iso-fail upper bound (QCD sideband)

  // ---- W: denominator (for efficiency-style diagnostics) ----
  double wDenPtMin_{15.0};
  double wDenEtaMax_{2.4};

  // ---- TagTight (universal) ----
  double tagPtMin_{3.0};
  double tagEtaMax_{2.4};
  double tagTrkChi2NdofMax_{3.0}; // POG recommendation (tracker chi2/ndof)
  int    tagMinTrkLayers_{4};
  int    tagMinPixLayers_{0};
  int    tagMinStaStations_{0};   // keep 0 by default (tightness controllable)
  int    tagMinStaHits_{0};

  // ---- Probe base (avoid "anything goes" for 1TAG etc) ----
  double probePtMin_{3.0};
  double probeEtaMax_{2.4};

  // ---- Global / muon-quality component (measured, not required in tag) ----
  double gloNormChi2Max_{999.0};  // only used if isGlobalMuon==true; keep high by default

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

  // W selection (monitor)
  double wPtMin_{27.0};
  double wEtaMax_{2.4};
  double wRelIsoMax_{0.15};
  double wMetMin_{30.0};
  double wDphiMin_{2.0};
  double wMtMax_{160.0};
  double wSecondMuVetoPt_{25.0};

  // PV selection (optional control)
  double vtxNdofMin_{4.0};
  double vtxAbsZMax_{24.0};
  double vtxRhoMax_{2.0};

  // ---- Services ----
  edm::Service<TFileService> fs_;

  // ---- Buffers (per LS) ----
  using Buf = std::map<unsigned, std::vector<double>>;

  struct PerRegionBuf {
    Buf bb, be, ee;
  };

  // ---- Output 2D histos (booked in endJob with observed LS range) ----
  struct PerRegionH2 {
    TH2D* bb{nullptr}; TH2D* be{nullptr}; TH2D* ee{nullptr};
  };

  static void fillPerRegion(PerRegionBuf& r, unsigned ls, Region reg, double val) {
    if (reg == Region::BB) r.bb[ls].push_back(val);
    else if (reg == Region::EE) r.ee[ls].push_back(val);
    else r.be[ls].push_back(val);
  }

  // ---- Flush buffered per-LS values into TH2s ----
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

  // Z: per-LS LS×mass categories (working window)
  PerRegionBuf z_2tag_, z_1tag_;
  PerRegionBuf z_id_pass_, z_id_fail_;
  PerRegionBuf z_sa_pass_, z_sa_fail_;
  PerRegionBuf z_glo_pass_, z_glo_fail_;
  PerRegionBuf z_trkSta_pass_, z_trkSta_fail_;
  PerRegionBuf z_trkGlo_pass_, z_trkGlo_fail_;

  // J/psi: per-LS categories for 3 classes
  struct JpsiCats {
    PerRegionBuf twoTag, oneTag;
    PerRegionBuf id_pass, id_fail;
    PerRegionBuf sa_pass, sa_fail;
    PerRegionBuf glo_pass, glo_fail;
    PerRegionBuf trkSta_pass, trkSta_fail;
    PerRegionBuf trkGlo_pass, trkGlo_fail;
  };
  JpsiCats j_all_, j_prompt_, j_nonprompt_;

  // W: per-LS LS×mt, LS×2pt
  Buf w_mt_, w_m2pt_;

  // PV per-LS (control)
  std::map<unsigned, std::vector<unsigned>> z_npv_, j_npv_, w_npv_;

  // ---- Control histos (1D only) ----
  // Z control
  TH1I* h_z_cutflow_{nullptr};
  TH1D* h_z_mass_ctrl_{nullptr};
  TH1D* h_z_tag_pt_[3]{nullptr,nullptr,nullptr};   // BB/BE/EE by pair region
  TH1D* h_z_tag_eta_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_relIso_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_trkChi2Ndof_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_nTrkLay_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_nPixLay_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_staStations_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_staHits_[3]{nullptr,nullptr,nullptr};
  TH1I* h_z_tag_isGlobal_[3]{nullptr,nullptr,nullptr};
  TH1D* h_z_tag_normChi2_[3]{nullptr,nullptr,nullptr};

  // J/psi control
  TH1I* h_j_cutflow_{nullptr};
  TH1D* h_j_lxy_{nullptr};
  TH1D* h_j_lxySig_{nullptr};

  // W control
  TH1I* h_w_cutflow_{nullptr};
  TH1D* h_w_met_{nullptr};
  TH1D* h_w_mt_{nullptr};
  TH1D* h_w_dphi_{nullptr};
  TH1D* h_w_mu_pt_{nullptr};
  TH1D* h_w_mu_relIso_{nullptr};
  TH1I* h_w_nmu_{nullptr};
  TH1D* h_w_secondMu_pt_{nullptr};

  // Z output hists
  struct ZH2Pack {
    PerRegionH2 h2tag, h1tag;
    PerRegionH2 hidp, hidf;
    PerRegionH2 hsap, hsaf;
    PerRegionH2 hglop, hglof;
    PerRegionH2 htrkStap, htrkStaf;
    PerRegionH2 htrkGlop, htrkGlof;
  } zH2_;

  // J/psi output hists per class
  struct JH2Pack {
    PerRegionH2 h2tag, h1tag;
    PerRegionH2 hidp, hidf;
    PerRegionH2 hsap, hsaf;
    PerRegionH2 hglop, hglof;
    PerRegionH2 htrkStap, htrkStaf;
    PerRegionH2 htrkGlop, htrkGlof;
  } jAllH2_, jPrH2_, jNpH2_;

  // W output
  TH2D* h2_w_mt_{nullptr};
  TH2D* h2_w_m2pt_{nullptr};

  // ---- W: 2D buffers per WP (Signal/QCD) ----
  struct WBufPack {
    PerRegionBuf mt_sig, mt_qcd;
    PerRegionBuf met_sig, met_qcd;
    PerRegionBuf pt_sig,  pt_qcd;
    PerRegionBuf dphi_sig, dphi_qcd;

    // efficiency-like: denom and component-pass pt distributions
    PerRegionBuf pt_denom;     // base denominator
    PerRegionBuf pt_id_pass;   // passTrackerQuality
    PerRegionBuf pt_sa_pass;   // passMuonSystem
    PerRegionBuf pt_glo_pass;  // passGlobalMeasured
    PerRegionBuf pt_iso_pass;  // relIso < wIsoPassMax
  };
  std::vector<WBufPack> wbuf_;  // size = nWP

  // ---- W: output hists per WP ----
  struct WH2Pack {
    PerRegionH2 h_mt_sig, h_mt_qcd;
    PerRegionH2 h_met_sig, h_met_qcd;
    PerRegionH2 h_pt_sig,  h_pt_qcd;
    PerRegionH2 h_dphi_sig, h_dphi_qcd;

    PerRegionH2 h_pt_denom;
    PerRegionH2 h_pt_id_pass;
    PerRegionH2 h_pt_sa_pass;
    PerRegionH2 h_pt_glo_pass;
    PerRegionH2 h_pt_iso_pass;
  };
  std::vector<WH2Pack> wH2_;

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
    // If not available in this CMSSW, we do not cut (keep inclusive) but still can monitor.
    if (chi2ndof >= 0.0 && chi2ndof > tagTrkChi2NdofMax_) return false;
    return true;
  }

  bool passMuonSystem_(const Run3ScoutingMuon& mu) const {
    if (tagMinStaStations_ > 0 && mu.nStandAloneMuonMatchedStations() < tagMinStaStations_) return false;
    if (tagMinStaHits_     > 0 && mu.nValidStandAloneMuonHits() < tagMinStaHits_) return false;
    return true;
  }

  bool passTagTight_(const Run3ScoutingMuon& mu) const {
    if (mu.pt() < tagPtMin_) return false;
    if (std::abs(mu.eta()) > tagEtaMax_) return false;
    if (!passTrackerQuality_(mu)) return false;
    if (!passMuonSystem_(mu)) return false;
    return true;
  }

  bool passGlobalMeasured_(const Run3ScoutingMuon& mu) const {
    if (!mu.isGlobalMuon()) return false;
    // normalizedChi2 is only meaningful for global; keep it as measured component, not tag requirement
    if (mu.normalizedChi2() > gloNormChi2Max_) return false;
    return true;
  }

  bool passTrackProbe_(const Run3ScoutingTrack& trk) const {
    if (trk.tk_pt() < trkPtMin_) return false;
    if (std::abs(trk.tk_eta()) > trkEtaMax_) return false;
    if (trk.tk_nTrackerLayersWithMeasurement() < trkMinLayers_) return false;

    const double chi2ndof = trkChi2OverNdof(trk);
    // if chi2ndof is actually chi2 (no ndof), user can tune trkChi2NdofMax accordingly
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

  TH2D* makeLSVar2D(TFileDirectory& dir,
                    const std::string& name, const std::string& title,
                    unsigned minLS, unsigned maxLS,
                    int ybins, double ymin, double ymax,
                    const std::string& ytitle) const {
    const int nLumiBins = (maxLS >= minLS) ? int(maxLS - minLS + 1) : 1;
    const double lmin = double(minLS) - 0.5;
    const double lmax = double(maxLS) + 0.5;
    auto* h = dir.make<TH2D>(name.c_str(), title.c_str(),
                            nLumiBins, lmin, lmax, ybins, ymin, ymax);
    h->GetXaxis()->SetTitle("luminosity section");
    h->GetYaxis()->SetTitle(ytitle.c_str());
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

  void bookPerRegionH2Var(PerRegionH2& out,
                          TFileDirectory& bbDir, TFileDirectory& beDir, TFileDirectory& eeDir,
                          const std::string& baseName, const std::string& baseTitle,
                          unsigned minLS, unsigned maxLS,
                          int ybins, double ymin, double ymax,
                          const std::string& ytitle) const {
    out.bb = makeLSVar2D(bbDir, baseName, baseTitle + " BB", minLS, maxLS, ybins, ymin, ymax, ytitle);
    out.be = makeLSVar2D(beDir, baseName, baseTitle + " BE", minLS, maxLS, ybins, ymin, ymax, ytitle);
    out.ee = makeLSVar2D(eeDir, baseName, baseTitle + " EE", minLS, maxLS, ybins, ymin, ymax, ytitle);
  }

  // ---- W selection ----
  bool passWTag_(const Run3ScoutingMuon& mu) const {
    if (!passTagTight_(mu)) return false;
    if (mu.pt() < wPtMin_) return false;
    if (std::abs(mu.eta()) > wEtaMax_) return false;
    if (relIso(mu) > (float)wRelIsoMax_) return false;
    return true;
  }

  inline Region muRegion_(const Run3ScoutingMuon& mu) const {
    return isCentral(mu.eta()) ? Region::BB : Region::EE; // BE unused for single-muon
  }

  bool passWDenomMuon_(const Run3ScoutingMuon& mu) const {
    // максимально “безопасный” denominator: только кинематика/геометрия
    if (mu.pt() < wDenPtMin_) return false;
    if (std::abs(mu.eta()) > wDenEtaMax_) return false;
    return true;
  }

  inline bool isIsoPass_(const Run3ScoutingMuon& mu) const {
    return relIso(mu) < (float)wIsoPassMax_;
  }
  inline bool isIsoFail_(const Run3ScoutingMuon& mu) const {
    const float iso = relIso(mu);
    return (iso >= (float)wIsoPassMax_) && (iso < (float)wIsoFailMax_);
  }
};

ScoutingSC::ScoutingSC(const edm::ParameterSet& cfg) {
  usesResource("TFileService");

  muonToken_  = consumes<std::vector<Run3ScoutingMuon>>(   cfg.getParameter<edm::InputTag>("muonCollection") );
  pvToken_    = consumes<std::vector<Run3ScoutingVertex>>( cfg.getParameter<edm::InputTag>("primaryVertexCollection") );
  trackToken_ = consumes<std::vector<Run3ScoutingTrack>>(  cfg.getParameter<edm::InputTag>("trackCollection") );

  if (cfg.existsAs<edm::InputTag>("muonVertexCollection")) {
    muonVtxToken_ = consumes<std::vector<Run3ScoutingVertex>>( cfg.getParameter<edm::InputTag>("muonVertexCollection") );
  }

  doZ_    = cfg.getUntrackedParameter<bool>("doZ", true);
  doJpsi_ = cfg.getUntrackedParameter<bool>("doJpsi", true);
  doW_    = cfg.getUntrackedParameter<bool>("doW", true);

  // Tag/probe
  tagPtMin_ = cfg.getUntrackedParameter<double>("tagPtMin", 3.0);
  tagEtaMax_ = cfg.getUntrackedParameter<double>("tagEtaMax", 2.4);
  tagTrkChi2NdofMax_ = cfg.getUntrackedParameter<double>("tagTrkChi2NdofMax", 3.0);
  tagMinTrkLayers_ = cfg.getUntrackedParameter<int>("tagMinTrackerLayers", 4);
  tagMinPixLayers_ = cfg.getUntrackedParameter<int>("tagMinPixelLayers", 0);
  tagMinStaStations_ = cfg.getUntrackedParameter<int>("tagMinStaStations", 0);
  tagMinStaHits_ = cfg.getUntrackedParameter<int>("tagMinStaHits", 0);

  probePtMin_ = cfg.getUntrackedParameter<double>("probePtMin", 3.0);
  probeEtaMax_ = cfg.getUntrackedParameter<double>("probeEtaMax", 2.4);

  gloNormChi2Max_ = cfg.getUntrackedParameter<double>("gloNormChi2Max", 999.0);

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

  // W
  wPtMin_ = cfg.getUntrackedParameter<double>("wPtMin", 27.0);
  wEtaMax_ = cfg.getUntrackedParameter<double>("wEtaMax", 2.4);
  wRelIsoMax_ = cfg.getUntrackedParameter<double>("wRelIsoMax", 0.15);
  wMetMin_ = cfg.getUntrackedParameter<double>("wMetMin", 30.0);
  wDphiMin_ = cfg.getUntrackedParameter<double>("wDphiMin", 2.0);
  wMtMax_ = cfg.getUntrackedParameter<double>("wMtMax", 160.0);
  wSecondMuVetoPt_ = cfg.getUntrackedParameter<double>("wSecondMuVetoPt", 25.0);

  // W WPs
  wWpNames_      = cfg.getUntrackedParameter<std::vector<std::string>>("wWpNames", wWpNames_);
  wWpMuonPtMin_  = cfg.getUntrackedParameter<std::vector<double>>("wWpMuonPtMin", wWpMuonPtMin_);
  wWpMetMin_     = cfg.getUntrackedParameter<std::vector<double>>("wWpMetMin", wWpMetMin_);
  wWpDphiMin_    = cfg.getUntrackedParameter<std::vector<double>>("wWpDphiMin", wWpDphiMin_);
  wWpMtMin_      = cfg.getUntrackedParameter<std::vector<double>>("wWpMtMin", wWpMtMin_);

  wIsoPassMax_   = cfg.getUntrackedParameter<double>("wIsoPassMax", wIsoPassMax_);
  wIsoFailMax_   = cfg.getUntrackedParameter<double>("wIsoFailMax", wIsoFailMax_);

  wDenPtMin_     = cfg.getUntrackedParameter<double>("wDenPtMin", wDenPtMin_);
  wDenEtaMax_    = cfg.getUntrackedParameter<double>("wDenEtaMax", wDenEtaMax_);

  // sanity: sizes must match
  const size_t nwp = wWpNames_.size();
  if (wWpMuonPtMin_.size() != nwp || wWpMetMin_.size() != nwp ||
      wWpDphiMin_.size() != nwp || wWpMtMin_.size() != nwp) {
    throw cms::Exception("Configuration")
      << "W WP vectors must have same size as wWpNames";
  }

  wbuf_.assign(nwp, WBufPack{});
  wH2_.assign(nwp, WH2Pack{});

  if (doW_) {
    metPtToken_  = consumes<double>( cfg.getUntrackedParameter<edm::InputTag>("pfMetPtValue",  edm::InputTag("hltScoutingPFPacker","pfMetPt")) );
    metPhiToken_ = consumes<double>( cfg.getUntrackedParameter<edm::InputTag>("pfMetPhiValue", edm::InputTag("hltScoutingPFPacker","pfMetPhi")) );
  }
}

void ScoutingSC::beginJob() {
  if (!fs_.isAvailable()) return;

  // Z/Control
  auto zdir = fs_->mkdir("Z");
  auto zc   = zdir.mkdir("Control");

  h_z_cutflow_ = zc.make<TH1I>("cutflow", "Z cutflow", 10, 0.5, 10.5);
  h_z_cutflow_->GetXaxis()->SetBinLabel(1, "events");
  h_z_cutflow_->GetXaxis()->SetBinLabel(2, ">=2 mu");
  h_z_cutflow_->GetXaxis()->SetBinLabel(3, ">=1 tag");
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
    h_z_tag_staStations_[ir] = zc.make<TH1D>(std::string("tag_staStations_"+std::string(rname[ir])).c_str(),
                                             ("tag SA stations "+std::string(rname[ir])+"; n").c_str(), 10, -0.5, 9.5);
    h_z_tag_staHits_[ir] = zc.make<TH1D>(std::string("tag_staHits_"+std::string(rname[ir])).c_str(),
                                         ("tag SA hits "+std::string(rname[ir])+"; n").c_str(), 60, -0.5, 59.5);
    h_z_tag_isGlobal_[ir] = zc.make<TH1I>(std::string("tag_isGlobal_"+std::string(rname[ir])).c_str(),
                                          ("tag isGlobal "+std::string(rname[ir])+"; flag").c_str(), 2, -0.5, 1.5);
    h_z_tag_normChi2_[ir] = zc.make<TH1D>(std::string("tag_normChi2_"+std::string(rname[ir])).c_str(),
                                          ("tag normalizedChi2 "+std::string(rname[ir])+"; #chi^{2}/ndof").c_str(), 200, 0, 200);
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

  // W/Control
  auto wdir = fs_->mkdir("W");
  auto wc   = wdir.mkdir("Control");
  h_w_cutflow_ = wc.make<TH1I>("cutflow", "W cutflow", 10, 0.5, 10.5);
  h_w_cutflow_->GetXaxis()->SetBinLabel(1, "events");
  h_w_cutflow_->GetXaxis()->SetBinLabel(2, "met avail");
  h_w_cutflow_->GetXaxis()->SetBinLabel(3, "exactly1 Wtag");
  h_w_cutflow_->GetXaxis()->SetBinLabel(4, "2nd mu veto");
  h_w_cutflow_->GetXaxis()->SetBinLabel(5, "met cut");
  h_w_cutflow_->GetXaxis()->SetBinLabel(6, "dphi cut");
  h_w_cutflow_->GetXaxis()->SetBinLabel(7, "mt win");

  h_w_met_  = wc.make<TH1D>("met", "MET; MET [GeV]", 200, 0, 200);
  h_w_mt_   = wc.make<TH1D>("mt", "m_{T}; m_{T} [GeV]", 200, 0, 200);
  h_w_dphi_ = wc.make<TH1D>("dphi", "|#Delta#phi(#mu,MET)|; rad", 64, 0, 3.2);

  h_w_mu_pt_     = wc.make<TH1D>("mu_pt", "W mu p_{T}; p_{T} [GeV]", 200, 0, 200);
  h_w_mu_relIso_ = wc.make<TH1D>("mu_relIso", "W mu relIso; relIso", 100, 0, 2);
  h_w_nmu_       = wc.make<TH1I>("nmu", "nMuons; n", 10, -0.5, 9.5);
  h_w_secondMu_pt_= wc.make<TH1D>("secondMu_pt", "2nd mu p_{T} (max); p_{T} [GeV]", 200, 0, 200);

}

void ScoutingSC::analyze(const edm::Event& ev, const edm::EventSetup&) {
  const unsigned ls = ev.luminosityBlock();

  edm::Handle<std::vector<Run3ScoutingMuon>> hMu;
  edm::Handle<std::vector<Run3ScoutingVertex>> hPV;
  edm::Handle<std::vector<Run3ScoutingTrack>> hTrk;
  edm::Handle<std::vector<Run3ScoutingVertex>> hMuVtx;

  ev.getByToken(muonToken_, hMu);
  ev.getByToken(pvToken_, hPV);
  ev.getByToken(trackToken_, hTrk);

  const bool haveMuVtx = !muonVtxToken_.isUninitialized() && ev.getByToken(muonVtxToken_, hMuVtx) && hMuVtx.isValid();

  if (!hMu.isValid() || !hPV.isValid() || !hTrk.isValid()) return;

  const unsigned npv = countGoodPVs_(*hPV);

  if (doZ_)  z_npv_[ls].push_back(npv);
  if (doJpsi_) j_npv_[ls].push_back(npv);
  if (doW_)  w_npv_[ls].push_back(npv);

  // avg PV position for Lxy
  double avgX = 0.0, avgY = 0.0;
  if (!hPV->empty()) {
    avgX = std::accumulate(hPV->begin(), hPV->end(), 0.0, [](double a, const auto& v){return a+v.x();}) / hPV->size();
    avgY = std::accumulate(hPV->begin(), hPV->end(), 0.0, [](double a, const auto& v){return a+v.y();}) / hPV->size();
  }

  const size_t nMu = hMu->size();
  if (doZ_ && h_z_cutflow_) h_z_cutflow_->Fill(1);
  if (doJpsi_ && h_j_cutflow_) h_j_cutflow_->Fill(1);
  if (doW_ && h_w_cutflow_) h_w_cutflow_->Fill(1);

  // =========================
// W: mu nu (store 2D for offline fits + W/Z ratio + efficiency-style diagnostics)
// =========================
if (doW_) {
    edm::Handle<double> hMetPt, hMetPhi;
    const bool okMet = ev.getByToken(metPtToken_, hMetPt) &&
                      ev.getByToken(metPhiToken_, hMetPhi) &&
                      hMetPt.isValid() && hMetPhi.isValid();

    // Для W просто скипаем, но не убиваем Z/Jpsi
    if (okMet) {
      if (h_w_cutflow_) h_w_cutflow_->Fill(2);

      if (h_w_nmu_) h_w_nmu_->Fill((int)nMu);
      if (nMu >= 1) {
        const double metPt  = *hMetPt;
        const double metPhi = *hMetPhi;

        // ---- choose leading denominator muon (loose kin only) ----
        int best = -1;
        double bestPt = -1.0;
        for (size_t i=0; i<nMu; ++i) {
          const auto& mu = hMu->at(i);
          if (!passWDenomMuon_(mu)) continue;
          if (mu.pt() > bestPt) { bestPt = mu.pt(); best = (int)i; }
        }
        if (best >= 0) {
          const auto& mu = hMu->at((size_t)best);
          const Region reg = muRegion_(mu);

          if (h_w_mu_pt_)     h_w_mu_pt_->Fill(mu.pt());
          if (h_w_mu_relIso_) h_w_mu_relIso_->Fill(relIso(mu));

          // ---- “efficiency-style” buffers: denom and component passes vs pt ----
          for (size_t iw=0; iw<wbuf_.size(); ++iw) {
            fillPerRegion(wbuf_[iw].pt_denom,   ls, reg, mu.pt());
            if (passTrackerQuality_(mu)) fillPerRegion(wbuf_[iw].pt_id_pass,  ls, reg, mu.pt());
            if (passMuonSystem_(mu))     fillPerRegion(wbuf_[iw].pt_sa_pass,  ls, reg, mu.pt());
            if (passGlobalMeasured_(mu)) fillPerRegion(wbuf_[iw].pt_glo_pass, ls, reg, mu.pt());
            if (isIsoPass_(mu))          fillPerRegion(wbuf_[iw].pt_iso_pass, ls, reg, mu.pt());
          }

          // ---- second muon veto diagnostic (not hard cut here) ----
          double secondMaxPt = 0.0;
          for (size_t i=0; i<nMu; ++i) {
            if ((int)i == best) continue;
            secondMaxPt = std::max(secondMaxPt, (double)hMu->at(i).pt());
          }
          if (h_w_secondMu_pt_) h_w_secondMu_pt_->Fill(secondMaxPt);

          // ---- compute kinematics for W-like fit observables ----
          const double dphi = std::abs(reco::deltaPhi(mu.phi(), metPhi));
          const double mt   = std::sqrt(2.0 * mu.pt() * metPt * (1.0 - std::cos(dphi)));

          if (h_w_met_)  h_w_met_->Fill(metPt);
          if (h_w_dphi_) h_w_dphi_->Fill(dphi);
          if (h_w_mt_)   h_w_mt_->Fill(mt);

          // ---- apply per-WP selection and fill 2D buffers ----
          const bool isoP = isIsoPass_(mu);
          const bool isoF = isIsoFail_(mu);

          for (size_t iw=0; iw<wbuf_.size(); ++iw) {
            // базовые WP пороги на μ и topology
            if (mu.pt() < wWpMuonPtMin_[iw]) continue;
            if (metPt  < wWpMetMin_[iw])    continue;
            if (dphi   < wWpDphiMin_[iw])   continue;
            if (mt     < wWpMtMin_[iw])     continue;
            if (mt     > wMtMax_)           continue;

            // Дополнительный Z-veto (как у тебя было), но только на WP-уровне:
            bool veto = false;
            for (size_t i=0; i<nMu; ++i) {
              if ((int)i == best) continue;
              const auto& mu2 = hMu->at(i);
              if (mu2.pt() > wSecondMuVetoPt_ && std::abs(mu2.eta()) < wEtaMax_) { veto = true; break; }
            }
            if (veto) continue;

            // Signal (iso-pass) and QCD template (iso-fail)
            if (isoP) {
              fillPerRegion(wbuf_[iw].mt_sig,   ls, reg, mt);
              fillPerRegion(wbuf_[iw].met_sig,  ls, reg, metPt);
              fillPerRegion(wbuf_[iw].pt_sig,   ls, reg, mu.pt());
              fillPerRegion(wbuf_[iw].dphi_sig, ls, reg, dphi);
            } else if (isoF) {
              fillPerRegion(wbuf_[iw].mt_qcd,   ls, reg, mt);
              fillPerRegion(wbuf_[iw].met_qcd,  ls, reg, metPt);
              fillPerRegion(wbuf_[iw].pt_qcd,   ls, reg, mu.pt());
              fillPerRegion(wbuf_[iw].dphi_qcd, ls, reg, dphi);
            }
          }
        }
      }
    }
  }

  if (nMu < 2) return;
  if (doZ_ && h_z_cutflow_) h_z_cutflow_->Fill(2);
  if (doJpsi_ && h_j_cutflow_) h_j_cutflow_->Fill(2);

  // =========================
  // Z: muon-muon T&P, tag=TagTight
  // =========================
  if (doZ_) {
    TLorentzVector v1, v2;

    // collect tag indices
    std::vector<size_t> tags;
    tags.reserve(nMu);
    for (size_t i=0; i<nMu; ++i) if (passTagTight_(hMu->at(i))) tags.push_back(i);
    if (!tags.empty() && h_z_cutflow_) h_z_cutflow_->Fill(3);

    for (size_t it=0; it<tags.size(); ++it) {
      const auto& tag = hMu->at(tags[it]);
      v1.SetPtEtaPhiM(tag.pt(), tag.eta(), tag.phi(), MUON_MASS);
      const bool cenTag = isCentral(tag.eta());

      for (size_t j=0; j<nMu; ++j) {
        if (j == tags[it]) continue;
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
        const int ir = (reg==Region::BB ? 0 : (reg==Region::BE ? 1 : 2));

        // control fill for tag (per pair region)
        if (h_z_tag_pt_[ir])  h_z_tag_pt_[ir]->Fill(tag.pt());
        if (h_z_tag_eta_[ir]) h_z_tag_eta_[ir]->Fill(tag.eta());
        if (h_z_tag_relIso_[ir]) h_z_tag_relIso_[ir]->Fill(relIso(tag));
        const double tchi = muTrkChi2OverNdof(tag);
        if (h_z_tag_trkChi2Ndof_[ir]) h_z_tag_trkChi2Ndof_[ir]->Fill(tchi >= 0.0 ? tchi : -0.1);
        if (h_z_tag_nTrkLay_[ir]) h_z_tag_nTrkLay_[ir]->Fill(tag.nTrackerLayersWithMeasurement());
        if (h_z_tag_nPixLay_[ir]) h_z_tag_nPixLay_[ir]->Fill(tag.nPixelLayersWithMeasurement());
        if (h_z_tag_staStations_[ir]) h_z_tag_staStations_[ir]->Fill(tag.nStandAloneMuonMatchedStations());
        if (h_z_tag_staHits_[ir]) h_z_tag_staHits_[ir]->Fill(tag.nValidStandAloneMuonHits());
        if (h_z_tag_isGlobal_[ir]) h_z_tag_isGlobal_[ir]->Fill(tag.isGlobalMuon() ? 1 : 0);
        if (h_z_tag_normChi2_[ir]) h_z_tag_normChi2_[ir]->Fill(tag.normalizedChi2());

        // 2TAG / 1TAG (probe tight measured)
        const bool probeTight = passTagTight_(probe);
        if (probeTight) fillPerRegion(z_2tag_, ls, reg, mass);
        else           fillPerRegion(z_1tag_, ls, reg, mass);

        // component efficiencies measured on probe
        const bool idok = passTrackerQuality_(probe);
        if (idok) fillPerRegion(z_id_pass_, ls, reg, mass);
        else      fillPerRegion(z_id_fail_, ls, reg, mass);

        const bool saok = passMuonSystem_(probe);
        if (saok) fillPerRegion(z_sa_pass_, ls, reg, mass);
        else      fillPerRegion(z_sa_fail_, ls, reg, mass);

        const bool glok = passGlobalMeasured_(probe);
        if (glok) fillPerRegion(z_glo_pass_, ls, reg, mass);
        else      fillPerRegion(z_glo_fail_, ls, reg, mass);
      }
    }

    // =========================
    // Z: muon-track branch (tag tight + track probe), mass in Z window
    // =========================
    for (size_t it=0; it<tags.size(); ++it) {
      const auto& tag = hMu->at(tags[it]);
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
          fillPerRegion(z_trkSta_pass_, ls, reg, mass);
          const bool glok = passGlobalMeasured_(hMu->at((size_t)match));
          if (glok) fillPerRegion(z_trkGlo_pass_, ls, reg, mass);
          else      fillPerRegion(z_trkGlo_fail_, ls, reg, mass);
        } else {
          fillPerRegion(z_trkSta_fail_, ls, reg, mass);
        }
      }
    }
  }

  // =========================
  // J/psi: all/prompt/nonprompt, tag tight, same categories
  // =========================
  if (doJpsi_) {
    TLorentzVector v1, v2;

    const bool canClassify = haveMuVtx && hMuVtx.isValid() && !hMuVtx->empty() && !hPV->empty();

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

        JpsiClass jc = JpsiClass::Unknown;
        if (canClassify) jc = classifyJpsi_(mu1, mu2, *hMuVtx, avgX, avgY);

        // allow either muon to be tag (standard T&P)
        auto fillPair = [&](const Run3ScoutingMuon& tag, const Run3ScoutingMuon& probe,
                            bool cenTag, bool cenProbe) {
          if (!passTagTight_(tag)) return;

          const Region reg = pairRegion(cenTag, cenProbe);

          auto doFill = [&](JpsiCats& cat) {
            const bool probeTight = passTagTight_(probe);
            if (probeTight) fillPerRegion(cat.twoTag, ls, reg, mass);
            else           fillPerRegion(cat.oneTag, ls, reg, mass);

            const bool idok = passTrackerQuality_(probe);
            if (idok) fillPerRegion(cat.id_pass, ls, reg, mass);
            else      fillPerRegion(cat.id_fail, ls, reg, mass);

            const bool saok = passMuonSystem_(probe);
            if (saok) fillPerRegion(cat.sa_pass, ls, reg, mass);
            else      fillPerRegion(cat.sa_fail, ls, reg, mass);

            const bool glok = passGlobalMeasured_(probe);
            if (glok) fillPerRegion(cat.glo_pass, ls, reg, mass);
            else      fillPerRegion(cat.glo_fail, ls, reg, mass);
          };

          doFill(j_all_);
          if (jc == JpsiClass::Prompt) doFill(j_prompt_);
          if (jc == JpsiClass::NonPrompt) doFill(j_nonprompt_);
        };

        const bool cen1 = isCentral(mu1.eta());
        const bool cen2 = isCentral(mu2.eta());
        fillPair(mu1, mu2, cen1, cen2);
        fillPair(mu2, mu1, cen2, cen1);
      }
    }

    // J/psi track branch: tag tight + track probe, match to muon
    // (kept symmetric with Z for later comparisons)
    // Build tag list once
    std::vector<size_t> tags;
    tags.reserve(nMu);
    for (size_t i=0; i<nMu; ++i) if (passTagTight_(hMu->at(i))) tags.push_back(i);

    for (size_t it=0; it<tags.size(); ++it) {
      const auto& tag = hMu->at(tags[it]);
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

        auto doFillTrack = [&](JpsiCats& cat) {
          if (match >= 0) {
            fillPerRegion(cat.trkSta_pass, ls, reg, mass);
            const bool glok = passGlobalMeasured_(hMu->at((size_t)match));
            if (glok) fillPerRegion(cat.trkGlo_pass, ls, reg, mass);
            else      fillPerRegion(cat.trkGlo_fail, ls, reg, mass);
          } else {
            fillPerRegion(cat.trkSta_fail, ls, reg, mass);
          }
        };

        doFillTrack(j_all_);
        // prompt/nonprompt for track branch needs per-event J/psi vertex; we do not have it here (tag+trk),
        // so keep only "all" for now. (This is intentional to avoid fake classification.)
      }
    }
  }
}

void ScoutingSC::endJob() {
  if (!fs_.isAvailable()) return;

  // Determine LS span from all buffers
  unsigned minLS = std::numeric_limits<unsigned>::max();
  unsigned maxLS = 0;

  auto scanAll = [&](const PerRegionBuf& r){ scanLS(r, minLS, maxLS); };
  auto scanW   = [&](const Buf& b){ scanLS(b, minLS, maxLS); };

  if (doZ_) {
    scanAll(z_2tag_); scanAll(z_1tag_);
    scanAll(z_id_pass_); scanAll(z_id_fail_);
    scanAll(z_sa_pass_); scanAll(z_sa_fail_);
    scanAll(z_glo_pass_); scanAll(z_glo_fail_);
    scanAll(z_trkSta_pass_); scanAll(z_trkSta_fail_);
    scanAll(z_trkGlo_pass_); scanAll(z_trkGlo_fail_);
  }
  if (doJpsi_) {
    auto scanJ = [&](const JpsiCats& c){
      scanAll(c.twoTag); scanAll(c.oneTag);
      scanAll(c.id_pass); scanAll(c.id_fail);
      scanAll(c.sa_pass); scanAll(c.sa_fail);
      scanAll(c.glo_pass); scanAll(c.glo_fail);
      scanAll(c.trkSta_pass); scanAll(c.trkSta_fail);
      scanAll(c.trkGlo_pass); scanAll(c.trkGlo_fail);
    };
    scanJ(j_all_);
    scanJ(j_prompt_);
    scanJ(j_nonprompt_);
  }
  if (doW_) {
    scanW(w_mt_);
    scanW(w_m2pt_);
    for (const auto& wb : wbuf_) {
      scanAll(wb.mt_sig);  scanAll(wb.mt_qcd);
      scanAll(wb.met_sig); scanAll(wb.met_qcd);
      scanAll(wb.pt_sig);  scanAll(wb.pt_qcd);
      scanAll(wb.dphi_sig);scanAll(wb.dphi_qcd);

      scanAll(wb.pt_denom);
      scanAll(wb.pt_id_pass);
      scanAll(wb.pt_sa_pass);
      scanAll(wb.pt_glo_pass);
      scanAll(wb.pt_iso_pass);
    }
  }

  if (minLS == std::numeric_limits<unsigned>::max()) return;

  // Book output folders and per-LS 2D hists with observed LS range
  // Z
  if (doZ_) {
    auto zdir = fs_->mkdir("Z");
    auto zBB  = zdir.mkdir("BB");
    auto zBE  = zdir.mkdir("BE");
    auto zEE  = zdir.mkdir("EE");

    bookPerRegionH2(zH2_.h2tag, zBB, zBE, zEE, "h_mass_2TAG",      "Z 2TAG",      minLS, maxLS, zMassBins_, zMassMin_, zMassMax_);
    bookPerRegionH2(zH2_.h1tag, zBB, zBE, zEE, "h_mass_1TAG",      "Z 1TAG",      minLS, maxLS, zMassBins_, zMassMin_, zMassMax_);

    bookPerRegionH2(zH2_.hidp,  zBB, zBE, zEE, "h_mass_ID_pass",   "Z ID pass",   minLS, maxLS, zMassBins_, zMassMin_, zMassMax_);
    bookPerRegionH2(zH2_.hidf,  zBB, zBE, zEE, "h_mass_ID_fail",   "Z ID fail",   minLS, maxLS, zMassBins_, zMassMin_, zMassMax_);

    bookPerRegionH2(zH2_.hsap,  zBB, zBE, zEE, "h_mass_SA_pass",   "Z SA pass",   minLS, maxLS, zMassBins_, zMassMin_, zMassMax_);
    bookPerRegionH2(zH2_.hsaf,  zBB, zBE, zEE, "h_mass_SA_fail",   "Z SA fail",   minLS, maxLS, zMassBins_, zMassMin_, zMassMax_);

    bookPerRegionH2(zH2_.hglop, zBB, zBE, zEE, "h_mass_Glo_pass",  "Z Global pass", minLS, maxLS, zMassBins_, zMassMin_, zMassMax_);
    bookPerRegionH2(zH2_.hglof, zBB, zBE, zEE, "h_mass_Glo_fail",  "Z Global fail", minLS, maxLS, zMassBins_, zMassMin_, zMassMax_);

    bookPerRegionH2(zH2_.htrkStap, zBB, zBE, zEE, "h_mass_TrkSta_pass", "Z TrkSta pass", minLS, maxLS, zMassBins_, zMassMin_, zMassMax_);
    bookPerRegionH2(zH2_.htrkStaf, zBB, zBE, zEE, "h_mass_TrkSta_fail", "Z TrkSta fail", minLS, maxLS, zMassBins_, zMassMin_, zMassMax_);
    bookPerRegionH2(zH2_.htrkGlop, zBB, zBE, zEE, "h_mass_TrkGlo_pass", "Z TrkGlo pass", minLS, maxLS, zMassBins_, zMassMin_, zMassMax_);
    bookPerRegionH2(zH2_.htrkGlof, zBB, zBE, zEE, "h_mass_TrkGlo_fail", "Z TrkGlo fail", minLS, maxLS, zMassBins_, zMassMin_, zMassMax_);

    flushPerRegion(zH2_.h2tag, z_2tag_);
    flushPerRegion(zH2_.h1tag, z_1tag_);
    flushPerRegion(zH2_.hidp, z_id_pass_);
    flushPerRegion(zH2_.hidf, z_id_fail_);
    flushPerRegion(zH2_.hsap, z_sa_pass_);
    flushPerRegion(zH2_.hsaf, z_sa_fail_);
    flushPerRegion(zH2_.hglop, z_glo_pass_);
    flushPerRegion(zH2_.hglof, z_glo_fail_);
    flushPerRegion(zH2_.htrkStap, z_trkSta_pass_);
    flushPerRegion(zH2_.htrkStaf, z_trkSta_fail_);
    flushPerRegion(zH2_.htrkGlop, z_trkGlo_pass_);
    flushPerRegion(zH2_.htrkGlof, z_trkGlo_fail_);
  }

  // J/psi: book per class (All/Prompt/NonPrompt)
  if (doJpsi_) {
    auto jdir = fs_->mkdir("Jpsi");

    auto bookJ = [&](const std::string& cls, JH2Pack& out, const JpsiCats& buf) {
      auto clsDir = jdir.mkdir(cls);
      auto jBB = clsDir.mkdir("BB");
      auto jBE = clsDir.mkdir("BE");
      auto jEE = clsDir.mkdir("EE");

      bookPerRegionH2(out.h2tag, jBB, jBE, jEE, "h_mass_2TAG", "J/psi 2TAG " + cls, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_);
      bookPerRegionH2(out.h1tag, jBB, jBE, jEE, "h_mass_1TAG", "J/psi 1TAG " + cls, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_);
      bookPerRegionH2(out.hidp,  jBB, jBE, jEE, "h_mass_ID_pass", "J/psi ID pass " + cls, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_);
      bookPerRegionH2(out.hidf,  jBB, jBE, jEE, "h_mass_ID_fail", "J/psi ID fail " + cls, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_);
      bookPerRegionH2(out.hsap,  jBB, jBE, jEE, "h_mass_SA_pass", "J/psi SA pass " + cls, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_);
      bookPerRegionH2(out.hsaf,  jBB, jBE, jEE, "h_mass_SA_fail", "J/psi SA fail " + cls, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_);
      bookPerRegionH2(out.hglop, jBB, jBE, jEE, "h_mass_Glo_pass","J/psi Glo pass " + cls, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_);
      bookPerRegionH2(out.hglof, jBB, jBE, jEE, "h_mass_Glo_fail","J/psi Glo fail " + cls, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_);

      bookPerRegionH2(out.htrkStap, jBB, jBE, jEE, "h_mass_TrkSta_pass", "J/psi TrkSta pass " + cls, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_);
      bookPerRegionH2(out.htrkStaf, jBB, jBE, jEE, "h_mass_TrkSta_fail", "J/psi TrkSta fail " + cls, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_);
      bookPerRegionH2(out.htrkGlop, jBB, jBE, jEE, "h_mass_TrkGlo_pass", "J/psi TrkGlo pass " + cls, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_);
      bookPerRegionH2(out.htrkGlof, jBB, jBE, jEE, "h_mass_TrkGlo_fail", "J/psi TrkGlo fail " + cls, minLS, maxLS, jMassBins_, jMassMin_, jMassMax_);

      flushPerRegion(out.h2tag, buf.twoTag);
      flushPerRegion(out.h1tag, buf.oneTag);
      flushPerRegion(out.hidp,  buf.id_pass);
      flushPerRegion(out.hidf,  buf.id_fail);
      flushPerRegion(out.hsap,  buf.sa_pass);
      flushPerRegion(out.hsaf,  buf.sa_fail);
      flushPerRegion(out.hglop, buf.glo_pass);
      flushPerRegion(out.hglof, buf.glo_fail);
      flushPerRegion(out.htrkStap, buf.trkSta_pass);
      flushPerRegion(out.htrkStaf, buf.trkSta_fail);
      flushPerRegion(out.htrkGlop, buf.trkGlo_pass);
      flushPerRegion(out.htrkGlof, buf.trkGlo_fail);
    };

    bookJ("All",       jAllH2_, j_all_);
    bookJ("Prompt",    jPrH2_,  j_prompt_);
    bookJ("NonPrompt", jNpH2_,  j_nonprompt_);
  }

  // W per-LS 2D
  if (doW_) {
    auto wdir = fs_->mkdir("W");

    for (size_t iw=0; iw<wWpNames_.size(); ++iw) {
      auto wpDir = wdir.mkdir(wWpNames_[iw]);

      auto sigDir = wpDir.mkdir("Signal");
      auto qcdDir = wpDir.mkdir("QCD");
      auto effDir = wpDir.mkdir("Eff"); // efficiency-style diagnostics

      auto mkRegions = [&](TFileDirectory& base){
        auto bb = base.mkdir("BB");
        auto be = base.mkdir("BE");
        auto ee = base.mkdir("EE");
        return std::tuple<TFileDirectory,TFileDirectory,TFileDirectory>(bb,be,ee);
      };

      // --- Signal ---
      {
        auto [bb,be,ee] = mkRegions(sigDir);

        bookPerRegionH2Var(wH2_[iw].h_mt_sig,  bb,be,ee, "h_mt",   "W m_{T} (iso-pass)", minLS, maxLS, 160, 0.0, 160.0, "m_{T} [GeV]");
        bookPerRegionH2Var(wH2_[iw].h_met_sig, bb,be,ee, "h_met",  "W MET (iso-pass)",   minLS, maxLS, 200, 0.0, 200.0,  "MET [GeV]");
        bookPerRegionH2Var(wH2_[iw].h_pt_sig,  bb,be,ee, "h_pt",   "W mu p_{T} (iso-pass)", minLS, maxLS, 200, 0.0, 200.0, "p_{T} [GeV]");
        bookPerRegionH2Var(wH2_[iw].h_dphi_sig,bb,be,ee, "h_dphi", "|#Delta#phi(#mu,MET)| (iso-pass)", minLS, maxLS, 64, 0.0, 3.2, "rad");

        flushPerRegion(wH2_[iw].h_mt_sig,   wbuf_[iw].mt_sig);
        flushPerRegion(wH2_[iw].h_met_sig,  wbuf_[iw].met_sig);
        flushPerRegion(wH2_[iw].h_pt_sig,   wbuf_[iw].pt_sig);
        flushPerRegion(wH2_[iw].h_dphi_sig, wbuf_[iw].dphi_sig);
      }

      // --- QCD (iso-fail) ---
      {
        auto [bb,be,ee] = mkRegions(qcdDir);

        bookPerRegionH2Var(wH2_[iw].h_mt_qcd,  bb,be,ee, "h_mt",   "W m_{T} (iso-fail)", minLS, maxLS, 160, 0.0, 160.0, "m_{T} [GeV]");
        bookPerRegionH2Var(wH2_[iw].h_met_qcd, bb,be,ee, "h_met",  "W MET (iso-fail)",   minLS, maxLS, 200, 0.0, 200.0,  "MET [GeV]");
        bookPerRegionH2Var(wH2_[iw].h_pt_qcd,  bb,be,ee, "h_pt",   "W mu p_{T} (iso-fail)", minLS, maxLS, 200, 0.0, 200.0, "p_{T} [GeV]");
        bookPerRegionH2Var(wH2_[iw].h_dphi_qcd,bb,be,ee, "h_dphi", "|#Delta#phi(#mu,MET)| (iso-fail)", minLS, maxLS, 64, 0.0, 3.2, "rad");

        flushPerRegion(wH2_[iw].h_mt_qcd,   wbuf_[iw].mt_qcd);
        flushPerRegion(wH2_[iw].h_met_qcd,  wbuf_[iw].met_qcd);
        flushPerRegion(wH2_[iw].h_pt_qcd,   wbuf_[iw].pt_qcd);
        flushPerRegion(wH2_[iw].h_dphi_qcd, wbuf_[iw].dphi_qcd);
      }

      // --- Eff (denom and component-pass) ---
      {
        auto [bb,be,ee] = mkRegions(effDir);

        bookPerRegionH2Var(wH2_[iw].h_pt_denom,   bb,be,ee, "h_pt_denom", "W denom mu p_{T}", minLS, maxLS, 200, 0.0, 200.0, "p_{T} [GeV]");
        bookPerRegionH2Var(wH2_[iw].h_pt_id_pass, bb,be,ee, "h_pt_IDpass","W ID pass p_{T}",  minLS, maxLS, 200, 0.0, 200.0, "p_{T} [GeV]");
        bookPerRegionH2Var(wH2_[iw].h_pt_sa_pass, bb,be,ee, "h_pt_SApass","W SA pass p_{T}",  minLS, maxLS, 200, 0.0, 200.0, "p_{T} [GeV]");
        bookPerRegionH2Var(wH2_[iw].h_pt_glo_pass,bb,be,ee, "h_pt_Glopass","W Glo pass p_{T}",minLS, maxLS, 200, 0.0, 200.0, "p_{T} [GeV]");
        bookPerRegionH2Var(wH2_[iw].h_pt_iso_pass,bb,be,ee, "h_pt_IsoPass","W Iso pass p_{T}",minLS, maxLS, 200, 0.0, 200.0, "p_{T} [GeV]");

        flushPerRegion(wH2_[iw].h_pt_denom,    wbuf_[iw].pt_denom);
        flushPerRegion(wH2_[iw].h_pt_id_pass,  wbuf_[iw].pt_id_pass);
        flushPerRegion(wH2_[iw].h_pt_sa_pass,  wbuf_[iw].pt_sa_pass);
        flushPerRegion(wH2_[iw].h_pt_glo_pass, wbuf_[iw].pt_glo_pass);
        flushPerRegion(wH2_[iw].h_pt_iso_pass, wbuf_[iw].pt_iso_pass);
      }
    }
  }
}


DEFINE_FWK_MODULE(ScoutingSC);
