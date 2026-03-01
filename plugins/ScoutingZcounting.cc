// ScoutingZCounting.cc
// Extended Z + J/psi (all/prompt/nonprompt) + W
// Tag-and-probe style, scouting-friendly cuts; tighten from cfg later.

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
#include "TH2D.h"

#include <map>
#include <unordered_set>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cctype>

namespace {
  constexpr double MUON_MASS   = 0.1056583755;
  constexpr double MUON_BOUND  = 0.9;      // |eta|<0.9 → barrel
  constexpr double DRMAX_IO2   = 0.09;     // ΔR^2 ~ 0.3 for matching track↔mu

  // IsoMu24-ish kinematics for Z tag
  constexpr double ISOMU24_PT_MIN  = 24.0;
  constexpr double ISOMU24_ETA_MAX = 2.4;

  // Soft quality defaults for Scouting
  constexpr int    MIN_TRK_LAYERS_DEF  = 4;
  constexpr int    MIN_PIX_LAYERS_DEF  = 0;
  constexpr double MAX_GLO_CHI2_DEF    = 20.0;

  // J/psi window & prompt threshold
  constexpr double JPSI_MMIN = 2.6;
  constexpr double JPSI_MMAX = 3.6;
  constexpr double LXY_SIG_PROMPT_MAX = 3.0;

  // ---- W selection (single-muon topology) ----
  constexpr double W_PT_MIN          = 27.0;   // pT(mu) threshold for W tag
  constexpr double W_ETA_MAX         = 2.4;    // |eta(mu)| < 2.4
  constexpr double W_RELISO_MAX      = 0.15;   // tighter relIso for W
  constexpr double W_MET_MIN         = 30.0;   // MET threshold
  constexpr double W_MT_MIN          = 0.0;   // mT window
  constexpr double W_MT_MAX          = 160.0;
  constexpr double W_DPHI_MIN        = 2.0;    // Δφ(mu,MET) > 2 rad
  constexpr double W_SECOND_MU_PT_VETO = 25.0; // veto on 2nd high-pT muon

  inline std::string sanitize(const std::string& s) {
    std::string r; r.reserve(s.size());
    for (char c : s) r.push_back(std::isalnum((unsigned char)c) ? c : '_');
    return r;
  }

  inline bool isCentral(double eta) {
    return std::abs(eta) < MUON_BOUND;
  }

  inline float muRelIso(const Run3ScoutingMuon& mu) {
    const float absIso = mu.trackIso() + mu.ecalIso() + mu.hcalIso();
    return absIso / std::max(1.f, (float)mu.pt());
  }

  enum class JpsiClass {
    Unknown = 0,
    Prompt  = 1,
    NonPrompt = 2
  };

  // Indices for J/psi T&P arrays
  enum JpsiTPType {
    kHLT2 = 0,
    kHLT1 = 1,
    kStaPass = 2,
    kStaFail = 3,
    kGloPass = 4,
    kGloFail = 5,
    kIDFail  = 6,
    kNTPTypes = 7
  };

  enum JpsiCat {
    kAllCat = 0,
    kPromptCat = 1,
    kNonPromptCat = 2,
    kNCats = 3
  };

  enum Region {
    kBB = 0,
    kBE = 1,
    kEE = 2,
    kNRegions = 3
  };
}

class ScoutingZCounting : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingZCounting(const edm::ParameterSet&);
  ~ScoutingZCounting() override = default;

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;

private:
  // ==== Inputs ====
  edm::EDGetTokenT<std::vector<Run3ScoutingMuon>>      muonToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>>    pvToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingTrack>>     trackToken_;
  // muon vertices for J/psi prompt/non-prompt
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>>    muonVtxToken_;

  // pfMET (optional, for W)
  edm::EDGetTokenT<double>                             metPtToken_;
  edm::EDGetTokenT<double>                             metPhiToken_;
  bool                                                 doW_{false};

  // ==== Core config (Z) ====
  double PtCutL2_;
  double EtaCutL2_;
  int    MassBin_;
  double MassMin_;
  double MassMax_;
  int    PVBin_;
  double PVMin_;
  double PVMax_;

  // PV cuts
  double VtxNTracksFitCut_;
  double VtxNdofCut_;
  double VtxAbsZCut_;
  double VtxRhoCut_;

  // Proxies parameters
  double      isoCutHLT_;           // relIso cut for IsoMu24 proxy
  int         minTrkLayers_;
  int         minPixLayers_;
  double      maxChi2_;

  // switches
  bool        doJpsi_{true};
  bool        doZ_{true};

  // === Services / histos ===
  edm::Service<TFileService> fs_;

  // Z 2D: LS × mass (original Z-counting style)
  TH2D *h_mass_2HLT_BB_{nullptr}, *h_mass_2HLT_BE_{nullptr}, *h_mass_2HLT_EE_{nullptr};
  TH2D *h_mass_1HLT_BB_{nullptr}, *h_mass_1HLT_BE_{nullptr}, *h_mass_1HLT_EE_{nullptr};
  TH2D *h_mass_Sta_pass_BB_{nullptr}, *h_mass_Sta_pass_BE_{nullptr}, *h_mass_Sta_pass_EE_{nullptr};
  TH2D *h_mass_Sta_fail_BB_{nullptr}, *h_mass_Sta_fail_BE_{nullptr}, *h_mass_Sta_fail_EE_{nullptr};
  TH2D *h_mass_Glo_pass_BB_{nullptr}, *h_mass_Glo_pass_BE_{nullptr}, *h_mass_Glo_pass_EE_{nullptr};
  TH2D *h_mass_Glo_fail_BB_{nullptr}, *h_mass_Glo_fail_BE_{nullptr}, *h_mass_Glo_fail_EE_{nullptr};
  TH2D *h_mass_ID_fail_BB_{nullptr},  *h_mass_ID_fail_BE_{nullptr},  *h_mass_ID_fail_EE_{nullptr};

  // PV
  TH2D *h_npv_{nullptr};

  // J/psi T&P style: 3 категории (all/prompt/nonprompt) × 7 типов T&P × 3 региона (BB/BE/EE)
  TH2D* h_jpsi_[kNCats][kNTPTypes][kNRegions]{{{nullptr}}};

  // J/psi 1D monitor (global)
  TH1D *h_lxy_{nullptr}, *h_lxySig_{nullptr};

  // W 2D: LS × mT (без разделения по заряду)
  TH2D *h_w_mt_{nullptr};
  TH2D *h_w_m2pt_{nullptr};

  // buffers per-LS (Z)
  std::map<unsigned, std::vector<double>> buf_mass_2HLT_BB_,    buf_mass_2HLT_BE_,    buf_mass_2HLT_EE_;
  std::map<unsigned, std::vector<double>> buf_mass_1HLT_BB_,    buf_mass_1HLT_BE_,    buf_mass_1HLT_EE_;
  std::map<unsigned, std::vector<double>> buf_mass_Sta_pass_BB_, buf_mass_Sta_pass_BE_, buf_mass_Sta_pass_EE_;
  std::map<unsigned, std::vector<double>> buf_mass_Sta_fail_BB_, buf_mass_Sta_fail_BE_, buf_mass_Sta_fail_EE_;
  std::map<unsigned, std::vector<double>> buf_mass_Glo_pass_BB_, buf_mass_Glo_pass_BE_, buf_mass_Glo_pass_EE_;
  std::map<unsigned, std::vector<double>> buf_mass_Glo_fail_BB_, buf_mass_Glo_fail_BE_, buf_mass_Glo_fail_EE_;
  std::map<unsigned, std::vector<double>> buf_mass_ID_fail_BB_,  buf_mass_ID_fail_BE_,  buf_mass_ID_fail_EE_;
  std::map<unsigned, std::vector<unsigned>> buf_npv_;

  // J/psi (3 категории × 7 T&P типов × 3 региона) — те же структуры, что и для Z, но отдельно
  std::map<unsigned, std::vector<double>> buf_jpsi_[kNCats][kNTPTypes][kNRegions];

  // W
  std::map<unsigned, std::vector<double>> buf_w_mt_;
  std::map<unsigned, std::vector<double>> buf_w_m2pt_;

  // ==== Helpers ====
  bool passPV_(const Run3ScoutingVertex& v) const;
  unsigned selectGoodPVs_(const std::vector<Run3ScoutingVertex>& pvs, int& firstGoodIdx) const;

  // Proxies (Scouting-friendly)
  bool passIsoMu24Proxy_(const Run3ScoutingMuon& mu) const;   // Z tag (IsoMu24-like)
  bool passGlobalMuonProxy_(const Run3ScoutingMuon& mu) const;
  bool passMuonIDProxy_(const Run3ScoutingMuon& mu) const;    // ID w/o isolation
  static bool passTrack_(const Run3ScoutingTrack& trk);

  // J/psi Tag definition (softer, low-pT)
  bool passJpsiTagProxy_(const Run3ScoutingMuon& mu) const;

  // J/psi prompt/non-prompt классификатор (возвращает класс и заполняет Lxy-гисты)
  JpsiClass classifyJpsiPromptNonPrompt_(const Run3ScoutingMuon& mu1,
                                         const Run3ScoutingMuon& mu2,
                                         const std::vector<Run3ScoutingVertex>& muonVertices,
                                         double avgX, double avgY) const;

  // Fill BB/BE/EE helper
  static inline void fillBBBE_(std::map<unsigned, std::vector<double>>& bb,
                               std::map<unsigned, std::vector<double>>& be,
                               std::map<unsigned, std::vector<double>>& ee,
                               unsigned ls, bool firstCentral, bool secondCentral, double mass) {
    if (firstCentral && secondCentral)        bb[ls].push_back(mass);
    else if (!firstCentral && !secondCentral) ee[ls].push_back(mass);
    else                                      be[ls].push_back(mass);
  }

  bool isWTagMuon_(const Run3ScoutingMuon& mu) const;

  void bookAllHistosWithLS_(unsigned minLS, unsigned maxLS);

  // J/psi helper to fill T&P arrays
  void fillJpsiTP_(JpsiCat cat, JpsiTPType tp, unsigned ls,
                   bool tagCentral, bool probeCentral, double mass);

  // misc
  int  runNumber_{-1};
  bool saveIntermediate_{true};

  // debug counters
  mutable uint64_t nEvt_=0, nMu_=0, nTrk_=0, nTagHLT_=0, nPairs_=0, n2HLT_=0, n1HLT_=0;
  mutable uint64_t nStaPass_=0, nStaFail_=0, nGloPass_=0, nGloFail_=0;
};

// --------------------------------------------------------------------------------------------

ScoutingZCounting::ScoutingZCounting(const edm::ParameterSet& iCfg) {
  usesResource("TFileService");

  // Core inputs
  muonToken_  = consumes<std::vector<Run3ScoutingMuon>>(   iCfg.getParameter<edm::InputTag>("muonCollection") );
  pvToken_    = consumes<std::vector<Run3ScoutingVertex>>( iCfg.getParameter<edm::InputTag>("primaryVertexCollection") );
  trackToken_ = consumes<std::vector<Run3ScoutingTrack>>(  iCfg.getParameter<edm::InputTag>("trackCollection") );

  // Muon vertices (для J/psi prompt/non-prompt)
  if (iCfg.existsAs<edm::InputTag>("muonVertexCollection")) {
    muonVtxToken_ = consumes<std::vector<Run3ScoutingVertex>>( iCfg.getParameter<edm::InputTag>("muonVertexCollection") );
  }

  // pfMET for W
  doW_ = iCfg.getUntrackedParameter<bool>("doW", true);
  if (doW_) {
    metPtToken_  = consumes<double>( iCfg.getUntrackedParameter<edm::InputTag>("pfMetPtValue",  edm::InputTag("hltScoutingPFPacker","pfMetPt")) );
    metPhiToken_ = consumes<double>( iCfg.getUntrackedParameter<edm::InputTag>("pfMetPhiValue", edm::InputTag("hltScoutingPFPacker","pfMetPhi")) );
  }

  // Z kinematics/binning
  PtCutL2_  = iCfg.getUntrackedParameter<double>("PtCutL2", 15.0);
  EtaCutL2_ = iCfg.getUntrackedParameter<double>("EtaCutL2", 2.4);
  MassBin_  = iCfg.getUntrackedParameter<int>("MassBin", 60);
  MassMin_  = iCfg.getUntrackedParameter<double>("MassMin", 60.0);
  MassMax_  = iCfg.getUntrackedParameter<double>("MassMax", 120.0);

  // PV histo binning
  PVBin_ = iCfg.getUntrackedParameter<int>("PVBin", 100);
  PVMin_ = iCfg.getUntrackedParameter<double>("PVMin", 0.);
  PVMax_ = iCfg.getUntrackedParameter<double>("PVMax", 100.);

  // PV cuts
  VtxNTracksFitCut_ = iCfg.getUntrackedParameter<double>("VtxNTracksFitMin", 0.);
  VtxNdofCut_       = iCfg.getUntrackedParameter<double>("VtxNdofMin", 4.);
  VtxAbsZCut_       = iCfg.getUntrackedParameter<double>("VtxAbsZMax", 24.);
  VtxRhoCut_        = iCfg.getUntrackedParameter<double>("VtxRhoMax", 2.);

  // switches
  doZ_    = iCfg.getUntrackedParameter<bool>("doZ", true);
  doJpsi_ = iCfg.getUntrackedParameter<bool>("doJpsi", true);

  // Misc
  runNumber_        = iCfg.getUntrackedParameter<int>("runNumber", -1);
  saveIntermediate_ = iCfg.getUntrackedParameter<bool>("saveIntermediate", true);

  // proxies tunables
  isoCutHLT_        = iCfg.getUntrackedParameter<double>("hltIsoRelCut", 0.25); // relIso cut

  minTrkLayers_     = iCfg.getUntrackedParameter<int>("minTrackerLayers", MIN_TRK_LAYERS_DEF);
  minPixLayers_     = iCfg.getUntrackedParameter<int>("minPixelLayers",   MIN_PIX_LAYERS_DEF);
  maxChi2_          = iCfg.getUntrackedParameter<double>("maxGlobalChi2", MAX_GLO_CHI2_DEF);
}

// ------------------------------ PV helpers ---------------------------------------------------

bool ScoutingZCounting::passPV_(const Run3ScoutingVertex& v) const {
  if (!v.isValidVtx()) return false;
  if (v.tracksSize() < VtxNTracksFitCut_) return false;
  if (v.ndof() < VtxNdofCut_) return false;
  if (std::abs(v.z()) > VtxAbsZCut_) return false;
  const double rho = std::hypot(v.x(), v.y());
  if (rho > VtxRhoCut_) return false;
  return true;
}

unsigned ScoutingZCounting::selectGoodPVs_(const std::vector<Run3ScoutingVertex>& pvs, int& firstGoodIdx) const {
  unsigned nvtx = 0; firstGoodIdx = -1;
  for (size_t i = 0; i < pvs.size(); ++i) {
    if (passPV_(pvs[i])) {
      if (firstGoodIdx < 0) firstGoodIdx = (int)i;
      ++nvtx;
    }
  }
  return nvtx;
}

// ------------------------------ Proxies ------------------------------------------------------

bool ScoutingZCounting::passIsoMu24Proxy_(const Run3ScoutingMuon& mu) const {
  const float pt  = mu.pt();
  const float eta = mu.eta();

  if (pt < ISOMU24_PT_MIN)    return false;
  if (std::abs(eta) > ISOMU24_ETA_MAX) return false;
  if (!mu.isGlobalMuon())     return false;
  if (mu.normalizedChi2() >= maxChi2_) return false;

  if (mu.nTrackerLayersWithMeasurement() < minTrkLayers_) return false;
  if (mu.nPixelLayersWithMeasurement()   < minPixLayers_)  return false;

  const float relIso = muRelIso(mu);
  if (relIso >= (float)isoCutHLT_) return false;

  return true;
}

bool ScoutingZCounting::passGlobalMuonProxy_(const Run3ScoutingMuon& mu) const {
  if (!mu.isGlobalMuon()) return false;
  if (mu.normalizedChi2() >= maxChi2_) return false;
  // оставим без жёстких требований к слоям, чтобы отделить ID от чисто global-флага
  return true;
}

bool ScoutingZCounting::passMuonIDProxy_(const Run3ScoutingMuon& mu) const {
  // ID без изоляции: глобальный мюон, разумная χ² и слои трекера
  if (!mu.isGlobalMuon()) return false;
  if (mu.normalizedChi2() >= maxChi2_) return false;
  if (mu.nTrackerLayersWithMeasurement() < minTrkLayers_) return false;
  if (mu.nPixelLayersWithMeasurement()   < minPixLayers_)  return false;
  return true;
}

bool ScoutingZCounting::passTrack_(const Run3ScoutingTrack& trk) {
  if (trk.tk_pt() < 3.f) return false;
  if (std::abs(trk.tk_eta()) > ISOMU24_ETA_MAX) return false;
  if (trk.tk_chi2() >= 10.0f) return false;
  if (trk.tk_nTrackerLayersWithMeasurement() < 4) return false;
  return true;
}

// Softer tag for J/psi T&P (scouting-friendly)
bool ScoutingZCounting::passJpsiTagProxy_(const Run3ScoutingMuon& mu) const {
  if (!mu.isGlobalMuon()) return false;
  if (mu.pt() < 6.0) return false;
  if (std::abs(mu.eta()) > 2.4) return false;

  if (!passMuonIDProxy_(mu)) return false;

  const float relIso = muRelIso(mu);
  if (relIso >= 0.30f) return false;

  return true;
}

// ------------------------------ W tag helper ---------------------------------

bool ScoutingZCounting::isWTagMuon_(const Run3ScoutingMuon& mu) const {
  // Базовое качество: global proxy (как для Z)
  if (!passGlobalMuonProxy_(mu)) return false;
  if (std::abs(mu.eta()) > W_ETA_MAX) return false;
  if (mu.pt() < W_PT_MIN) return false;

  // Более строгая изоляция, чем для Z-tag
  const float absIso = mu.trackIso() + mu.ecalIso() + mu.hcalIso();
  const float relIso = absIso / std::max(1.f, (float)mu.pt());
  if (relIso >= (float)W_RELISO_MAX) return false;

  return true;
}

// ------------------------------ J/psi prompt/non-prompt helper ------------------------------

JpsiClass ScoutingZCounting::classifyJpsiPromptNonPrompt_(
  const Run3ScoutingMuon& mu1,
  const Run3ScoutingMuon& mu2,
  const std::vector<Run3ScoutingVertex>& muonVertices,
  double avgX, double avgY) const
{
  const auto& idx1 = mu1.vtxIndx();
  const auto& idx2 = mu2.vtxIndx();
  if (idx1.empty() || idx2.empty()) return JpsiClass::Unknown;

  std::unordered_set<int> set1(idx1.begin(), idx1.end());
  for (int vidx : idx2) {
    if (vidx < 0 || static_cast<size_t>(vidx) >= muonVertices.size()) continue;
    if (!set1.count(vidx)) continue;

    const auto& vtx = muonVertices[static_cast<size_t>(vidx)];

    // dx, dy относительно среднего PV (ровно как в рабочем коде-ntuplizer)
    const double dx = vtx.x() - avgX;
    const double dy = vtx.y() - avgY;
    const double r2 = dx*dx + dy*dy;
    if (r2 <= 0.0) break;

    const double Lxy = std::sqrt(r2);

    // Ошибка только из muon-vertex
    const double sx2 = double(vtx.xError()) * double(vtx.xError());
    const double sy2 = double(vtx.yError()) * double(vtx.yError());
    const double num = dx*dx*sx2 + dy*dy*sy2;
    const double LxyErr = (num > 0.0) ? std::sqrt(num) / Lxy : -1.0;
    const double LxySig = (LxyErr > 0.0) ? (Lxy / LxyErr) : -1.0;

    if (h_lxy_)    h_lxy_->Fill(Lxy);
    if (h_lxySig_) h_lxySig_->Fill(LxySig);

    if (LxyErr <= 0.0) return JpsiClass::Unknown;

    if (LxySig < LXY_SIG_PROMPT_MAX) return JpsiClass::Prompt;
    else                             return JpsiClass::NonPrompt;
  }
  return JpsiClass::Unknown;
}

// ------------------------------ beginJob -----------------------------------------------------

void ScoutingZCounting::beginJob() {
  // 1D-мониторы для Lxy и LxySig бронируем сразу,
  // чтобы classifyJpsiPromptNonPrompt_ мог их заполнять во время analyze().
  if (fs_.isAvailable()) {
    if (!h_lxy_) {
      h_lxy_ = fs_->make<TH1D>("h_lxy", "L_{xy}", 1000, -1.0, 1.0);
    }
    if (!h_lxySig_) {
      h_lxySig_ = fs_->make<TH1D>("h_lxySig", "L_{xy}/#sigma", 1000, -10.0, 100.0);
    }
  }
}

// ------------------------------ J/psi T&P helper --------------------------------------------

void ScoutingZCounting::fillJpsiTP_(JpsiCat cat, JpsiTPType tp, unsigned ls,
                                    bool tagCentral, bool probeCentral, double mass) {
  auto& bb = buf_jpsi_[cat][tp][kBB];
  auto& be = buf_jpsi_[cat][tp][kBE];
  auto& ee = buf_jpsi_[cat][tp][kEE];
  fillBBBE_(bb, be, ee, ls, tagCentral, probeCentral, mass);
}

// ------------------------------ Analyze ------------------------------------------------------

void ScoutingZCounting::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  ++nEvt_;
  const unsigned ls = iEvent.luminosityBlock();

  // Inputs
  edm::Handle<std::vector<Run3ScoutingMuon>>   hMu;
  edm::Handle<std::vector<Run3ScoutingVertex>> hPV;
  edm::Handle<std::vector<Run3ScoutingTrack>>  hTrk;
  edm::Handle<std::vector<Run3ScoutingVertex>> hMuVtx;

  iEvent.getByToken(muonToken_,  hMu);
  iEvent.getByToken(pvToken_,    hPV);
  iEvent.getByToken(trackToken_, hTrk);

  if (!muonVtxToken_.isUninitialized()) {
    iEvent.getByToken(muonVtxToken_, hMuVtx);
  }

  if (!hMu.isValid() || !hPV.isValid() || !hTrk.isValid()) return;
  if (hMu->empty()) return;
  nMu_  += hMu->size();
  nTrk_ += hTrk->size();

  // PV per-event
  int firstPV = -1;
  const unsigned nvtx = selectGoodPVs_(*hPV, firstPV);
  buf_npv_[ls].push_back(nvtx);

  // PV for pseudodecay Lxy
  double avgX = 0.0, avgY = 0.0;
  if (!hPV->empty()) {
    avgX = std::accumulate(hPV->begin(), hPV->end(), 0.0,
            [](double acc, const auto& v){ return acc + v.x(); }) / hPV->size();
    avgY = std::accumulate(hPV->begin(), hPV->end(), 0.0,
            [](double acc, const auto& v){ return acc + v.y(); }) / hPV->size();
  }

  // ---------------- Z: tag & probe on muons ----------------
  if (doZ_) {
    TLorentzVector vi, vj;
    const size_t n = hMu->size();
    for (size_t i = 0; i < n; ++i) {
      const auto& mu_i = hMu->at(i);
      vi.SetPtEtaPhiM(mu_i.pt(), mu_i.eta(), mu_i.phi(), MUON_MASS);
      const bool tagPassHLT = passIsoMu24Proxy_(mu_i);
      if (!tagPassHLT) continue;
      ++nTagHLT_;

      const bool cen_i  = isCentral(mu_i.eta());
      const int  qTag   = mu_i.charge();

      for (size_t j = i + 1; j < n; ++j) {
        const auto& mu_j = hMu->at(j);
        if (qTag == mu_j.charge()) continue;
        if (mu_j.pt() < PtCutL2_ || std::abs(mu_j.eta()) > EtaCutL2_) continue;

        vj.SetPtEtaPhiM(mu_j.pt(), mu_j.eta(), mu_j.phi(), MUON_MASS);
        const double mass = (vi + vj).M();
        if (mass < MassMin_ || mass > MassMax_) continue;
        ++nPairs_;

        const bool cen_j  = isCentral(mu_j.eta());
        const bool probePassHLT = passIsoMu24Proxy_(mu_j);

        if (probePassHLT) {
          ++n2HLT_;
          fillBBBE_(buf_mass_2HLT_BB_, buf_mass_2HLT_BE_, buf_mass_2HLT_EE_, ls, cen_i, cen_j, mass);
        } else {
          ++n1HLT_;
          fillBBBE_(buf_mass_1HLT_BB_, buf_mass_1HLT_BE_, buf_mass_1HLT_EE_, ls, cen_i, cen_j, mass);
        }

        const bool probeIsGlobal = passGlobalMuonProxy_(mu_j);
        const bool probeIDok     = passMuonIDProxy_(mu_j);
        if (probeIsGlobal && !probeIDok) {
          fillBBBE_(buf_mass_ID_fail_BB_, buf_mass_ID_fail_BE_, buf_mass_ID_fail_EE_, ls, cen_i, cen_j, mass);
        }
      }
    }
  }

  // ---------------- Track probes for standalone/global (Z T&P: mu tag + track probe) --------
  {
    TLorentzVector vtag, vtrk;
    const size_t nMu = hMu->size();
    const size_t nTk = hTrk->size();

    for (size_t i = 0; i < nMu; ++i) {
      const auto& mu_i = hMu->at(i);
      if (!passIsoMu24Proxy_(mu_i)) continue;
      vtag.SetPtEtaPhiM(mu_i.pt(), mu_i.eta(), mu_i.phi(), MUON_MASS);
      const bool cen_i = isCentral(mu_i.eta());
      const int  qTag  = mu_i.charge();

      for (size_t k = 0; k < nTk; ++k) {
        const auto& trk = hTrk->at(k);
        if (qTag == trk.tk_charge()) continue;
        if (!passTrack_(trk)) continue;

        vtrk.SetPtEtaPhiM(trk.tk_pt(), trk.tk_eta(), trk.tk_phi(), MUON_MASS);
        const double mass = (vtag + vtrk).M();
        if (mass < MassMin_ || mass > MassMax_) continue;

        const bool cen_k = isCentral(trk.tk_eta());

        int muMatch = -1;
        for (size_t j = 0; j < nMu; ++j) {
          const auto& mu_j = hMu->at(j);
          if (reco::deltaR2(mu_j.eta(), mu_j.phi(), trk.tk_eta(), trk.tk_phi()) < DRMAX_IO2) {
            muMatch = (int)j; break;
          }
        }

        if (muMatch >= 0) {
          ++nStaPass_;
          fillBBBE_(buf_mass_Sta_pass_BB_, buf_mass_Sta_pass_BE_, buf_mass_Sta_pass_EE_, ls, cen_i, cen_k, mass);

          const bool gloOK = passGlobalMuonProxy_(hMu->at(muMatch));
          if (gloOK) {
            ++nGloPass_;
            fillBBBE_(buf_mass_Glo_pass_BB_, buf_mass_Glo_pass_BE_, buf_mass_Glo_pass_EE_, ls, cen_i, cen_k, mass);
          } else {
            ++nGloFail_;
            fillBBBE_(buf_mass_Glo_fail_BB_, buf_mass_Glo_fail_BE_, buf_mass_Glo_fail_EE_, ls, cen_i, cen_k, mass);
          }
        } else {
          ++nStaFail_;
          fillBBBE_(buf_mass_Sta_fail_BB_, buf_mass_Sta_fail_BE_, buf_mass_Sta_fail_EE_, ls, cen_i, cen_k, mass);
        }
      }
    }
  }

  // ---------------- J/psi counting (all/prompt/nonprompt, Z-like T&P) -----------------------
  if (doJpsi_) {
    TLorentzVector vi, vj;
    const size_t n = hMu->size();
    const bool haveMuVtx = hMuVtx.isValid() && !hMuVtx->empty() && !hPV->empty();

    for (size_t i = 0; i < n; ++i) {
      const auto& mu_i = hMu->at(i);
      vi.SetPtEtaPhiM(mu_i.pt(), mu_i.eta(), mu_i.phi(), MUON_MASS);

      for (size_t j = i + 1; j < n; ++j) {
        const auto& mu_j = hMu->at(j);
        vj.SetPtEtaPhiM(mu_j.pt(), mu_j.eta(), mu_j.phi(), MUON_MASS);

        const double mass = (vi + vj).M();
        if (mass < JPSI_MMIN || mass > JPSI_MMAX) continue;
        if (mu_i.charge() == mu_j.charge()) continue; // OS only

        // prompt/non-prompt классификатор
        JpsiClass jclass = JpsiClass::Unknown;
        if (haveMuVtx) {
          jclass = classifyJpsiPromptNonPrompt_(mu_i, mu_j, *hMuVtx, avgX, avgY);
        }

        const bool cen_i = isCentral(mu_i.eta());
        const bool cen_j = isCentral(mu_j.eta());

        auto fillAllCats = [&](const Run3ScoutingMuon& tag, const Run3ScoutingMuon& probe,
                               bool cenTag, bool cenProbe) {
          if (!passJpsiTagProxy_(tag)) return;

          const bool probeHLT = passJpsiTagProxy_(probe);    // J/psi "HLT-like" proxy
          const bool probeGlo = passGlobalMuonProxy_(probe);
          const bool probeID  = passMuonIDProxy_(probe);

          // Standalone-подобный флаг: есть ли hits+matched stations в мюонной системе
          const bool probeSta = (probe.nValidStandAloneMuonHits() > 0 &&
                                 probe.nStandAloneMuonMatchedStations() > 0);

          const double m = mass;

          // Все события
          if (probeHLT) fillJpsiTP_(kAllCat,    kHLT2,    ls, cenTag, cenProbe, m);
          else          fillJpsiTP_(kAllCat,    kHLT1,    ls, cenTag, cenProbe, m);

          if (probeSta) fillJpsiTP_(kAllCat,    kStaPass, ls, cenTag, cenProbe, m);
          else          fillJpsiTP_(kAllCat,    kStaFail, ls, cenTag, cenProbe, m);

          if (probeGlo) fillJpsiTP_(kAllCat,    kGloPass, ls, cenTag, cenProbe, m);
          else          fillJpsiTP_(kAllCat,    kGloFail, ls, cenTag, cenProbe, m);

          if (probeGlo && !probeID)
                        fillJpsiTP_(kAllCat,    kIDFail,  ls, cenTag, cenProbe, m);

          // prompt
          if (jclass == JpsiClass::Prompt) {
            if (probeHLT) fillJpsiTP_(kPromptCat, kHLT2,    ls, cenTag, cenProbe, m);
            else          fillJpsiTP_(kPromptCat, kHLT1,    ls, cenTag, cenProbe, m);

            if (probeSta) fillJpsiTP_(kPromptCat, kStaPass, ls, cenTag, cenProbe, m);
            else          fillJpsiTP_(kPromptCat, kStaFail, ls, cenTag, cenProbe, m);

            if (probeGlo) fillJpsiTP_(kPromptCat, kGloPass, ls, cenTag, cenProbe, m);
            else          fillJpsiTP_(kPromptCat, kGloFail, ls, cenTag, cenProbe, m);

            if (probeGlo && !probeID)
                           fillJpsiTP_(kPromptCat, kIDFail,  ls, cenTag, cenProbe, m);
          }

          // non-prompt
          if (jclass == JpsiClass::NonPrompt) {
            if (probeHLT) fillJpsiTP_(kNonPromptCat, kHLT2,    ls, cenTag, cenProbe, m);
            else          fillJpsiTP_(kNonPromptCat, kHLT1,    ls, cenTag, cenProbe, m);

            if (probeSta) fillJpsiTP_(kNonPromptCat, kStaPass, ls, cenTag, cenProbe, m);
            else          fillJpsiTP_(kNonPromptCat, kStaFail, ls, cenTag, cenProbe, m);

            if (probeGlo) fillJpsiTP_(kNonPromptCat, kGloPass, ls, cenTag, cenProbe, m);
            else          fillJpsiTP_(kNonPromptCat, kGloFail, ls, cenTag, cenProbe, m);

            if (probeGlo && !probeID)
                           fillJpsiTP_(kNonPromptCat, kIDFail,  ls, cenTag, cenProbe, m);
          }
        };

        // Разрешаем любой из мюонов быть tag, как в стандартном T&P
        fillAllCats(mu_i, mu_j, cen_i, cen_j);
        fillAllCats(mu_j, mu_i, cen_j, cen_i);
      }
    }
  }

  // ---------------- W→μν (single-muon topology) ----------------
  if (doW_) {
    edm::Handle<double> hMetPt, hMetPhi;
    iEvent.getByToken(metPtToken_,  hMetPt);
    iEvent.getByToken(metPhiToken_, hMetPhi);
    if (!hMetPt.isValid() || !hMetPhi.isValid()) {
      // Нет MET — пропускаем W-анализ
    } else {
      const double metPt  = *hMetPt;
      const double metPhi = *hMetPhi;

      const size_t nMu = hMu->size();

      // Найти W-tag мюоны
      std::vector<size_t> wTagIdx;
      wTagIdx.reserve(nMu);
      for (size_t i = 0; i < nMu; ++i) {
        const auto& mu = hMu->at(i);
        if (isWTagMuon_(mu)) wTagIdx.push_back(i);
      }

      // Требуем ровно один хороший мюон-кандидат (single-muon topology)
      if (wTagIdx.size() != 1) {
        // либо нет W-tag, либо их несколько -> фон (Z/ttbar/QCD)
      } else {
        const size_t idxTag = wTagIdx.front();
        const auto&  muTag  = hMu->at(idxTag);

        // Veto: второй "тяжёлый" мюон (подавление Z→μμ, tt̄, T&P)
        bool hasSecondHighPtMuon = false;
        for (size_t i = 0; i < nMu; ++i) {
          if (i == idxTag) continue;
          const auto& mu = hMu->at(i);
          if (mu.pt() > W_SECOND_MU_PT_VETO && std::abs(mu.eta()) < W_ETA_MAX) {
            hasSecondHighPtMuon = true;
            break;
          }
        }
        if (!hasSecondHighPtMuon) {
          // MET cut
          if (metPt >= W_MET_MIN) {
            const double dphi = std::abs(reco::deltaPhi(muTag.phi(), metPhi));
            // Требуем "не коллинеарности" MET и мюона (отбрасываем плохой MET)
            if (dphi >= W_DPHI_MIN) {
              const double mt = std::sqrt(2.0 * muTag.pt() * metPt * (1.0 - std::cos(dphi)));
              if (mt >= W_MT_MIN && mt <= W_MT_MAX) {
                // Записываем две величины: настоящий mT и proxy 2*pT(mu)
                buf_w_mt_[ls].push_back(mt);
                buf_w_m2pt_[ls].push_back(2.0 * muTag.pt());
              }
            }
          }
        }
      }
    }
  }
}

// ------------------------------ Book & Flush --------------------------------------------------

void ScoutingZCounting::bookAllHistosWithLS_(unsigned minLS, unsigned maxLS) {
  const int    nLumiBins = (maxLS >= minLS) ? int(maxLS - minLS + 1) : 1;
  const double lmin      = std::max(0, (int)minLS) - 0.5;
  const double lmax      = (int)maxLS + 0.5;

  auto make2D = [&](const char* name, const char* title, int ybins, double ymin, double ymax) -> TH2D* {
    return fs_->make<TH2D>(name, title, nLumiBins, lmin, lmax, ybins, ymin, ymax);
  };
  auto setAxes = [&](TH2D* h, const char* ytitle) {
    if (!h) return;
    h->GetXaxis()->SetTitle("luminosity section");
    h->GetYaxis()->SetTitle(ytitle);
  };

  // Z mass-vs-LS
  h_mass_2HLT_BB_ = make2D("h_mass_2HLT_BB", "2HLT BB", MassBin_, MassMin_, MassMax_);
  h_mass_2HLT_BE_ = make2D("h_mass_2HLT_BE", "2HLT BE", MassBin_, MassMin_, MassMax_);
  h_mass_2HLT_EE_ = make2D("h_mass_2HLT_EE", "2HLT EE", MassBin_, MassMin_, MassMax_);
  h_mass_1HLT_BB_ = make2D("h_mass_1HLT_BB", "1HLT BB", MassBin_, MassMin_, MassMax_);
  h_mass_1HLT_BE_ = make2D("h_mass_1HLT_BE", "1HLT BE", MassBin_, MassMin_, MassMax_);
  h_mass_1HLT_EE_ = make2D("h_mass_1HLT_EE", "1HLT EE", MassBin_, MassMin_, MassMax_);
  h_mass_Sta_pass_BB_ = make2D("h_mass_Sta_pass_BB", "Sta pass BB", MassBin_, MassMin_, MassMax_);
  h_mass_Sta_pass_BE_ = make2D("h_mass_Sta_pass_BE", "Sta pass BE", MassBin_, MassMin_, MassMax_);
  h_mass_Sta_pass_EE_ = make2D("h_mass_Sta_pass_EE", "Sta pass EE", MassBin_, MassMin_, MassMax_);
  h_mass_Sta_fail_BB_ = make2D("h_mass_Sta_fail_BB", "Sta fail BB", MassBin_, MassMin_, MassMax_);
  h_mass_Sta_fail_BE_ = make2D("h_mass_Sta_fail_BE", "Sta fail BE", MassBin_, MassMin_, MassMax_);
  h_mass_Sta_fail_EE_ = make2D("h_mass_Sta_fail_EE", "Sta fail EE", MassBin_, MassMin_, MassMax_);
  h_mass_Glo_pass_BB_ = make2D("h_mass_Glo_pass_BB", "Glo pass BB", MassBin_, MassMin_, MassMax_);
  h_mass_Glo_pass_BE_ = make2D("h_mass_Glo_pass_BE", "Glo pass BE", MassBin_, MassMin_, MassMax_);
  h_mass_Glo_pass_EE_ = make2D("h_mass_Glo_pass_EE", "Glo pass EE", MassBin_, MassMin_, MassMax_);
  h_mass_Glo_fail_BB_ = make2D("h_mass_Glo_fail_BB", "Glo fail BB", MassBin_, MassMin_, MassMax_);
  h_mass_Glo_fail_BE_ = make2D("h_mass_Glo_fail_BE", "Glo fail BE", MassBin_, MassMin_, MassMax_);
  h_mass_Glo_fail_EE_ = make2D("h_mass_Glo_fail_EE", "Glo fail EE", MassBin_, MassMin_, MassMax_);
  h_mass_ID_fail_BB_  = make2D("h_mass_ID_fail_BB",  "ID fail BB",  MassBin_, MassMin_, MassMax_);
  h_mass_ID_fail_BE_  = make2D("h_mass_ID_fail_BE",  "ID fail BE",  MassBin_, MassMin_, MassMax_);
  h_mass_ID_fail_EE_  = make2D("h_mass_ID_fail_EE",  "ID fail EE",  MassBin_, MassMin_, MassMax_);

  for (auto* h : {h_mass_2HLT_BB_,h_mass_2HLT_BE_,h_mass_2HLT_EE_,
                  h_mass_1HLT_BB_,h_mass_1HLT_BE_,h_mass_1HLT_EE_,
                  h_mass_Sta_pass_BB_,h_mass_Sta_pass_BE_,h_mass_Sta_pass_EE_,
                  h_mass_Sta_fail_BB_,h_mass_Sta_fail_BE_,h_mass_Sta_fail_EE_,
                  h_mass_Glo_pass_BB_,h_mass_Glo_pass_BE_,h_mass_Glo_pass_EE_,
                  h_mass_Glo_fail_BB_,h_mass_Glo_fail_BE_,h_mass_Glo_fail_EE_,
                  h_mass_ID_fail_BB_, h_mass_ID_fail_BE_, h_mass_ID_fail_EE_})
    setAxes(h, "tag+probe mass");

  // PV
  h_npv_ = fs_->make<TH2D>("h_npv", "valid PV per event", nLumiBins, lmin, lmax, PVBin_, PVMin_, PVMax_);
  setAxes(h_npv_, "number of primary vertices");

  // J/psi T&P: 3 категории × 7 типов × 3 региона
  const char* catNames[kNCats]     = {"all", "prompt", "nonprompt"};
  const char* catTitles[kNCats]    = {"all", "prompt", "non-prompt"};
  const char* tpNames[kNTPTypes]   = {"2HLT", "1HLT", "Sta_pass", "Sta_fail", "Glo_pass", "Glo_fail", "ID_fail"};
  const char* tpTitles[kNTPTypes]  = {"2HLT", "1HLT", "Sta pass", "Sta fail", "Glo pass", "Glo fail", "ID fail"};
  const char* regNames[kNRegions]  = {"BB", "BE", "EE"};
  const char* regTitles[kNRegions] = {"BB", "BE", "EE"};

  for (int ic = 0; ic < kNCats; ++ic) {
    for (int it = 0; it < kNTPTypes; ++it) {
      for (int ir = 0; ir < kNRegions; ++ir) {
        std::string name  = std::string("h_jpsi_") + catNames[ic] + "_" + tpNames[it] + "_" + regNames[ir];
        std::string title = std::string("J/psi ") + catTitles[ic] + " " +
                            tpTitles[it] + " " + regTitles[ir];
        h_jpsi_[ic][it][ir] = fs_->make<TH2D>(name.c_str(), title.c_str(),
                                              nLumiBins, lmin, lmax,
                                              50, JPSI_MMIN, JPSI_MMAX);
        setAxes(h_jpsi_[ic][it][ir], "dimuon mass");
      }
    }
  }

  // W (mT vs LS)
  h_w_mt_ = fs_->make<TH2D>("h_w_mt", "W m_{T}", nLumiBins, lmin, lmax, 80, 0.0, 160.0);
  setAxes(h_w_mt_,  "m_{T} [GeV]");
  h_w_m2pt_ = fs_->make<TH2D>("h_w_m2pt", "W m2pt proxy",
                            nLumiBins, lmin, lmax,
                            80, 0.0, 160.0);
  setAxes(h_w_m2pt_, "m_{2pt} = 2·p_{T}^{μ} [GeV]");
}

void ScoutingZCounting::endJob() {
  // Determine LS range from all buffers
  unsigned minLS = std::numeric_limits<unsigned>::max();
  unsigned maxLS = 0;
  auto scan = [&](const auto& m) {
    for (const auto& kv : m) {
      minLS = std::min(minLS, kv.first);
      maxLS = std::max(maxLS, kv.first);
    }
  };

  // Z
  scan(buf_mass_2HLT_BB_);    scan(buf_mass_2HLT_BE_);    scan(buf_mass_2HLT_EE_);
  scan(buf_mass_1HLT_BB_);    scan(buf_mass_1HLT_BE_);    scan(buf_mass_1HLT_EE_);
  scan(buf_mass_Sta_pass_BB_); scan(buf_mass_Sta_pass_BE_); scan(buf_mass_Sta_pass_EE_);
  scan(buf_mass_Sta_fail_BB_); scan(buf_mass_Sta_fail_BE_); scan(buf_mass_Sta_fail_EE_);
  scan(buf_mass_Glo_pass_BB_); scan(buf_mass_Glo_pass_BE_); scan(buf_mass_Glo_pass_EE_);
  scan(buf_mass_Glo_fail_BB_); scan(buf_mass_Glo_fail_BE_); scan(buf_mass_Glo_fail_EE_);
  scan(buf_mass_ID_fail_BB_);  scan(buf_mass_ID_fail_BE_);  scan(buf_mass_ID_fail_EE_);
  scan(buf_npv_);
  // J/psi
  for (int ic = 0; ic < kNCats; ++ic)
    for (int it = 0; it < kNTPTypes; ++it)
      for (int ir = 0; ir < kNRegions; ++ir)
        scan(buf_jpsi_[ic][it][ir]);
  // W
  scan(buf_w_mt_);

  if (minLS == std::numeric_limits<unsigned>::max()) return; // no data

  // Book with exact LS span
  bookAllHistosWithLS_(minLS, maxLS);

  auto fill2D = [&](TH2D* h, const std::map<unsigned, std::vector<double>>& m) {
    if (!h) return;
    for (const auto& [ls, vals] : m)
      for (double v : vals)
        h->Fill(ls, v);
  };

  // Z flush
  fill2D(h_mass_2HLT_BB_, buf_mass_2HLT_BB_);
  fill2D(h_mass_2HLT_BE_, buf_mass_2HLT_BE_);
  fill2D(h_mass_2HLT_EE_, buf_mass_2HLT_EE_);
  fill2D(h_mass_1HLT_BB_, buf_mass_1HLT_BB_);
  fill2D(h_mass_1HLT_BE_, buf_mass_1HLT_BE_);
  fill2D(h_mass_1HLT_EE_, buf_mass_1HLT_EE_);
  fill2D(h_mass_Sta_pass_BB_, buf_mass_Sta_pass_BB_);
  fill2D(h_mass_Sta_pass_BE_, buf_mass_Sta_pass_BE_);
  fill2D(h_mass_Sta_pass_EE_, buf_mass_Sta_pass_EE_);
  fill2D(h_mass_Sta_fail_BB_, buf_mass_Sta_fail_BB_);
  fill2D(h_mass_Sta_fail_BE_, buf_mass_Sta_fail_BE_);
  fill2D(h_mass_Sta_fail_EE_, buf_mass_Sta_fail_EE_);
  fill2D(h_mass_Glo_pass_BB_, buf_mass_Glo_pass_BB_);
  fill2D(h_mass_Glo_pass_BE_, buf_mass_Glo_pass_BE_);
  fill2D(h_mass_Glo_pass_EE_, buf_mass_Glo_pass_EE_);
  fill2D(h_mass_Glo_fail_BB_, buf_mass_Glo_fail_BB_);
  fill2D(h_mass_Glo_fail_BE_, buf_mass_Glo_fail_BE_);
  fill2D(h_mass_Glo_fail_EE_, buf_mass_Glo_fail_EE_);
  fill2D(h_mass_ID_fail_BB_,  buf_mass_ID_fail_BB_);
  fill2D(h_mass_ID_fail_BE_,  buf_mass_ID_fail_BE_);
  fill2D(h_mass_ID_fail_EE_,  buf_mass_ID_fail_EE_);

  // PV
  if (h_npv_) {
    for (const auto& [ls, nvs] : buf_npv_)
      for (auto nv : nvs)
        h_npv_->Fill(ls, nv);
  }

  // J/psi flush
  for (int ic = 0; ic < kNCats; ++ic) {
    for (int it = 0; it < kNTPTypes; ++it) {
      for (int ir = 0; ir < kNRegions; ++ir) {
        fill2D(h_jpsi_[ic][it][ir], buf_jpsi_[ic][it][ir]);
      }
    }
  }

  // W
  fill2D(h_w_mt_, buf_w_mt_);
  fill2D(h_w_m2pt_, buf_w_m2pt_);

  edm::LogInfo("ScoutingZCounting")
    << "Evt="<<nEvt_<<" mu="<<nMu_<<" trk="<<nTrk_
    << " tagHLT="<<nTagHLT_<<" pairs="<<nPairs_
    << " 2HLT="<<n2HLT_<<" 1HLT="<<n1HLT_
    << " StaPass="<<nStaPass_<<" StaFail="<<nStaFail_
    << " GloPass="<<nGloPass_<<" GloFail="<<nGloFail_;

  // clear buffers
  buf_mass_2HLT_BB_.clear();    buf_mass_2HLT_BE_.clear();    buf_mass_2HLT_EE_.clear();
  buf_mass_1HLT_BB_.clear();    buf_mass_1HLT_BE_.clear();    buf_mass_1HLT_EE_.clear();
  buf_mass_Sta_pass_BB_.clear(); buf_mass_Sta_pass_BE_.clear(); buf_mass_Sta_pass_EE_.clear();
  buf_mass_Sta_fail_BB_.clear(); buf_mass_Sta_fail_BE_.clear(); buf_mass_Sta_fail_EE_.clear();
  buf_mass_Glo_pass_BB_.clear(); buf_mass_Glo_pass_BE_.clear(); buf_mass_Glo_pass_EE_.clear();
  buf_mass_Glo_fail_BB_.clear(); buf_mass_Glo_fail_BE_.clear(); buf_mass_Glo_fail_EE_.clear();
  buf_mass_ID_fail_BB_.clear();  buf_mass_ID_fail_BE_.clear();  buf_mass_ID_fail_EE_.clear();
  buf_npv_.clear();

  for (int ic = 0; ic < kNCats; ++ic)
    for (int it = 0; it < kNTPTypes; ++it)
      for (int ir = 0; ir < kNRegions; ++ir)
        buf_jpsi_[ic][it][ir].clear();

  buf_w_mt_.clear();
  buf_w_m2pt_.clear();

  edm::LogInfo("ScoutingZCounting")
    << "Finished. runNumber=" << runNumber_
    << ", saveIntermediate=" << (saveIntermediate_ ? "true" : "false");
}

DEFINE_FWK_MODULE(ScoutingZCounting);