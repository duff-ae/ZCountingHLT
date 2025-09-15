// ScoutingZCounting.cc
//
// EDAnalyzer for HLT Scouting data that mirrors the main logic of
// DQMOffline/Lumi/plugins/ZCounting.cc, adapted to Run3Scouting* formats.
//
// Ported elements (analog to original):
//  - PV selection (ndof, |z|, rho, tracksSize)
//  - Tag-and-Probe on muons with BB/BE/EE categories
//  - 2HLT / 1HLT logic (HLT event pass via TriggerResults, per-object HLT match is approximated)
//  - ID-fail category
//  - 2D histos vs lumisection and mass (booked after seeing LS span)
//  - h_npv vs LS
//
// Adaptations / limitations in Scouting:
//  - Exact TriggerObject matching is not stored: we use TriggerResults (event-level pass) + a kinematic/iso proxy
//    for per-object HLT matching. Replace passHLTObjApprox() with a proper matcher if you add TrigObj info.
//  - No separate standalone-muon track collection: we keep Global/ID categories; Sta_* categories are omitted.
//
// Configuration parameters mirror the original names as much as possible.
//
// Author: your base + ChatGPT

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

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/RegexMatch.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "TLorentzVector.h"
#include "TH2D.h"

#include <map>
#include <vector>
#include <string>
#include <limits>
#include <cmath>

namespace {
  constexpr double MUON_MASS  = 0.1056583755;
  constexpr double MUON_BOUND = 0.9; // barrel |eta| < 0.9
}

class ScoutingZCounting : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingZCounting(const edm::ParameterSet&);
  ~ScoutingZCounting() override = default;

  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;

private:
  // ==== Config ====
  // Collections
  edm::EDGetTokenT<std::vector<Run3ScoutingMuon>>   muonToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> pvToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingTrack>>  trackToken_;
  edm::EDGetTokenT<edm::TriggerResults>             trgResultsToken_;

  // Trigger path patterns (regex), like original MuonTriggerNames
  std::vector<std::string> muonTriggerPatterns_;

  // Kinematic cuts
  double PtCutL1_;
  double PtCutL2_;
  double EtaCutL1_;
  double EtaCutL2_;

  // Mass / PV binning
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

  // ID / ISO types
  std::string IDTypestr_;
  std::string IsoTypestr_;
  double      IsoCut_;

  enum IDType  { NoneID=0, LooseID, MediumID, TightID, CustomTightID } IDType_{NoneID};
  enum IsoType { NoneIso=0, TrackerIso, PFIso } IsoType_{NoneIso};

  edm::Service<TFileService> fs_;

  // ==== Histos (created in endJob when LS span is known) ====
  TH2D *h_mass_2HLT_BB_{nullptr}, *h_mass_2HLT_BE_{nullptr}, *h_mass_2HLT_EE_{nullptr};
  TH2D *h_mass_1HLT_BB_{nullptr}, *h_mass_1HLT_BE_{nullptr}, *h_mass_1HLT_EE_{nullptr};
  TH2D *h_mass_ID_fail_BB_{nullptr}, *h_mass_ID_fail_BE_{nullptr}, *h_mass_ID_fail_EE_{nullptr};
  TH2D *h_npv_{nullptr};

  // Buffers per LS (so we can book exact LS axis later)
  std::map<unsigned, std::vector<double>> buf_mass_2HLT_BB_, buf_mass_2HLT_BE_, buf_mass_2HLT_EE_;
  std::map<unsigned, std::vector<double>> buf_mass_1HLT_BB_, buf_mass_1HLT_BE_, buf_mass_1HLT_EE_;
  std::map<unsigned, std::vector<double>> buf_mass_ID_fail_BB_, buf_mass_ID_fail_BE_, buf_mass_ID_fail_EE_;
  std::map<unsigned, unsigned>            buf_npv_;

  // ==== Helpers ====
  static inline bool isCentral(double eta) { return std::abs(eta) < MUON_BOUND; }
  static inline bool oppositeCharge(int q1, int q2) { return (q1 * q2) < 0; }

  bool eventPassesMuonHLT(const edm::Event& iEvent, const edm::TriggerResults& res) const;

  bool passPV(const Run3ScoutingVertex& v) const;
  unsigned selectGoodPVs(const std::vector<Run3ScoutingVertex>& pvs, int& firstGoodIdx) const;

  bool passGlobalMuon(const Run3ScoutingMuon& mu) const; // scouting proxy
  bool passMuonID   (const Run3ScoutingMuon& mu) const;  // scouting proxy
  bool passMuonIso  (const Run3ScoutingMuon& mu) const;  // scouting proxy

  // Per-object HLT proxy (no TriggerObjects in Scouting): kinematic + relIso
  bool passHLTObjApprox(const Run3ScoutingMuon& mu) const;

  // Track quality for track-based probe (optional use)
  bool passTrack(const Run3ScoutingTrack& trk) const;

  void fillBuffered(std::map<unsigned, std::vector<double>>& m, unsigned ls, double mass) {
    m[ls].push_back(mass);
  }

  void bookAllHistosWithLS(unsigned minLS, unsigned maxLS);
};

ScoutingZCounting::ScoutingZCounting(const edm::ParameterSet& iCfg) {
  usesResource("TFileService");

  muonToken_        = consumes<std::vector<Run3ScoutingMuon>>( iCfg.getParameter<edm::InputTag>("muonCollection") );
  pvToken_          = consumes<std::vector<Run3ScoutingVertex>>( iCfg.getParameter<edm::InputTag>("primaryVertexCollection") );
  trackToken_       = consumes<std::vector<Run3ScoutingTrack>>( iCfg.getParameter<edm::InputTag>("trackCollection") );
  trgResultsToken_  = consumes<edm::TriggerResults>( iCfg.getParameter<edm::InputTag>("triggerResults") );

  muonTriggerPatterns_ = iCfg.getParameter<std::vector<std::string>>("MuonTriggerNames");

  PtCutL1_  = iCfg.getUntrackedParameter<double>("PtCutL1", 24.0);
  PtCutL2_  = iCfg.getUntrackedParameter<double>("PtCutL2", 15.0);
  EtaCutL1_ = iCfg.getUntrackedParameter<double>("EtaCutL1", 2.4);
  EtaCutL2_ = iCfg.getUntrackedParameter<double>("EtaCutL2", 2.4);

  MassBin_ = iCfg.getUntrackedParameter<int>("MassBin", 60);
  MassMin_ = iCfg.getUntrackedParameter<double>("MassMin", 60.0);
  MassMax_ = iCfg.getUntrackedParameter<double>("MassMax", 120.0);

  PVBin_ = iCfg.getUntrackedParameter<int>("PVBin", 100);
  PVMin_ = iCfg.getUntrackedParameter<double>("PVMin", 0.);
  PVMax_ = iCfg.getUntrackedParameter<double>("PVMax", 100.);

  VtxNTracksFitCut_ = iCfg.getUntrackedParameter<double>("VtxNTracksFitMin", 0.);
  VtxNdofCut_       = iCfg.getUntrackedParameter<double>("VtxNdofMin", 4.);
  VtxAbsZCut_       = iCfg.getUntrackedParameter<double>("VtxAbsZMax", 24.);
  VtxRhoCut_        = iCfg.getUntrackedParameter<double>("VtxRhoMax", 2.);

  IDTypestr_  = iCfg.getUntrackedParameter<std::string>("IDType", "CustomTight");
  IsoTypestr_ = iCfg.getUntrackedParameter<std::string>("IsoType", "PF-based");
  IsoCut_     = iCfg.getUntrackedParameter<double>("IsoCut", 10.0);

  if      (IDTypestr_ == "Loose")       IDType_ = LooseID;
  else if (IDTypestr_ == "Medium")      IDType_ = MediumID;
  else if (IDTypestr_ == "Tight")       IDType_ = TightID;
  else if (IDTypestr_ == "CustomTight") IDType_ = CustomTightID;
  else                                    IDType_ = NoneID;

  if      (IsoTypestr_ == "Tracker-based") IsoType_ = TrackerIso;
  else if (IsoTypestr_ == "PF-based")      IsoType_ = PFIso;
  else                                       IsoType_ = NoneIso;
}

void ScoutingZCounting::beginJob() {}

bool ScoutingZCounting::eventPassesMuonHLT(const edm::Event& iEvent, const edm::TriggerResults& res) const {
  const edm::TriggerNames& names = iEvent.triggerNames(res);
  for (size_t i = 0; i < names.size(); ++i) {
    const std::string& path = names.triggerName(i);
    // Try each user pattern against this single path
    for (const auto& pat : muonTriggerPatterns_) {
      std::vector<std::string> one{path};
      if (!edm::regexMatch(one, pat).empty()) {
        if (res.accept(i)) return true; // any matching accepted path is enough
        break; // matched a pattern; no need to test others for this path
      }
    }
  }
  return false;
}

bool ScoutingZCounting::passPV(const Run3ScoutingVertex& v) const {
  if (!v.isValidVtx()) return false;
  if (v.tracksSize() < VtxNTracksFitCut_) return false;
  if (v.ndof() < VtxNdofCut_) return false;
  if (std::abs(v.z()) > VtxAbsZCut_) return false;
  const double rho = std::hypot(v.x(), v.y());
  if (rho > VtxRhoCut_) return false;
  return true;
}

unsigned ScoutingZCounting::selectGoodPVs(const std::vector<Run3ScoutingVertex>& pvs, int& firstGoodIdx) const {
  unsigned nvtx = 0; firstGoodIdx = -1;
  for (size_t i = 0; i < pvs.size(); ++i) {
    if (passPV(pvs[i])) { if (firstGoodIdx < 0) firstGoodIdx = (int)i; ++nvtx; }
  }
  return nvtx;
}

bool ScoutingZCounting::passGlobalMuon(const Run3ScoutingMuon& mu) const {
  // Scouting proxy (no inner/outer tracks available):
  return mu.isGlobalMuon();
}

bool ScoutingZCounting::passMuonID(const Run3ScoutingMuon& mu) const {
  switch (IDType_) {
    case CustomTightID:
      return (mu.isGlobalMuon() && mu.nTrackerLayersWithMeasurement() > 5 && mu.normalizedChi2() < 10.0);
    case TightID:
      return (mu.isGlobalMuon() && mu.nTrackerLayersWithMeasurement() > 6 && mu.normalizedChi2() < 5.0);
    case MediumID:
      return (mu.nTrackerLayersWithMeasurement() > 5 && mu.normalizedChi2() < 10.0);
    case LooseID:
      return (mu.nTrackerLayersWithMeasurement() >= 1);
    case NoneID:
    default:
      return true;
  }
}

bool ScoutingZCounting::passMuonIso(const Run3ScoutingMuon& mu) const {
  switch (IsoType_) {
    case TrackerIso: {
      return mu.trackIso() < IsoCut_;
    }
    case PFIso: {
      const double absIso = mu.trackIso() + mu.ecalIso() + mu.hcalIso();
      return absIso < IsoCut_;
    }
    case NoneIso:
    default:
      return true;
  }
}

bool ScoutingZCounting::passHLTObjApprox(const Run3ScoutingMuon& mu) const {
  if (mu.pt() < PtCutL1_) return false;
  if (std::abs(mu.eta()) > EtaCutL1_) return false;
  if (IsoType_ != NoneIso) {
    const double relIso = (mu.trackIso() + mu.ecalIso() + mu.hcalIso()) / std::max(1e-6f, mu.pt());
    if (relIso > 0.15) return false; // typical single-muon HLT proxy
  }
  return true;
}

bool ScoutingZCounting::passTrack(const Run3ScoutingTrack& trk) const {
  // Quality proxy using available Run3ScoutingTrack accessors
  const int nLayers = trk.tk_nTrackerLayersWithMeasurement();
  const int nPix    = trk.tk_nValidPixelHits();
  const float chi2  = trk.tk_chi2();
  const float ndof  = trk.tk_ndof();
  const bool chi2ok = (ndof > 0.f) ? (chi2/ndof < 5.0f) : true;
  if (nLayers < 6) return false;
  if (nPix    < 1) return false;
  if (!chi2ok)     return false;
  return true;
}

void ScoutingZCounting::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  const unsigned ls = iEvent.luminosityBlock();

  // Collections
  edm::Handle<std::vector<Run3ScoutingMuon>>   hMu;
  edm::Handle<std::vector<Run3ScoutingVertex>> hPV;
  edm::Handle<std::vector<Run3ScoutingTrack>>  hTrk;
  edm::Handle<edm::TriggerResults>             hTrg;

  iEvent.getByToken(muonToken_, hMu);
  iEvent.getByToken(pvToken_,   hPV);
  iEvent.getByToken(trackToken_,hTrk);
  iEvent.getByToken(trgResultsToken_, hTrg);

  if (!hMu.isValid() || !hPV.isValid()) { // || !hTrg.isValid()) {
    edm::LogWarning("ScoutingZCounting") << "Invalid handle(s): muons/PV/TriggerResults. Skipping.";
    return;
  }
  if (hMu->empty()) return;

  std::cout << "New Event" << std::endl;

  // Event-level HLT requirement (instead of TriggerObj per-object gate)
  //if (!eventPassesMuonHLT(iEvent, *hTrg)) return;

  // PV selection & count
  int firstPVidx = -1;
  const unsigned nvtx = selectGoodPVs(*hPV, firstPVidx);
  buf_npv_[ls] += nvtx; // per-event fill like original DQM
  const Run3ScoutingVertex* pv = (firstPVidx >= 0 ? &hPV->at(firstPVidx) : nullptr);
  (void)pv; // kept for future extensions if you want tight-with-PV proxy

  // Tag-and-Probe (muon-muon), categories 2HLT / 1HLT / ID_fail, split by BB/BE/EE
  TLorentzVector vTag, vProbe;

  for (size_t i = 0; i < hMu->size(); ++i) {
    const auto& mu1 = hMu->at(i);

    const float pt1  = mu1.pt();
    const float eta1 = mu1.eta();
    const float phi1 = mu1.phi();
    const int   q1   = mu1.charge();

    // Tag selection: kinematics + Global + ID + ISO + HLT proxy
    if (pt1 < PtCutL1_)             continue;
    //if (std::abs(eta1) > EtaCutL1_) continue;
    //if (!(passGlobalMuon(mu1) && passMuonID(mu1) && passMuonIso(mu1))) continue;
    //if (!passHLTObjApprox(mu1))     continue;

    vTag.SetPtEtaPhiM(pt1, eta1, phi1, MUON_MASS);
    const bool tagCentral = isCentral(eta1);

    for (size_t j = 0; j < hMu->size(); ++j) {
      if (j == i) continue;
      const auto& mu2 = hMu->at(j);

      const float pt2  = mu2.pt();
      const float eta2 = mu2.eta();
      const float phi2 = mu2.phi();
      const int   q2   = mu2.charge();

      if (pt2 < PtCutL2_)             continue;
      //if (std::abs(eta2) > EtaCutL2_) continue;
      //if (!oppositeCharge(q1, q2))    continue;

      vProbe.SetPtEtaPhiM(pt2, eta2, phi2, MUON_MASS);
      const double mass = (vTag + vProbe).M();
      if (mass < MassMin_ || mass > MassMax_) continue;

      const bool probeCentral = isCentral(eta2);
      const bool probeGlobal  = passGlobalMuon(mu2);
      const bool probeIDISO   = passMuonID(mu2) && passMuonIso(mu2);
      const bool probeHLT     = passHLTObjApprox(mu2);

      if (probeGlobal && probeIDISO) {
        if (probeHLT) {
          if (tagCentral && probeCentral)          fillBuffered(buf_mass_2HLT_BB_, ls, mass);
          else if (!tagCentral && !probeCentral)   fillBuffered(buf_mass_2HLT_EE_, ls, mass);
          else                                      fillBuffered(buf_mass_2HLT_BE_, ls, mass);
        } else {
          if (tagCentral && probeCentral)          fillBuffered(buf_mass_1HLT_BB_, ls, mass);
          else if (!tagCentral && !probeCentral)   fillBuffered(buf_mass_1HLT_EE_, ls, mass);
          else                                      fillBuffered(buf_mass_1HLT_BE_, ls, mass);
        }
      } else if (probeGlobal) {
        if (tagCentral && probeCentral)            fillBuffered(buf_mass_ID_fail_BB_, ls, mass);
        else if (!tagCentral && !probeCentral)     fillBuffered(buf_mass_ID_fail_EE_, ls, mass);
        else                                        fillBuffered(buf_mass_ID_fail_BE_, ls, mass);
      }
    }
  }

  // Optional: track-based probe (Sta_* like). Without SA muons in Scouting, we skip this section by default.
  // If needed, you can enable a block that loops over hTrk and matches to any global muon within ΔR, then
  // fill pass/fail equivalents with separate histos.
  (void)hTrk;
}

void ScoutingZCounting::bookAllHistosWithLS(unsigned minLS, unsigned maxLS) {
  const int nLumiBins = (maxLS >= minLS) ? int(maxLS - minLS + 1) : 1;
  const double lmin = std::max(0, (int)minLS) - 0.5;
  const double lmax = (int)maxLS + 0.5;

  auto make2D = [&](const char* name, const char* title) -> TH2D* {
    return fs_->make<TH2D>(name, title, nLumiBins, lmin, lmax, MassBin_, MassMin_, MassMax_);
  };

  h_mass_2HLT_BB_ = make2D("h_mass_2HLT_BB", "Both muons pass HLT proxy (BB)");
  h_mass_2HLT_BE_ = make2D("h_mass_2HLT_BE", "Both muons pass HLT proxy (BE)");
  h_mass_2HLT_EE_ = make2D("h_mass_2HLT_EE", "Both muons pass HLT proxy (EE)");

  h_mass_1HLT_BB_ = make2D("h_mass_1HLT_BB", "Only tag passes HLT proxy (BB)");
  h_mass_1HLT_BE_ = make2D("h_mass_1HLT_BE", "Only tag passes HLT proxy (BE)");
  h_mass_1HLT_EE_ = make2D("h_mass_1HLT_EE", "Only tag passes HLT proxy (EE)");

  h_mass_ID_fail_BB_ = make2D("h_mass_ID_fail_BB", "Probe global but fails ID/ISO (BB)");
  h_mass_ID_fail_BE_ = make2D("h_mass_ID_fail_BE", "Probe global but fails ID/ISO (BE)");
  h_mass_ID_fail_EE_ = make2D("h_mass_ID_fail_EE", "Probe global but fails ID/ISO (EE)");

  h_npv_ = fs_->make<TH2D>("h_npv", "Events with valid primary vertex",
                           nLumiBins, lmin, lmax, PVBin_, PVMin_, PVMax_);

  auto setAxes = [&](TH2D* h, const char* y) {
    if (!h) return;
    h->GetXaxis()->SetTitle("luminosity section");
    h->GetYaxis()->SetTitle(y);
  };
  setAxes(h_mass_2HLT_BB_, "tag and probe mass");
  setAxes(h_mass_2HLT_BE_, "tag and probe mass");
  setAxes(h_mass_2HLT_EE_, "tag and probe mass");
  setAxes(h_mass_1HLT_BB_, "tag and probe mass");
  setAxes(h_mass_1HLT_BE_, "tag and probe mass");
  setAxes(h_mass_1HLT_EE_, "tag and probe mass");

  if (h_npv_) {
    h_npv_->GetXaxis()->SetTitle("luminosity section");
    h_npv_->GetYaxis()->SetTitle("number of primary vertices");
  }
}

void ScoutingZCounting::endJob() {
  unsigned minLS = std::numeric_limits<unsigned>::max();
  unsigned maxLS = 0;
  auto scan = [&](const auto& m) { for (const auto& kv : m) { minLS = std::min(minLS, kv.first); maxLS = std::max(maxLS, kv.first); } };
  scan(buf_mass_2HLT_BB_); scan(buf_mass_2HLT_BE_); scan(buf_mass_2HLT_EE_);
  scan(buf_mass_1HLT_BB_); scan(buf_mass_1HLT_BE_); scan(buf_mass_1HLT_EE_);
  scan(buf_mass_ID_fail_BB_); scan(buf_mass_ID_fail_BE_); scan(buf_mass_ID_fail_EE_);
  scan(buf_npv_);

  if (minLS == std::numeric_limits<unsigned>::max()) return; // no data

  bookAllHistosWithLS(minLS, maxLS);

  auto fill2D = [&](TH2D* h, const std::map<unsigned, std::vector<double>>& m) {
    if (!h) return;
    for (const auto& [ls, masses] : m) for (double mval : masses) h->Fill(ls, mval);
  };

  fill2D(h_mass_2HLT_BB_, buf_mass_2HLT_BB_);
  fill2D(h_mass_2HLT_BE_, buf_mass_2HLT_BE_);
  fill2D(h_mass_2HLT_EE_, buf_mass_2HLT_EE_);

  fill2D(h_mass_1HLT_BB_, buf_mass_1HLT_BB_);
  fill2D(h_mass_1HLT_BE_, buf_mass_1HLT_BE_);
  fill2D(h_mass_1HLT_EE_, buf_mass_1HLT_EE_);

  fill2D(h_mass_ID_fail_BB_, buf_mass_ID_fail_BB_);
  fill2D(h_mass_ID_fail_BE_, buf_mass_ID_fail_BE_);
  fill2D(h_mass_ID_fail_EE_, buf_mass_ID_fail_EE_);

  if (h_npv_) { for (const auto& [ls, nv] : buf_npv_) h_npv_->Fill(ls, nv); }

  // clear buffers
  buf_mass_2HLT_BB_.clear(); buf_mass_2HLT_BE_.clear(); buf_mass_2HLT_EE_.clear();
  buf_mass_1HLT_BB_.clear(); buf_mass_1HLT_BE_.clear(); buf_mass_1HLT_EE_.clear();
  buf_mass_ID_fail_BB_.clear(); buf_mass_ID_fail_BE_.clear(); buf_mass_ID_fail_EE_.clear();
  buf_npv_.clear();
}

DEFINE_FWK_MODULE(ScoutingZCounting);
