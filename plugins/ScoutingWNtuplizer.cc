// ScoutingWNtuplizer.cc
// MC ntuple: ONLY fiducial GEN W->mu nu events.
// Save GEN truth + event-level scouting + ALL scouting muons (+ DR to GEN mu).
// No PAT. No histograms. No internal indices in output.

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "TTree.h"

#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>

namespace {
  inline float dR2(float eta1, float phi1, float eta2, float phi2) {
    float dphi = reco::deltaPhi(phi1, phi2);
    float deta = eta1 - eta2;
    return dphi*dphi + deta*deta;
  }

  inline const reco::GenParticle* findHardW(const std::vector<reco::GenParticle>& gps) {
    const reco::GenParticle* w = nullptr;
    for (const auto& gp : gps) {
      if (std::abs(gp.pdgId()) != 24) continue;
      if (!gp.statusFlags().fromHardProcess()) continue;
      if (!w) w = &gp;
      else {
        bool w_isLast  = w->statusFlags().isLastCopy();
        bool gp_isLast = gp.statusFlags().isLastCopy();
        if (!w_isLast && gp_isLast) w = &gp;
      }
    }
    return w;
  }

  struct GenWmunu {
    bool ok{false};
    int  w_charge{0};
    float w_pt{0.f}, w_eta{0.f}, w_phi{0.f}, w_mass{0.f};
    float mu_pt{0.f}, mu_eta{0.f}, mu_phi{0.f}; int mu_charge{0};
    float nu_pt{0.f}, nu_eta{0.f}, nu_phi{0.f};
  };

  inline GenWmunu selectGenWmunu(const std::vector<reco::GenParticle>& gps) {
    GenWmunu out;
    const reco::GenParticle* w = findHardW(gps);
    if (!w) return out;

    const reco::GenParticle* mu = nullptr;
    const reco::GenParticle* nu = nullptr;

    for (const auto& gp : gps) {
      const auto& f = gp.statusFlags();
      if (!f.fromHardProcess() || !f.isLastCopy()) continue;
      int id = std::abs(gp.pdgId());
      if (id == 13) {
        if (!mu || gp.pt() > mu->pt()) mu = &gp;
      } else if (id == 14) {
        if (!nu || gp.pt() > nu->pt()) nu = &gp;
      }
    }
    if (!(mu && nu)) return out;

    out.ok = true;
    out.w_charge = (w->pdgId() > 0) ? +1 : -1;
    out.w_pt = w->pt(); out.w_eta = w->eta(); out.w_phi = w->phi(); out.w_mass = w->mass();

    out.mu_pt = mu->pt(); out.mu_eta = mu->eta(); out.mu_phi = mu->phi();
    // pdgId: mu- = +13, mu+ = -13
    out.mu_charge = (mu->pdgId() > 0) ? -1 : +1;

    out.nu_pt = nu->pt(); out.nu_eta = nu->eta(); out.nu_phi = nu->phi();
    return out;
  }
} // namespace

class ScoutingWNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingWNtuplizer(const edm::ParameterSet&);
  ~ScoutingWNtuplizer() override = default;

  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  void reset_();

  float fid_mu_pt_min_{3.f};
  float fid_mu_abs_eta_max_{2.4f};
  float dr_match_tight_{0.10f};
  float dr_match_loose_{0.20f};

  // Scouting
  edm::EDGetTokenT<std::vector<Run3ScoutingMuon>>    scMuToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> scPvToken_;
  edm::EDGetTokenT<double> scMetPtToken_;
  edm::EDGetTokenT<double> scMetPhiToken_;
  edm::EDGetTokenT<double> rhoToken_;

  // MC
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>>  puToken_;
  edm::EDGetTokenT<GenEventInfoProduct>             genEvtToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>>  genToken_;
  edm::EDGetTokenT<std::vector<reco::GenMET>>       genMetToken_;

  edm::Service<TFileService> fs_;
  TTree* t_{nullptr};

  // Event/meta
  uint32_t run_{0}, lumi_{0};
  uint64_t event_{0};
  float genWeightRaw_{1.f}; // optional diagnostic
  float nTrueInt_{-1.f};
  int   nPV_{0};
  float rho_{-1.f};

  // GEN (W->munu)
  int   gen_w_charge_{0};
  float gen_w_pt_{0.f}, gen_w_eta_{0.f}, gen_w_phi_{0.f}, gen_w_mass_{0.f};
  float gen_mu_pt_{0.f}, gen_mu_eta_{0.f}, gen_mu_phi_{0.f}; int gen_mu_charge_{0};
  float gen_nu_pt_{0.f}, gen_nu_eta_{0.f}, gen_nu_phi_{0.f};
  float gen_met_pt_{-1.f}, gen_met_phi_{0.f};

  // Scouting event
  float sc_met_pt_{-1.f}, sc_met_phi_{0.f};

  // Scouting muons (ALL)
  int sc_nMuon_{0};

  std::vector<float> sc_mu_pt_;
  std::vector<float> sc_mu_eta_;
  std::vector<float> sc_mu_phi_;
  std::vector<int>   sc_mu_charge_;

  std::vector<int>   sc_mu_isGlobal_;
  std::vector<int>   sc_mu_isTracker_;
  std::vector<float> sc_mu_normChi2_;
  std::vector<float> sc_mu_dxy_;
  std::vector<float> sc_mu_dz_;
  std::vector<float> sc_mu_trackIso_;
  std::vector<float> sc_mu_ecalIso_;
  std::vector<float> sc_mu_hcalIso_;

  // Matching per muon to GEN mu
  std::vector<float> sc_mu_drToGen_;
  std::vector<int>   sc_mu_isMatchedTight_;
  std::vector<int>   sc_mu_isMatchedLoose_;
};

ScoutingWNtuplizer::ScoutingWNtuplizer(const edm::ParameterSet& cfg) {
  usesResource("TFileService");

  fid_mu_pt_min_      = cfg.getParameter<double>("fid_mu_pt_min");
  fid_mu_abs_eta_max_ = cfg.getParameter<double>("fid_mu_abs_eta_max");
  dr_match_tight_     = cfg.getParameter<double>("dr_match_tight");
  dr_match_loose_     = cfg.getParameter<double>("dr_match_loose");

  scMuToken_     = consumes<std::vector<Run3ScoutingMuon>>(cfg.getParameter<edm::InputTag>("scoutingMuons"));
  scPvToken_     = consumes<std::vector<Run3ScoutingVertex>>(cfg.getParameter<edm::InputTag>("scoutingPrimaryVertices"));
  scMetPtToken_  = consumes<double>(cfg.getParameter<edm::InputTag>("pfMetPtValue"));
  scMetPhiToken_ = consumes<double>(cfg.getParameter<edm::InputTag>("pfMetPhiValue"));
  rhoToken_      = consumes<double>(cfg.getParameter<edm::InputTag>("rhoValue"));

  puToken_     = consumes<std::vector<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("pileup"));
  genEvtToken_ = consumes<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("genEventInfo"));
  genToken_    = consumes<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"));
  genMetToken_ = consumes<std::vector<reco::GenMET>>(cfg.getParameter<edm::InputTag>("genMetTrue"));
}

void ScoutingWNtuplizer::beginJob() {
  t_ = fs_->make<TTree>("Events", "Fiducial GEN W->munu + scouting snapshot");

  t_->Branch("run", &run_, "run/i");
  t_->Branch("lumi", &lumi_, "lumi/i");
  t_->Branch("event", &event_, "event/l");

  t_->Branch("genWeightRaw", &genWeightRaw_, "genWeightRaw/F");
  t_->Branch("nTrueInt", &nTrueInt_, "nTrueInt/F");
  t_->Branch("nPV", &nPV_, "nPV/I");
  t_->Branch("rho", &rho_, "rho/F");

  t_->Branch("gen_w_charge", &gen_w_charge_, "gen_w_charge/I");
  t_->Branch("gen_w_pt", &gen_w_pt_, "gen_w_pt/F");
  t_->Branch("gen_w_eta", &gen_w_eta_, "gen_w_eta/F");
  t_->Branch("gen_w_phi", &gen_w_phi_, "gen_w_phi/F");
  t_->Branch("gen_w_mass", &gen_w_mass_, "gen_w_mass/F");

  t_->Branch("gen_mu_pt", &gen_mu_pt_, "gen_mu_pt/F");
  t_->Branch("gen_mu_eta", &gen_mu_eta_, "gen_mu_eta/F");
  t_->Branch("gen_mu_phi", &gen_mu_phi_, "gen_mu_phi/F");
  t_->Branch("gen_mu_charge", &gen_mu_charge_, "gen_mu_charge/I");

  t_->Branch("gen_nu_pt", &gen_nu_pt_, "gen_nu_pt/F");
  t_->Branch("gen_nu_eta", &gen_nu_eta_, "gen_nu_eta/F");
  t_->Branch("gen_nu_phi", &gen_nu_phi_, "gen_nu_phi/F");

  t_->Branch("gen_met_pt", &gen_met_pt_, "gen_met_pt/F");
  t_->Branch("gen_met_phi", &gen_met_phi_, "gen_met_phi/F");

  t_->Branch("sc_met_pt", &sc_met_pt_, "sc_met_pt/F");
  t_->Branch("sc_met_phi", &sc_met_phi_, "sc_met_phi/F");

  t_->Branch("sc_nMuon", &sc_nMuon_, "sc_nMuon/I");

  t_->Branch("sc_mu_pt", &sc_mu_pt_);
  t_->Branch("sc_mu_eta", &sc_mu_eta_);
  t_->Branch("sc_mu_phi", &sc_mu_phi_);
  t_->Branch("sc_mu_charge", &sc_mu_charge_);

  t_->Branch("sc_mu_isGlobal", &sc_mu_isGlobal_);
  t_->Branch("sc_mu_isTracker", &sc_mu_isTracker_);
  t_->Branch("sc_mu_normChi2", &sc_mu_normChi2_);
  t_->Branch("sc_mu_dxy", &sc_mu_dxy_);
  t_->Branch("sc_mu_dz", &sc_mu_dz_);
  t_->Branch("sc_mu_trackIso", &sc_mu_trackIso_);
  t_->Branch("sc_mu_ecalIso", &sc_mu_ecalIso_);
  t_->Branch("sc_mu_hcalIso", &sc_mu_hcalIso_);

  t_->Branch("sc_mu_drToGen", &sc_mu_drToGen_);
  t_->Branch("sc_mu_isMatchedTight", &sc_mu_isMatchedTight_);
  t_->Branch("sc_mu_isMatchedLoose", &sc_mu_isMatchedLoose_);
}

void ScoutingWNtuplizer::reset_() {
  run_=lumi_=0; event_=0;
  genWeightRaw_=1.f; nTrueInt_=-1.f; nPV_=0; rho_=-1.f;

  gen_w_charge_=0; gen_w_pt_=gen_w_eta_=gen_w_phi_=gen_w_mass_=0.f;
  gen_mu_pt_=gen_mu_eta_=gen_mu_phi_=0.f; gen_mu_charge_=0;
  gen_nu_pt_=gen_nu_eta_=gen_nu_phi_=0.f;
  gen_met_pt_=-1.f; gen_met_phi_=0.f;

  sc_met_pt_=-1.f; sc_met_phi_=0.f;

  sc_nMuon_=0;

  sc_mu_pt_.clear();
  sc_mu_eta_.clear();
  sc_mu_phi_.clear();
  sc_mu_charge_.clear();

  sc_mu_isGlobal_.clear();
  sc_mu_isTracker_.clear();
  sc_mu_normChi2_.clear();
  sc_mu_dxy_.clear();
  sc_mu_dz_.clear();
  sc_mu_trackIso_.clear();
  sc_mu_ecalIso_.clear();
  sc_mu_hcalIso_.clear();

  sc_mu_drToGen_.clear();
  sc_mu_isMatchedTight_.clear();
  sc_mu_isMatchedLoose_.clear();
}

void ScoutingWNtuplizer::analyze(const edm::Event& ev, const edm::EventSetup&) {
  reset_();

  run_ = ev.id().run();
  lumi_ = ev.luminosityBlock();
  event_ = ev.id().event();

  // weight (diagnostic only)
  edm::Handle<GenEventInfoProduct> hGenEvt;
  ev.getByToken(genEvtToken_, hGenEvt);
  if (hGenEvt.isValid()) genWeightRaw_ = static_cast<float>(hGenEvt->weight());

  // PU
  edm::Handle<std::vector<PileupSummaryInfo>> hPU;
  ev.getByToken(puToken_, hPU);
  if (hPU.isValid()) {
    for (const auto& pu : *hPU) {
      if (pu.getBunchCrossing() == 0) { nTrueInt_ = pu.getTrueNumInteractions(); break; }
    }
  }

  // scouting event: PV, rho, MET
  edm::Handle<std::vector<Run3ScoutingVertex>> hPV;
  ev.getByToken(scPvToken_, hPV);
  if (hPV.isValid()) nPV_ = (int)hPV->size();

  edm::Handle<double> hRho;
  ev.getByToken(rhoToken_, hRho);
  if (hRho.isValid()) rho_ = (float)(*hRho);

  edm::Handle<double> hMetPt, hMetPhi;
  ev.getByToken(scMetPtToken_, hMetPt);
  ev.getByToken(scMetPhiToken_, hMetPhi);
  if (hMetPt.isValid())  sc_met_pt_  = (float)(*hMetPt);
  if (hMetPhi.isValid()) sc_met_phi_ = (float)(*hMetPhi);

  // GEN selection: W->munu
  edm::Handle<std::vector<reco::GenParticle>> hGen;
  ev.getByToken(genToken_, hGen);
  if (!hGen.isValid()) return;

  GenWmunu g = selectGenWmunu(*hGen);
  if (!g.ok) return;

  // fiducial GEN mu
  if (!(g.mu_pt > fid_mu_pt_min_ && std::abs(g.mu_eta) < fid_mu_abs_eta_max_)) return;

  gen_w_charge_ = g.w_charge;
  gen_w_pt_ = g.w_pt; gen_w_eta_ = g.w_eta; gen_w_phi_ = g.w_phi; gen_w_mass_ = g.w_mass;

  gen_mu_pt_ = g.mu_pt; gen_mu_eta_ = g.mu_eta; gen_mu_phi_ = g.mu_phi; gen_mu_charge_ = g.mu_charge;
  gen_nu_pt_ = g.nu_pt; gen_nu_eta_ = g.nu_eta; gen_nu_phi_ = g.nu_phi;

  edm::Handle<std::vector<reco::GenMET>> hGenMet;
  ev.getByToken(genMetToken_, hGenMet);
  if (hGenMet.isValid() && !hGenMet->empty()) {
    gen_met_pt_ = hGenMet->at(0).pt();
    gen_met_phi_ = hGenMet->at(0).phi();
  }

  // scouting muons: save ALL + per-muon DR to GEN mu
  edm::Handle<std::vector<Run3ScoutingMuon>> hMu;
  ev.getByToken(scMuToken_, hMu);
  if (hMu.isValid()) {
    sc_nMuon_ = (int)hMu->size();

    sc_mu_pt_.reserve(sc_nMuon_);
    sc_mu_eta_.reserve(sc_nMuon_);
    sc_mu_phi_.reserve(sc_nMuon_);
    sc_mu_charge_.reserve(sc_nMuon_);

    sc_mu_isGlobal_.reserve(sc_nMuon_);
    sc_mu_isTracker_.reserve(sc_nMuon_);
    sc_mu_normChi2_.reserve(sc_nMuon_);
    sc_mu_dxy_.reserve(sc_nMuon_);
    sc_mu_dz_.reserve(sc_nMuon_);
    sc_mu_trackIso_.reserve(sc_nMuon_);
    sc_mu_ecalIso_.reserve(sc_nMuon_);
    sc_mu_hcalIso_.reserve(sc_nMuon_);

    sc_mu_drToGen_.reserve(sc_nMuon_);
    sc_mu_isMatchedTight_.reserve(sc_nMuon_);
    sc_mu_isMatchedLoose_.reserve(sc_nMuon_);

    for (const auto& m : *hMu) {
      sc_mu_pt_.push_back(m.pt());
      sc_mu_eta_.push_back(m.eta());
      sc_mu_phi_.push_back(m.phi());
      sc_mu_charge_.push_back(m.charge());

      sc_mu_isGlobal_.push_back(m.isGlobalMuon() ? 1 : 0);
      sc_mu_isTracker_.push_back(m.isTrackerMuon() ? 1 : 0);
      sc_mu_normChi2_.push_back(m.normalizedChi2());
      sc_mu_dxy_.push_back(m.trk_dxy());
      sc_mu_dz_.push_back(m.trk_dz());
      sc_mu_trackIso_.push_back(m.trackIso());
      sc_mu_ecalIso_.push_back(m.ecalIso());
      sc_mu_hcalIso_.push_back(m.hcalIso());

      float dr = 999.f;
      if (m.pt() > 0.f) dr = std::sqrt(dR2(gen_mu_eta_, gen_mu_phi_, m.eta(), m.phi()));
      sc_mu_drToGen_.push_back(dr);
      sc_mu_isMatchedTight_.push_back((dr < dr_match_tight_) ? 1 : 0);
      sc_mu_isMatchedLoose_.push_back((dr < dr_match_loose_) ? 1 : 0);
    }
  }

  t_->Fill();
}

DEFINE_FWK_MODULE(ScoutingWNtuplizer);