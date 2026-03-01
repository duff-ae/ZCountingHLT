#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

#include "TTree.h"

#include <vector>
#include <cmath>

// ====================================================================
// Flat ntuplizer for Run3 Scouting + PAT + GEN
//   - no J/psi logic, no efficiencies
//   - just dump everything needed for later reconstruction in Python
// ====================================================================

class ScoutingMCJpsiNtuplizer : public edm::one::EDAnalyzer<> {
public:
    explicit ScoutingMCJpsiNtuplizer(const edm::ParameterSet& cfg);
    ~ScoutingMCJpsiNtuplizer() override = default;

    void beginJob() override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
    // ---- Tokens ----
    edm::EDGetTokenT<std::vector<reco::GenParticle>>  genParticlesToken_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>>  pileupToken_;
    edm::EDGetTokenT<reco::BeamSpot>                  beamspotToken_;

    edm::EDGetTokenT<std::vector<Run3ScoutingMuon>>      scoutingMuonsToken_;
    edm::EDGetTokenT<std::vector<Run3ScoutingTrack>>     scoutingTracksToken_;
    edm::EDGetTokenT<std::vector<Run3ScoutingVertex>>    scoutingVerticesToken_;
    edm::EDGetTokenT<std::vector<Run3ScoutingElectron>>  scoutingElectronsToken_;
    edm::EDGetTokenT<std::vector<Run3ScoutingPhoton>>    scoutingPhotonsToken_;
    edm::EDGetTokenT<std::vector<Run3ScoutingPFJet>>     scoutingPFJetsToken_;
    edm::EDGetTokenT<std::vector<Run3ScoutingParticle>>  scoutingParticlesToken_;

    edm::EDGetTokenT<std::vector<pat::Muon>>          patMuonsToken_;

    // ---- Trees ----
    TTree* t_event_;        // 1 entry per event
    TTree* t_genMu_;        // 1 entry per GEN muon
    TTree* t_scoutMu_;      // 1 entry per scouting muon
    TTree* t_scoutTk_;      // 1 entry per scouting track
    TTree* t_scoutVtx_;     // 1 entry per scouting vertex
    TTree* t_patMu_;        // 1 entry per PAT muon
    TTree* t_scoutEle_;     // 1 entry per scouting electron
    TTree* t_scoutPho_;     // 1 entry per scouting photon
    TTree* t_scoutPFJet_;   // 1 entry per scouting PF jet
    TTree* t_spart_;        // 1 entry per generic Run3ScoutingParticle

    // ---- Event-level branches (shared) ----
    unsigned int ev_run_;
    unsigned int ev_lumi_;
    ULong64_t    ev_event_;
    float        ev_weight_;
    float        ev_nTrueInt_;
    int          ev_nPV_scout_;

    // ---- GEN muon branches ----
    unsigned int g_run_;
    unsigned int g_lumi_;
    ULong64_t    g_event_;
    float        g_nTrueInt_;
    float        g_pt_;
    float        g_eta_;
    float        g_phi_;
    int          g_q_;
    int          g_pdgId_;

    // ---- Scouting muon branches ----
    unsigned int sm_run_;
    unsigned int sm_lumi_;
    ULong64_t    sm_event_;
    float        sm_nTrueInt_;

    float sm_pt_;
    float sm_eta_;
    float sm_phi_;
    int   sm_charge_;
    float sm_normalizedChi2_;
    float sm_ecalIso_;
    float sm_hcalIso_;
    float sm_trackIso_;
    int   sm_nStandAloneMuonMatchedStations_;
    int   sm_nRecoMuonChambers_;
    int   sm_nPixelLayersWithMeasurement_;
    int   sm_nTrackerLayersWithMeasurement_;

    // ---- Scouting track branches ----
    unsigned int st_run_;
    unsigned int st_lumi_;
    ULong64_t    st_event_;
    float        st_nTrueInt_;

    float st_pt_;
    float st_eta_;
    float st_phi_;
    float st_chi2_;
    float st_ndof_;
    int   st_charge_;
    float st_dxy_;
    float st_dz_;
    int   st_nValidPixelHits_;
    int   st_nTrackerLayersWithMeasurement_;
    int   st_nValidStripHits_;

    // ---- Scouting vertex branches ----
    unsigned int sv_run_;
    unsigned int sv_lumi_;
    ULong64_t    sv_event_;
    float        sv_nTrueInt_;

    int   sv_index_;
    float sv_x_;
    float sv_y_;
    float sv_z_;
    float sv_xErr_;
    float sv_yErr_;
    float sv_zErr_;
    float sv_chi2_;
    float sv_ndof_;
    int   sv_nTracks_;
    float sv_Lxy_;
    float sv_LxyErr_;
    float sv_LxySig_;

    // ---- PAT muon branches ----
    unsigned int pm_run_;
    unsigned int pm_lumi_;
    ULong64_t    pm_event_;
    float        pm_nTrueInt_;

    float pm_pt_;
    float pm_eta_;
    float pm_phi_;
    int   pm_charge_;
    int   pm_isGlobal_;
    int   pm_isTracker_;
    int   pm_isStandalone_;
    float pm_normChi2_;
    float pm_dxy_;
    float pm_dz_;
    int   pm_nPixelLayers_;
    int   pm_nTrackerLayers_;
    float pm_relIso04_;

    // ---- Scouting electron branches ----
    unsigned int se_run_;
    unsigned int se_lumi_;
    ULong64_t    se_event_;
    float        se_nTrueInt_;

    float se_pt_;
    float se_eta_;
    float se_phi_;

    // ---- Scouting photon branches ----
    unsigned int spho_run_;
    unsigned int spho_lumi_;
    ULong64_t    spho_event_;
    float        spho_nTrueInt_;

    float spho_pt_;
    float spho_eta_;
    float spho_phi_;

    // ---- Scouting PF jet branches ----
    unsigned int jf_run_;
    unsigned int jf_lumi_;
    ULong64_t    jf_event_;
    float        jf_nTrueInt_;

    float jf_pt_;
    float jf_eta_;
    float jf_phi_;

    // ---- Generic Run3ScoutingParticle branches ----
    unsigned int sp_run_;
    unsigned int sp_lumi_;
    ULong64_t    sp_event_;
    float        sp_weight_;
    float        sp_nTrueInt_;

    float sp_pt_;
    float sp_eta_;
    float sp_phi_;

    int   sp_pdgId_;
    int   sp_vertex_;
    float sp_normChi2_;
    float sp_dz_;
    float sp_dxy_;
    float sp_dzsig_;
    float sp_dxysig_;
    int   sp_lostInnerHits_;
    int   sp_quality_;

    float sp_trk_pt_;
    float sp_trk_eta_;
    float sp_trk_phi_;
    int   sp_relativeTrkVars_;
};


// ====================================================================
// Constructor: take InputTags from cfg, with reasonable defaults
// ====================================================================

ScoutingMCJpsiNtuplizer::ScoutingMCJpsiNtuplizer(const edm::ParameterSet& cfg)
{
    auto getTag = [&](const std::string& name, const edm::InputTag& def) {
        if (cfg.exists(name))
            return cfg.getParameter<edm::InputTag>(name);
        else
            return def;
    };

    // GEN + PU + beam spot
    edm::InputTag genTag = getTag("genParticles", edm::InputTag("prunedGenParticles"));
    genParticlesToken_ = consumes<std::vector<reco::GenParticle>>(genTag);

    edm::InputTag puTag  = getTag("pileup",       edm::InputTag("slimmedAddPileupInfo"));
    pileupToken_ = consumes<std::vector<PileupSummaryInfo>>(puTag);

    edm::InputTag bsTag  = getTag("beamSpot",     edm::InputTag("offlineBeamSpot"));
    beamspotToken_ = consumes<reco::BeamSpot>(bsTag);

    // Scouting muons / tracks / vertices / EGM / PF
    edm::InputTag scoutMuTag = getTag("scoutingMuons",
                                      edm::InputTag("hltScoutingMuonPacker"));
    scoutingMuonsToken_ = consumes<std::vector<Run3ScoutingMuon>>(scoutMuTag);

    edm::InputTag scoutTkTag = getTag("scoutingTracks",
                                      edm::InputTag("hltScoutingTrackPacker"));
    scoutingTracksToken_ = consumes<std::vector<Run3ScoutingTrack>>(scoutTkTag);

    edm::InputTag scoutVtxTag = getTag("scoutingVertices",
                                       edm::InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"));
    scoutingVerticesToken_ = consumes<std::vector<Run3ScoutingVertex>>(scoutVtxTag);

    edm::InputTag scoutEleTag = getTag("scoutingElectrons",
                                       edm::InputTag("hltScoutingEgammaPacker", "electrons"));
    scoutingElectronsToken_ = consumes<std::vector<Run3ScoutingElectron>>(scoutEleTag);

    edm::InputTag scoutPhoTag = getTag("scoutingPhotons",
                                       edm::InputTag("hltScoutingEgammaPacker", "photons"));
    scoutingPhotonsToken_ = consumes<std::vector<Run3ScoutingPhoton>>(scoutPhoTag);

    edm::InputTag scoutPFJetTag = getTag("scoutingPFJets",
                                         edm::InputTag("hltScoutingPFPacker", "pfJets"));
    scoutingPFJetsToken_ = consumes<std::vector<Run3ScoutingPFJet>>(scoutPFJetTag);

    // Generic PF candidates as Run3ScoutingParticle (from PFPacker)
    edm::InputTag scoutPartTag = getTag(
        "scoutingParticles",
        edm::InputTag("hltScoutingPFPacker", "pfCandidates")
    );
    scoutingParticlesToken_ = consumes<std::vector<Run3ScoutingParticle>>(scoutPartTag);

    // PAT muons
    edm::InputTag patMuTag = getTag("patMuons", edm::InputTag("slimmedMuons"));
    patMuonsToken_ = consumes<std::vector<pat::Muon>>(patMuTag);
}


// ====================================================================
// beginJob: book trees
// ====================================================================

void ScoutingMCJpsiNtuplizer::beginJob() {
    edm::Service<TFileService> fs;

    // ---- Event tree ----
    t_event_ = fs->make<TTree>("t_event", "event-level info");
    t_event_->Branch("run",       &ev_run_,      "run/i");
    t_event_->Branch("lumi",      &ev_lumi_,     "lumi/i");
    t_event_->Branch("event",     &ev_event_,    "event/l");
    t_event_->Branch("weight",    &ev_weight_,   "weight/F");
    t_event_->Branch("nTrueInt",  &ev_nTrueInt_, "nTrueInt/F");
    t_event_->Branch("nPV_scout", &ev_nPV_scout_,"nPV_scout/I");

    // ---- GEN muon tree ----
    t_genMu_ = fs->make<TTree>("t_gen_mu", "GEN muons");
    t_genMu_->Branch("run",      &g_run_,      "run/i");
    t_genMu_->Branch("lumi",     &g_lumi_,     "lumi/i");
    t_genMu_->Branch("event",    &g_event_,    "event/l");
    t_genMu_->Branch("nTrueInt", &g_nTrueInt_, "nTrueInt/F");
    t_genMu_->Branch("pt",       &g_pt_,       "pt/F");
    t_genMu_->Branch("eta",      &g_eta_,      "eta/F");
    t_genMu_->Branch("phi",      &g_phi_,      "phi/F");
    t_genMu_->Branch("q",        &g_q_,        "q/I");
    t_genMu_->Branch("pdgId",    &g_pdgId_,    "pdgId/I");

    // ---- Scouting muon tree ----
    t_scoutMu_ = fs->make<TTree>("t_scout_mu", "Run3ScoutingMuon");
    t_scoutMu_->Branch("run",      &sm_run_,      "run/i");
    t_scoutMu_->Branch("lumi",     &sm_lumi_,     "lumi/i");
    t_scoutMu_->Branch("event",    &sm_event_,    "event/l");
    t_scoutMu_->Branch("nTrueInt", &sm_nTrueInt_, "nTrueInt/F");

    t_scoutMu_->Branch("pt",   &sm_pt_,   "pt/F");
    t_scoutMu_->Branch("eta",  &sm_eta_,  "eta/F");
    t_scoutMu_->Branch("phi",  &sm_phi_,  "phi/F");
    t_scoutMu_->Branch("charge", &sm_charge_, "charge/I");

    t_scoutMu_->Branch("normalizedChi2", &sm_normalizedChi2_, "normalizedChi2/F");
    t_scoutMu_->Branch("ecalIso",        &sm_ecalIso_,        "ecalIso/F");
    t_scoutMu_->Branch("hcalIso",        &sm_hcalIso_,        "hcalIso/F");
    t_scoutMu_->Branch("trackIso",       &sm_trackIso_,       "trackIso/F");
    t_scoutMu_->Branch("nStandAloneMuonMatchedStations",
                       &sm_nStandAloneMuonMatchedStations_, "nStandAloneMuonMatchedStations/I");
    t_scoutMu_->Branch("nRecoMuonChambers",
                       &sm_nRecoMuonChambers_, "nRecoMuonChambers/I");
    t_scoutMu_->Branch("nPixelLayersWithMeasurement",
                       &sm_nPixelLayersWithMeasurement_, "nPixelLayersWithMeasurement/I");
    t_scoutMu_->Branch("nTrackerLayersWithMeasurement",
                       &sm_nTrackerLayersWithMeasurement_, "nTrackerLayersWithMeasurement/I");

    // ---- Scouting track tree ----
    t_scoutTk_ = fs->make<TTree>("t_scout_track", "Run3ScoutingTrack");
    t_scoutTk_->Branch("run",      &st_run_,      "run/i");
    t_scoutTk_->Branch("lumi",     &st_lumi_,     "lumi/i");
    t_scoutTk_->Branch("event",    &st_event_,    "event/l");
    t_scoutTk_->Branch("nTrueInt", &st_nTrueInt_, "nTrueInt/F");

    t_scoutTk_->Branch("pt",   &st_pt_,   "pt/F");
    t_scoutTk_->Branch("eta",  &st_eta_,  "eta/F");
    t_scoutTk_->Branch("phi",  &st_phi_,  "phi/F");
    t_scoutTk_->Branch("chi2", &st_chi2_, "chi2/F");
    t_scoutTk_->Branch("ndof", &st_ndof_, "ndof/F");
    t_scoutTk_->Branch("charge", &st_charge_, "charge/I");
    t_scoutTk_->Branch("dxy",     &st_dxy_,    "dxy/F");
    t_scoutTk_->Branch("dz",      &st_dz_,     "dz/F");
    t_scoutTk_->Branch("nValidPixelHits",
                       &st_nValidPixelHits_, "nValidPixelHits/I");
    t_scoutTk_->Branch("nTrackerLayersWithMeasurement",
                       &st_nTrackerLayersWithMeasurement_, "nTrackerLayersWithMeasurement/I");
    t_scoutTk_->Branch("nValidStripHits",
                       &st_nValidStripHits_, "nValidStripHits/I");

    // ---- Scouting vertex tree ----
    t_scoutVtx_ = fs->make<TTree>("t_scout_vtx", "Run3ScoutingVertex");
    t_scoutVtx_->Branch("run",      &sv_run_,      "run/i");
    t_scoutVtx_->Branch("lumi",     &sv_lumi_,     "lumi/i");
    t_scoutVtx_->Branch("event",    &sv_event_,    "event/l");
    t_scoutVtx_->Branch("nTrueInt", &sv_nTrueInt_, "nTrueInt/F");

    t_scoutVtx_->Branch("index", &sv_index_, "index/I");
    t_scoutVtx_->Branch("x",     &sv_x_,     "x/F");
    t_scoutVtx_->Branch("y",     &sv_y_,     "y/F");
    t_scoutVtx_->Branch("z",     &sv_z_,     "z/F");
    t_scoutVtx_->Branch("xErr",  &sv_xErr_,  "xErr/F");
    t_scoutVtx_->Branch("yErr",  &sv_yErr_,  "yErr/F");
    t_scoutVtx_->Branch("zErr",  &sv_zErr_,  "zErr/F");
    t_scoutVtx_->Branch("chi2",  &sv_chi2_,  "chi2/F");
    t_scoutVtx_->Branch("ndof",  &sv_ndof_,  "ndof/F");
    t_scoutVtx_->Branch("nTracks", &sv_nTracks_, "nTracks/I");
    t_scoutVtx_->Branch("Lxy",     &sv_Lxy_,     "Lxy/F");
    t_scoutVtx_->Branch("LxyErr",  &sv_LxyErr_,  "LxyErr/F");
    t_scoutVtx_->Branch("LxySig",  &sv_LxySig_,  "LxySig/F");

    // ---- PAT muon tree ----
    t_patMu_ = fs->make<TTree>("t_pat_mu", "pat::Muon");
    t_patMu_->Branch("run",      &pm_run_,      "run/i");
    t_patMu_->Branch("lumi",     &pm_lumi_,     "lumi/i");
    t_patMu_->Branch("event",    &pm_event_,    "event/l");
    t_patMu_->Branch("nTrueInt", &pm_nTrueInt_, "nTrueInt/F");

    t_patMu_->Branch("pt",   &pm_pt_,   "pt/F");
    t_patMu_->Branch("eta",  &pm_eta_,  "eta/F");
    t_patMu_->Branch("phi",  &pm_phi_,  "phi/F");
    t_patMu_->Branch("charge", &pm_charge_, "charge/I");
    t_patMu_->Branch("isGlobal",     &pm_isGlobal_,     "isGlobal/I");
    t_patMu_->Branch("isTracker",    &pm_isTracker_,    "isTracker/I");
    t_patMu_->Branch("isStandalone", &pm_isStandalone_, "isStandalone/I");
    t_patMu_->Branch("normChi2",     &pm_normChi2_,     "normChi2/F");
    t_patMu_->Branch("dxy",          &pm_dxy_,          "dxy/F");
    t_patMu_->Branch("dz",           &pm_dz_,           "dz/F");
    t_patMu_->Branch("nPixelLayers",   &pm_nPixelLayers_,   "nPixelLayers/I");
    t_patMu_->Branch("nTrackerLayers", &pm_nTrackerLayers_, "nTrackerLayers/I");
    t_patMu_->Branch("relIso04",       &pm_relIso04_,       "relIso04/F");

    // ---- Scouting electron tree ----
    t_scoutEle_ = fs->make<TTree>("t_scout_ele", "Run3ScoutingElectron");
    t_scoutEle_->Branch("run",      &se_run_,      "run/i");
    t_scoutEle_->Branch("lumi",     &se_lumi_,     "lumi/i");
    t_scoutEle_->Branch("event",    &se_event_,    "event/l");
    t_scoutEle_->Branch("nTrueInt", &se_nTrueInt_, "nTrueInt/F");

    t_scoutEle_->Branch("pt",     &se_pt_,     "pt/F");
    t_scoutEle_->Branch("eta",    &se_eta_,    "eta/F");
    t_scoutEle_->Branch("phi",    &se_phi_,    "phi/F");

    // ---- Scouting photon tree ----
    t_scoutPho_ = fs->make<TTree>("t_scout_pho", "Run3ScoutingPhoton");
    t_scoutPho_->Branch("run",      &spho_run_,      "run/i");
    t_scoutPho_->Branch("lumi",     &spho_lumi_,     "lumi/i");
    t_scoutPho_->Branch("event",    &spho_event_,    "event/l");
    t_scoutPho_->Branch("nTrueInt", &spho_nTrueInt_, "nTrueInt/F");

    t_scoutPho_->Branch("pt",   &spho_pt_,   "pt/F");
    t_scoutPho_->Branch("eta",  &spho_eta_,  "eta/F");
    t_scoutPho_->Branch("phi",  &spho_phi_,  "phi/F");

    // ---- Scouting PF jet tree ----
    t_scoutPFJet_ = fs->make<TTree>("t_scout_pfjet", "Run3ScoutingPFJet");
    t_scoutPFJet_->Branch("run",      &jf_run_,      "run/i");
    t_scoutPFJet_->Branch("lumi",     &jf_lumi_,     "lumi/i");
    t_scoutPFJet_->Branch("event",    &jf_event_,    "event/l");
    t_scoutPFJet_->Branch("nTrueInt", &jf_nTrueInt_, "nTrueInt/F");

    t_scoutPFJet_->Branch("pt",   &jf_pt_,   "pt/F");
    t_scoutPFJet_->Branch("eta",  &jf_eta_,  "eta/F");
    t_scoutPFJet_->Branch("phi",  &jf_phi_,  "phi/F");

    // ---- Generic Run3ScoutingParticle tree ----
    t_spart_ = fs->make<TTree>("t_spart", "Run3ScoutingParticle ntuple");

    t_spart_->Branch("run",      &sp_run_,      "run/i");
    t_spart_->Branch("lumi",     &sp_lumi_,     "lumi/i");
    t_spart_->Branch("event",    &sp_event_,    "event/l");
    t_spart_->Branch("weight",   &sp_weight_,   "weight/F");
    t_spart_->Branch("nTrueInt", &sp_nTrueInt_, "nTrueInt/F");

    t_spart_->Branch("pt",    &sp_pt_,    "pt/F");
    t_spart_->Branch("eta",   &sp_eta_,   "eta/F");
    t_spart_->Branch("phi",   &sp_phi_,   "phi/F");

    t_spart_->Branch("pdgId",       &sp_pdgId_,       "pdgId/I");
    t_spart_->Branch("vertex",      &sp_vertex_,      "vertex/I");
    t_spart_->Branch("normChi2",    &sp_normChi2_,    "normChi2/F");
    t_spart_->Branch("dz",          &sp_dz_,          "dz/F");
    t_spart_->Branch("dxy",         &sp_dxy_,         "dxy/F");
    t_spart_->Branch("dzsig",       &sp_dzsig_,       "dzsig/F");
    t_spart_->Branch("dxysig",      &sp_dxysig_,      "dxysig/F");
    t_spart_->Branch("lostInnerHits",&sp_lostInnerHits_,"lostInnerHits/I");
    t_spart_->Branch("quality",     &sp_quality_,     "quality/I");

    t_spart_->Branch("trk_pt",      &sp_trk_pt_,      "trk_pt/F");
    t_spart_->Branch("trk_eta",     &sp_trk_eta_,     "trk_eta/F");
    t_spart_->Branch("trk_phi",     &sp_trk_phi_,     "trk_phi/F");
    t_spart_->Branch("relativeTrkVars", &sp_relativeTrkVars_, "relativeTrkVars/I");
}

// ====================================================================
// analyze: one event
// ====================================================================

void ScoutingMCJpsiNtuplizer::analyze(const edm::Event& iEvent,
                                      const edm::EventSetup&) {
    // ---- input handles ----
    edm::Handle<std::vector<reco::GenParticle>>  h_gen;
    edm::Handle<std::vector<PileupSummaryInfo>>  h_pu;
    edm::Handle<reco::BeamSpot>                  h_bs;

    edm::Handle<std::vector<Run3ScoutingMuon>>      h_scoutMu;
    edm::Handle<std::vector<Run3ScoutingTrack>>     h_scoutTk;
    edm::Handle<std::vector<Run3ScoutingVertex>>    h_scoutVtx;
    edm::Handle<std::vector<Run3ScoutingElectron>>  h_scoutEle;
    edm::Handle<std::vector<Run3ScoutingPhoton>>    h_scoutPho;
    edm::Handle<std::vector<Run3ScoutingPFJet>>     h_scoutPFJet;
    edm::Handle<std::vector<Run3ScoutingParticle>>  h_spart;

    edm::Handle<std::vector<pat::Muon>>          h_patMu;

    iEvent.getByToken(genParticlesToken_,  h_gen);
    iEvent.getByToken(pileupToken_,        h_pu);
    iEvent.getByToken(beamspotToken_,      h_bs);

    iEvent.getByToken(scoutingMuonsToken_,    h_scoutMu);
    iEvent.getByToken(scoutingTracksToken_,   h_scoutTk);
    iEvent.getByToken(scoutingVerticesToken_, h_scoutVtx);
    iEvent.getByToken(scoutingElectronsToken_,h_scoutEle);
    iEvent.getByToken(scoutingPhotonsToken_,  h_scoutPho);
    iEvent.getByToken(scoutingPFJetsToken_,   h_scoutPFJet);
    iEvent.getByToken(scoutingParticlesToken_, h_spart);

    iEvent.getByToken(patMuonsToken_,      h_patMu);

    // ---- Event / PU info ----
    ev_run_   = iEvent.id().run();
    ev_lumi_  = iEvent.luminosityBlock();
    ev_event_ = iEvent.id().event();
    ev_weight_ = 1.0f;

    double nTrue = 0.0;
    if (h_pu.isValid()) {
        for (const auto& pu : *h_pu) {
            if (pu.getBunchCrossing() == 0) {
                nTrue = pu.getTrueNumInteractions();
                break;
            }
        }
    }
    ev_nTrueInt_ = static_cast<float>(nTrue);

    ev_nPV_scout_ = h_scoutVtx.isValid() ? int(h_scoutVtx->size()) : 0;

    // event tree
    t_event_->Fill();

    // convenience for per-object trees:
    auto fillCommon = [&](unsigned int& run,
                          unsigned int& lumi,
                          ULong64_t&    evt,
                          float&        nTrueOut) {
        run   = ev_run_;
        lumi  = ev_lumi_;
        evt   = ev_event_;
        nTrueOut = ev_nTrueInt_;
    };

    // ---- GEN muons ----
    if (h_gen.isValid()) {
        for (const auto& gp : *h_gen) {
            int id = gp.pdgId();
            if (std::abs(id) != 13) continue;

            fillCommon(g_run_, g_lumi_, g_event_, g_nTrueInt_);
            g_pt_    = gp.pt();
            g_eta_   = gp.eta();
            g_phi_   = gp.phi();
            g_pdgId_ = id;
            g_q_     = (id > 0 ? -1 : +1); // 13 -> -1, -13 -> +1

            t_genMu_->Fill();
        }
    }

    // ---- Scouting muons ----
    if (h_scoutMu.isValid()) {
        for (const auto& mu : *h_scoutMu) {
            fillCommon(sm_run_, sm_lumi_, sm_event_, sm_nTrueInt_);

            sm_pt_   = mu.pt();
            sm_eta_  = mu.eta();
            sm_phi_  = mu.phi();
            sm_charge_ = mu.charge();

            sm_normalizedChi2_ = mu.normalizedChi2();
            sm_ecalIso_  = mu.ecalIso();
            sm_hcalIso_  = mu.hcalIso();
            sm_trackIso_ = mu.trackIso();
            sm_nStandAloneMuonMatchedStations_ = mu.nStandAloneMuonMatchedStations();
            sm_nRecoMuonChambers_              = mu.nRecoMuonChambers();
            sm_nPixelLayersWithMeasurement_    = mu.nPixelLayersWithMeasurement();
            sm_nTrackerLayersWithMeasurement_  = mu.nTrackerLayersWithMeasurement();

            t_scoutMu_->Fill();
        }
    }

    // ---- Scouting tracks ----
    if (h_scoutTk.isValid()) {
        for (const auto& tk : *h_scoutTk) {
            fillCommon(st_run_, st_lumi_, st_event_, st_nTrueInt_);

            st_pt_   = tk.tk_pt();
            st_eta_  = tk.tk_eta();
            st_phi_  = tk.tk_phi();
            st_chi2_ = tk.tk_chi2();
            st_ndof_ = tk.tk_ndof();
            st_charge_ = tk.tk_charge();
            st_dxy_  = tk.tk_dxy();
            st_dz_   = tk.tk_dz();
            st_nValidPixelHits_ = tk.tk_nValidPixelHits();
            st_nTrackerLayersWithMeasurement_ = tk.tk_nTrackerLayersWithMeasurement();
            st_nValidStripHits_ = tk.tk_nValidStripHits();

            t_scoutTk_->Fill();
        }
    }

    // ---- Scouting vertices ----
    if (h_scoutVtx.isValid()) {
        double bsX = 0.0, bsY = 0.0;
        bool haveBS = h_bs.isValid();
        if (haveBS) {
            bsX = h_bs->x0();
            bsY = h_bs->y0();
        }

        for (size_t i = 0; i < h_scoutVtx->size(); ++i) {
            const auto& v = h_scoutVtx->at(i);

            fillCommon(sv_run_, sv_lumi_, sv_event_, sv_nTrueInt_);

            sv_index_ = static_cast<int>(i);
            sv_x_   = v.x();
            sv_y_   = v.y();
            sv_z_   = v.z();
            sv_xErr_ = v.xError();
            sv_yErr_ = v.yError();
            sv_zErr_ = v.zError();
            sv_chi2_ = v.chi2();
            sv_ndof_ = v.ndof();
            sv_nTracks_ = v.tracksSize();

            double dx = v.x() - (haveBS ? bsX : 0.0);
            double dy = v.y() - (haveBS ? bsY : 0.0);
            double Lxy = std::sqrt(dx*dx + dy*dy);
            double LxyErr = 0.0;
            if (Lxy > 0.) {
                double var = dx*dx*v.xError()*v.xError() +
                             dy*dy*v.yError()*v.yError();
                LxyErr = std::sqrt(var) / Lxy;
            }
            double LxySig = (LxyErr > 0.) ? (Lxy / LxyErr) : -1.0;

            sv_Lxy_    = static_cast<float>(Lxy);
            sv_LxyErr_ = static_cast<float>(LxyErr);
            sv_LxySig_ = static_cast<float>(LxySig);

            t_scoutVtx_->Fill();
        }
    }

    // ---- PAT muons ----
    if (h_patMu.isValid()) {
        for (const auto& mu : *h_patMu) {
            fillCommon(pm_run_, pm_lumi_, pm_event_, pm_nTrueInt_);

            pm_pt_   = mu.pt();
            pm_eta_  = mu.eta();
            pm_phi_  = mu.phi();
            pm_charge_ = mu.charge();
            pm_isGlobal_     = mu.isGlobalMuon()    ? 1 : 0;
            pm_isTracker_    = mu.isTrackerMuon()   ? 1 : 0;
            pm_isStandalone_ = mu.isStandAloneMuon()? 1 : 0;

            pm_normChi2_   = -1.f;
            pm_dxy_        = 0.f;
            pm_dz_         = 0.f;
            pm_nPixelLayers_   = 0;
            pm_nTrackerLayers_ = 0;

            if (mu.globalTrack().isNonnull()) {
                pm_normChi2_ = float(mu.globalTrack()->normalizedChi2());
            }
            if (mu.innerTrack().isNonnull()) {
                pm_dxy_ = float(mu.innerTrack()->dxy());
                pm_dz_  = float(mu.innerTrack()->dz());
                pm_nPixelLayers_   =
                    mu.innerTrack()->hitPattern().pixelLayersWithMeasurement();
                pm_nTrackerLayers_ =
                    mu.innerTrack()->hitPattern().trackerLayersWithMeasurement();
            }

            // standard PF relIso (if available)
            pm_relIso04_ = float(
                (mu.pfIsolationR04().sumChargedHadronPt +
                 std::max(0.f, mu.pfIsolationR04().sumNeutralHadronEt +
                              mu.pfIsolationR04().sumPhotonEt -
                              0.5f * mu.pfIsolationR04().sumPUPt)) / mu.pt()
            );

            t_patMu_->Fill();
        }
    }

    // ---- Scouting electrons ----
    if (h_scoutEle.isValid()) {
        for (const auto& e : *h_scoutEle) {
            fillCommon(se_run_, se_lumi_, se_event_, se_nTrueInt_);

            se_pt_   = e.pt();
            se_eta_  = e.eta();
            se_phi_  = e.phi();

            t_scoutEle_->Fill();
        }
    }

    // ---- Scouting photons ----
    if (h_scoutPho.isValid()) {
        for (const auto& g : *h_scoutPho) {
            fillCommon(spho_run_, spho_lumi_, spho_event_, spho_nTrueInt_);

            spho_pt_   = g.pt();
            spho_eta_  = g.eta();
            spho_phi_  = g.phi();

            t_scoutPho_->Fill();
        }
    }

    // ---- Scouting PF jets ----
    if (h_scoutPFJet.isValid()) {
        for (const auto& j : *h_scoutPFJet) {
            fillCommon(jf_run_, jf_lumi_, jf_event_, jf_nTrueInt_);

            jf_pt_   = j.pt();
            jf_eta_  = j.eta();
            jf_phi_  = j.phi();

            t_scoutPFJet_->Fill();
        }
    }

    // ---- Generic Run3ScoutingParticle ntuple (PF candidates etc.) ----
    if (h_spart.isValid()) {
        for (const auto& p : *h_spart) {
            // event-level
            sp_run_      = ev_run_;
            sp_lumi_     = ev_lumi_;
            sp_event_    = ev_event_;
            sp_weight_   = 1.0f;
            sp_nTrueInt_ = ev_nTrueInt_;

            // kinematics
            sp_pt_   = p.pt();
            sp_eta_  = p.eta();
            sp_phi_  = p.phi();

            // ID / track info
            sp_pdgId_    = p.pdgId();
            sp_vertex_   = p.vertex();
            sp_normChi2_ = p.normchi2();
            sp_dz_       = p.dz();
            sp_dxy_      = p.dxy();
            sp_dzsig_    = p.dzsig();
            sp_dxysig_   = p.dxysig();
            sp_lostInnerHits_ = static_cast<int>(p.lostInnerHits());
            sp_quality_       = static_cast<int>(p.quality());

            sp_trk_pt_  = p.trk_pt();
            sp_trk_eta_ = p.trk_eta();
            sp_trk_phi_ = p.trk_phi();
            sp_relativeTrkVars_ = p.relative_trk_vars() ? 1 : 0;

            t_spart_->Fill();
        }
    }
}

DEFINE_FWK_MODULE(ScoutingMCJpsiNtuplizer);