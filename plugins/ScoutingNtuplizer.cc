// ScoutingMuonTrackNtuplizer.cc
// Extended ntuplizer: Run3ScoutingMuon / Track / Electron / Photon / PFJet / Particle / Vertex
// + MET, rho, TriggerResults -> ROOT TTree "Events".

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"

#include "DataFormats/Common/interface/TriggerResults.h"

#include "TTree.h"

#include <vector>

class ScoutingNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingNtuplizer(const edm::ParameterSet&);

  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  // --- input tokens
  edm::EDGetTokenT<std::vector<Run3ScoutingMuon>>     muonToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingTrack>>    trackToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingElectron>> eleToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingPhoton>>   phoToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingPFJet>>    jetToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> pfcandToken_;

  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>>   pvToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>>   dvNoVtxToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>>   dvVtxToken_;

  edm::EDGetTokenT<double>                           metPtToken_;
  edm::EDGetTokenT<double>                           metPhiToken_;
  edm::EDGetTokenT<double>                           rhoToken_;

  edm::EDGetTokenT<edm::TriggerResults>              trigToken_;

  // toggles for optional collections
  bool storeMuons_;
  bool storeTracks_;
  bool storeElectrons_;
  bool storePhotons_;
  bool storePFJets_;
  bool storePFCands_;
  bool storeVertices_;
  bool storeMET_;
  bool storeHLT_;

  // --- output
  edm::Service<TFileService> fs_;
  TTree* tree_{nullptr};

  // event ids
  unsigned int       run_{0}, lumi_{0};
  unsigned long long event_{0};

  // ===== MUONS =====
  std::vector<float> mu_pt_, mu_eta_, mu_phi_;
  std::vector<int>   mu_charge_;
  std::vector<float> mu_normChi2_, mu_ecalIso_, mu_hcalIso_, mu_trackIso_;
  std::vector<int>   mu_nPixLayMeas_, mu_nTrkLayMeas_;
  std::vector<int>   mu_nStandAlone_, mu_nRecoStations_;
  std::vector<float> mu_m_;
  std::vector<int>   mu_type_;
  std::vector<int>   mu_nValidStandAloneHits_;
  std::vector<int>   mu_nValidRecoHits_;
  std::vector<int>   mu_nValidPixelHits_;
  std::vector<int>   mu_nValidStripHits_;
  std::vector<int>   mu_nRecoMuonChambers_;
  std::vector<int>   mu_nRecoMuonChambersCSCorDT_;
  std::vector<int>   mu_nRecoMuonMatches_;
  std::vector<int>   mu_nRecoMuonExpectedMatchedStations_;
  std::vector<int>   mu_recoMuonStationMask_;
  std::vector<int>   mu_nRecoMuonMatchedRPCLayers_;
  std::vector<int>   mu_recoMuonRPCLayerMask_;

  // ===== TRACKS =====
  std::vector<float> tk_pt_, tk_eta_, tk_phi_;
  std::vector<float> tk_chi2_, tk_ndof_;
  std::vector<int>   tk_charge_;
  std::vector<float> tk_dxy_, tk_dz_;
  std::vector<int>   tk_nValidPixelHits_, tk_nTrkLayMeas_;
  std::vector<int>   tk_nValidStripHits_;
  std::vector<float> tk_qoverp_;
  std::vector<float> tk_lambda_;
  std::vector<float> tk_dxyErr_;
  std::vector<float> tk_dzErr_;
  std::vector<float> tk_qoverpErr_;
  std::vector<float> tk_lambdaErr_;
  std::vector<float> tk_phiErr_;
  std::vector<float> tk_dsz_;
  std::vector<float> tk_dszErr_;

  // ===== ELECTRONS =====
  std::vector<float> el_pt_, el_eta_, el_phi_, el_m_;
  std::vector<float> el_rawE_, el_preshowerE_, el_corrEcalEnergyError_;
  std::vector<float> el_dEtaIn_, el_dPhiIn_, el_sigmaIetaIeta_, el_hOverE_, el_ooEMOop_;
  std::vector<int>   el_missingHits_;
  std::vector<float> el_ecalIso_, el_hcalIso_, el_trackIso_;
  std::vector<float> el_r9_, el_sMin_, el_sMaj_;
  std::vector<int>   el_seedId_, el_nClusters_, el_nCrystals_;

  // ===== PHOTONS =====
  std::vector<float> ph_pt_, ph_eta_, ph_phi_, ph_m_;
  std::vector<float> ph_rawE_, ph_preshowerE_, ph_corrEcalEnergyError_;
  std::vector<float> ph_sigmaIetaIeta_, ph_hOverE_;
  std::vector<float> ph_ecalIso_, ph_hcalIso_, ph_trkIso_;
  std::vector<float> ph_r9_, ph_sMin_, ph_sMaj_;
  std::vector<int>   ph_seedId_, ph_nClusters_, ph_nCrystals_;

  // ===== PF JETS =====
  std::vector<float> jet_pt_, jet_eta_, jet_phi_, jet_m_;
  std::vector<float> jet_area_;
  std::vector<float> jet_chHadE_, jet_neuHadE_, jet_photonE_, jet_electronE_, jet_muonE_;
  std::vector<float> jet_HFHadE_, jet_HFEME_, jet_HOEnergy_;
  std::vector<int>   jet_chHadMult_, jet_neuHadMult_;
  std::vector<int>   jet_photonMult_, jet_electronMult_, jet_muonMult_;
  std::vector<int>   jet_HFHadMult_, jet_HFEMMult_;
  std::vector<float> jet_csv_, jet_mvaDisc_;

  // ===== PF PARTICLES =====
  std::vector<float> pfc_pt_, pfc_eta_, pfc_phi_;
  std::vector<int>   pfc_pdgId_;
  std::vector<float> pfc_normchi2_;
  std::vector<float> pfc_dz_, pfc_dxy_, pfc_dzsig_, pfc_dxysig_;
  std::vector<int>   pfc_lostInnerHits_, pfc_quality_;

  // ===== VERTICES =====
  // primary vertices
  std::vector<float> pv_x_, pv_y_, pv_z_;
  std::vector<float> pv_zError_, pv_xError_, pv_yError_;
  std::vector<int>   pv_tracksSize_;
  std::vector<float> pv_chi2_, pv_ndof_;
  std::vector<int>   pv_isValid_;
  std::vector<float> pv_xyCov_, pv_xzCov_, pv_yzCov_;

  // displaced vertices from MuonPackerNoVtx
  std::vector<float> dvNoVtx_x_, dvNoVtx_y_, dvNoVtx_z_;
  std::vector<float> dvNoVtx_zError_, dvNoVtx_xError_, dvNoVtx_yError_;
  std::vector<int>   dvNoVtx_tracksSize_;
  std::vector<float> dvNoVtx_chi2_, dvNoVtx_ndof_;
  std::vector<int>   dvNoVtx_isValid_;
  std::vector<float> dvNoVtx_xyCov_, dvNoVtx_xzCov_, dvNoVtx_yzCov_;

  // displaced vertices from MuonPackerVtx
  std::vector<float> dvVtx_x_, dvVtx_y_, dvVtx_z_;
  std::vector<float> dvVtx_zError_, dvVtx_xError_, dvVtx_yError_;
  std::vector<int>   dvVtx_tracksSize_;
  std::vector<float> dvVtx_chi2_, dvVtx_ndof_;
  std::vector<int>   dvVtx_isValid_;
  std::vector<float> dvVtx_xyCov_, dvVtx_xzCov_, dvVtx_yzCov_;

  // ===== MET & RHO =====
  std::vector<double> met_pt_, met_phi_;
  double              rho_{0.};

  // ===== Trigger bits =====
  std::vector<unsigned char> hlt_bits_;
};


// ======================================================================
// ctor
// ======================================================================

ScoutingNtuplizer::ScoutingNtuplizer(const edm::ParameterSet& cfg) {
  usesResource("TFileService");

  storeMuons_     = cfg.getUntrackedParameter<bool>("storeMuons",     true);
  storeTracks_    = cfg.getUntrackedParameter<bool>("storeTracks",    true);
  storeElectrons_ = cfg.getUntrackedParameter<bool>("storeElectrons", true);
  storePhotons_   = cfg.getUntrackedParameter<bool>("storePhotons",   true);
  storePFJets_    = cfg.getUntrackedParameter<bool>("storePFJets",    true);
  storePFCands_   = cfg.getUntrackedParameter<bool>("storePFCands",   true);
  storeVertices_  = cfg.getUntrackedParameter<bool>("storeVertices",  true);
  storeMET_       = cfg.getUntrackedParameter<bool>("storeMET",       true);
  storeHLT_       = cfg.getUntrackedParameter<bool>("storeHLT",       true);

  muonToken_  = consumes<std::vector<Run3ScoutingMuon>>(
      cfg.getParameter<edm::InputTag>("muonCollection"));

  trackToken_ = consumes<std::vector<Run3ScoutingTrack>>(
      cfg.getParameter<edm::InputTag>("trackCollection"));

  eleToken_   = consumes<std::vector<Run3ScoutingElectron>>(
      cfg.getParameter<edm::InputTag>("electronCollection"));

  phoToken_   = consumes<std::vector<Run3ScoutingPhoton>>(
      cfg.getParameter<edm::InputTag>("photonCollection"));

  jetToken_   = consumes<std::vector<Run3ScoutingPFJet>>(
      cfg.getParameter<edm::InputTag>("pfJetCollection"));

  pfcandToken_ = consumes<std::vector<Run3ScoutingParticle>>(
      cfg.getParameter<edm::InputTag>("pfCandidateCollection"));

  pvToken_ = consumes<std::vector<Run3ScoutingVertex>>(
      cfg.getParameter<edm::InputTag>("primaryVertexCollection"));

  dvNoVtxToken_ = consumes<std::vector<Run3ScoutingVertex>>(
      cfg.getParameter<edm::InputTag>("displacedVertexCollectionNoVtx"));

  dvVtxToken_ = consumes<std::vector<Run3ScoutingVertex>>(
      cfg.getParameter<edm::InputTag>("displacedVertexCollectionVtx"));

  metPtToken_  = consumes<double>(
      cfg.getUntrackedParameter<edm::InputTag>("pfMetPtValue",
                                               edm::InputTag("hltScoutingPFPacker","pfMetPt")));

  metPhiToken_ = consumes<double>(
      cfg.getUntrackedParameter<edm::InputTag>("pfMetPhiValue",
                                               edm::InputTag("hltScoutingPFPacker","pfMetPhi")));

  rhoToken_ = consumes<double>(
      cfg.getUntrackedParameter<edm::InputTag>("rhoValue",
                                               edm::InputTag("hltScoutingPFPacker","rho")));

  trigToken_ = consumes<edm::TriggerResults>(
      cfg.getUntrackedParameter<edm::InputTag>("triggerResults",
                                               edm::InputTag("TriggerResults","","HLT")));
}


// ======================================================================
// beginJob: define branches
// ======================================================================

void ScoutingNtuplizer::beginJob() {
  tree_ = fs_->make<TTree>("Events", "Extended HLT Scouting ntuple");

  // event ids
  tree_->Branch("run",   &run_,   "run/i");
  tree_->Branch("lumi",  &lumi_,  "lumi/i");
  tree_->Branch("event", &event_, "event/l");

  // ----- muons -----
  if (storeMuons_) {
    tree_->Branch("mu_pt",     &mu_pt_);
    tree_->Branch("mu_eta",    &mu_eta_);
    tree_->Branch("mu_phi",    &mu_phi_);
    tree_->Branch("mu_charge", &mu_charge_);

    tree_->Branch("mu_normalizedChi2", &mu_normChi2_);
    tree_->Branch("mu_ecalIso",        &mu_ecalIso_);
    tree_->Branch("mu_hcalIso",        &mu_hcalIso_);
    tree_->Branch("mu_trackIso",       &mu_trackIso_);

    tree_->Branch("mu_nPixelLayersWithMeasurement",   &mu_nPixLayMeas_);
    tree_->Branch("mu_nTrackerLayersWithMeasurement", &mu_nTrkLayMeas_);
    tree_->Branch("mu_nStandAloneMuonMatchedStations",&mu_nStandAlone_);
    tree_->Branch("mu_nRecoMuonMatchedStations",      &mu_nRecoStations_);

    // extra muon branches
    tree_->Branch("mu_m",                            &mu_m_);
    tree_->Branch("mu_type",                         &mu_type_);
    tree_->Branch("mu_nValidStandAloneMuonHits",     &mu_nValidStandAloneHits_);
    tree_->Branch("mu_nValidRecoMuonHits",           &mu_nValidRecoHits_);
    tree_->Branch("mu_nValidPixelHits",              &mu_nValidPixelHits_);
    tree_->Branch("mu_nValidStripHits",              &mu_nValidStripHits_);
    tree_->Branch("mu_nRecoMuonChambers",            &mu_nRecoMuonChambers_);
    tree_->Branch("mu_nRecoMuonChambersCSCorDT",     &mu_nRecoMuonChambersCSCorDT_);
    tree_->Branch("mu_nRecoMuonMatches",             &mu_nRecoMuonMatches_);
    tree_->Branch("mu_nRecoMuonExpectedMatchedStations", &mu_nRecoMuonExpectedMatchedStations_);
    tree_->Branch("mu_recoMuonStationMask",          &mu_recoMuonStationMask_);
    tree_->Branch("mu_nRecoMuonMatchedRPCLayers",    &mu_nRecoMuonMatchedRPCLayers_);
    tree_->Branch("mu_recoMuonRPCLayerMask",         &mu_recoMuonRPCLayerMask_);
  }

  if (storeTracks_) {
    // ----- tracks -----
    tree_->Branch("tk_pt",     &tk_pt_);
    tree_->Branch("tk_eta",    &tk_eta_);
    tree_->Branch("tk_phi",    &tk_phi_);
    tree_->Branch("tk_chi2",   &tk_chi2_);
    tree_->Branch("tk_ndof",   &tk_ndof_);
    tree_->Branch("tk_charge", &tk_charge_);
    tree_->Branch("tk_dxy",    &tk_dxy_);
    tree_->Branch("tk_dz",     &tk_dz_);
    tree_->Branch("tk_nValidPixelHits",            &tk_nValidPixelHits_);
    tree_->Branch("tk_nTrackerLayersWithMeasurement", &tk_nTrkLayMeas_);
    // extra track branches
    tree_->Branch("tk_nValidStripHits", &tk_nValidStripHits_);
    tree_->Branch("tk_qoverp",          &tk_qoverp_);
    tree_->Branch("tk_lambda",          &tk_lambda_);
    tree_->Branch("tk_dxyErr",          &tk_dxyErr_);
    tree_->Branch("tk_dzErr",           &tk_dzErr_);
    tree_->Branch("tk_qoverpErr",       &tk_qoverpErr_);
    tree_->Branch("tk_lambdaErr",       &tk_lambdaErr_);
    tree_->Branch("tk_phiErr",          &tk_phiErr_);
    tree_->Branch("tk_dsz",             &tk_dsz_);
    tree_->Branch("tk_dszErr",          &tk_dszErr_);
  }

  // ----- electrons -----
  if (storeElectrons_) {
    tree_->Branch("el_pt",  &el_pt_);
    tree_->Branch("el_eta", &el_eta_);
    tree_->Branch("el_phi", &el_phi_);
    tree_->Branch("el_m",   &el_m_);

    tree_->Branch("el_rawEnergy",           &el_rawE_);
    tree_->Branch("el_preshowerEnergy",     &el_preshowerE_);
    tree_->Branch("el_corrEcalEnergyError", &el_corrEcalEnergyError_);

    tree_->Branch("el_dEtaIn",        &el_dEtaIn_);
    tree_->Branch("el_dPhiIn",        &el_dPhiIn_);
    tree_->Branch("el_sigmaIetaIeta", &el_sigmaIetaIeta_);
    tree_->Branch("el_hOverE",        &el_hOverE_);
    tree_->Branch("el_ooEMOop",       &el_ooEMOop_);

    tree_->Branch("el_missingHits", &el_missingHits_);

    tree_->Branch("el_ecalIso",  &el_ecalIso_);
    tree_->Branch("el_hcalIso",  &el_hcalIso_);
    tree_->Branch("el_trackIso", &el_trackIso_);

    tree_->Branch("el_r9",      &el_r9_);
    tree_->Branch("el_sMin",    &el_sMin_);
    tree_->Branch("el_sMaj",    &el_sMaj_);
    tree_->Branch("el_seedId",  &el_seedId_);
    tree_->Branch("el_nClusters",&el_nClusters_);
    tree_->Branch("el_nCrystals",&el_nCrystals_);
  }

  // ----- photons -----
  if (storePhotons_) {
    tree_->Branch("ph_pt",  &ph_pt_);
    tree_->Branch("ph_eta", &ph_eta_);
    tree_->Branch("ph_phi", &ph_phi_);
    tree_->Branch("ph_m",   &ph_m_);

    tree_->Branch("ph_rawEnergy",           &ph_rawE_);
    tree_->Branch("ph_preshowerEnergy",     &ph_preshowerE_);
    tree_->Branch("ph_corrEcalEnergyError", &ph_corrEcalEnergyError_);

    tree_->Branch("ph_sigmaIetaIeta", &ph_sigmaIetaIeta_);
    tree_->Branch("ph_hOverE",        &ph_hOverE_);
    tree_->Branch("ph_ecalIso",       &ph_ecalIso_);
    tree_->Branch("ph_hcalIso",       &ph_hcalIso_);
    tree_->Branch("ph_trkIso",        &ph_trkIso_);

    tree_->Branch("ph_r9",      &ph_r9_);
    tree_->Branch("ph_sMin",    &ph_sMin_);
    tree_->Branch("ph_sMaj",    &ph_sMaj_);
    tree_->Branch("ph_seedId",  &ph_seedId_);
    tree_->Branch("ph_nClusters",&ph_nClusters_);
    tree_->Branch("ph_nCrystals",&ph_nCrystals_);
  }

  // ----- PF jets -----
  if (storePFJets_) {
    tree_->Branch("jet_pt",   &jet_pt_);
    tree_->Branch("jet_eta",  &jet_eta_);
    tree_->Branch("jet_phi",  &jet_phi_);
    tree_->Branch("jet_m",    &jet_m_);
    tree_->Branch("jet_area", &jet_area_);

    tree_->Branch("jet_chargedHadronEnergy", &jet_chHadE_);
    tree_->Branch("jet_neutralHadronEnergy", &jet_neuHadE_);
    tree_->Branch("jet_photonEnergy",        &jet_photonE_);
    tree_->Branch("jet_electronEnergy",      &jet_electronE_);
    tree_->Branch("jet_muonEnergy",          &jet_muonE_);
    tree_->Branch("jet_HFHadronEnergy",      &jet_HFHadE_);
    tree_->Branch("jet_HFEMEnergy",          &jet_HFEME_);
    tree_->Branch("jet_HOEnergy",            &jet_HOEnergy_);

    tree_->Branch("jet_chargedHadronMultiplicity", &jet_chHadMult_);
    tree_->Branch("jet_neutralHadronMultiplicity", &jet_neuHadMult_);
    tree_->Branch("jet_photonMultiplicity",        &jet_photonMult_);
    tree_->Branch("jet_electronMultiplicity",      &jet_electronMult_);
    tree_->Branch("jet_muonMultiplicity",          &jet_muonMult_);
    tree_->Branch("jet_HFHadronMultiplicity",      &jet_HFHadMult_);
    tree_->Branch("jet_HFEMMultiplicity",          &jet_HFEMMult_);

    tree_->Branch("jet_csv",         &jet_csv_);
    tree_->Branch("jet_mvaDiscr",    &jet_mvaDisc_);
  }

  // ----- PF candidates -----
  if (storePFCands_) {
    tree_->Branch("pfc_pt",   &pfc_pt_);
    tree_->Branch("pfc_eta",  &pfc_eta_);
    tree_->Branch("pfc_phi",  &pfc_phi_);
    tree_->Branch("pfc_pdgId",&pfc_pdgId_);

    tree_->Branch("pfc_normchi2",    &pfc_normchi2_);
    tree_->Branch("pfc_dz",          &pfc_dz_);
    tree_->Branch("pfc_dxy",         &pfc_dxy_);
    tree_->Branch("pfc_dzsig",       &pfc_dzsig_);
    tree_->Branch("pfc_dxysig",      &pfc_dxysig_);
    tree_->Branch("pfc_lostInnerHits",&pfc_lostInnerHits_);
    tree_->Branch("pfc_quality",     &pfc_quality_);
  }

  // ----- primary vertices -----
  if (storeVertices_) {
    tree_->Branch("pv_x",         &pv_x_);
    tree_->Branch("pv_y",         &pv_y_);
    tree_->Branch("pv_z",         &pv_z_);
    tree_->Branch("pv_zError",    &pv_zError_);
    tree_->Branch("pv_xError",    &pv_xError_);
    tree_->Branch("pv_yError",    &pv_yError_);
    tree_->Branch("pv_tracksSize",&pv_tracksSize_);
    tree_->Branch("pv_chi2",      &pv_chi2_);
    tree_->Branch("pv_ndof",      &pv_ndof_);
    tree_->Branch("pv_isValid",   &pv_isValid_);
    tree_->Branch("pv_xyCov",     &pv_xyCov_);
    tree_->Branch("pv_xzCov",     &pv_xzCov_);
    tree_->Branch("pv_yzCov",     &pv_yzCov_);

    // ----- displaced vertices (NoVtx) -----
    tree_->Branch("dvNoVtx_x",         &dvNoVtx_x_);
    tree_->Branch("dvNoVtx_y",         &dvNoVtx_y_);
    tree_->Branch("dvNoVtx_z",         &dvNoVtx_z_);
    tree_->Branch("dvNoVtx_zError",    &dvNoVtx_zError_);
    tree_->Branch("dvNoVtx_xError",    &dvNoVtx_xError_);
    tree_->Branch("dvNoVtx_yError",    &dvNoVtx_yError_);
    tree_->Branch("dvNoVtx_tracksSize",&dvNoVtx_tracksSize_);
    tree_->Branch("dvNoVtx_chi2",      &dvNoVtx_chi2_);
    tree_->Branch("dvNoVtx_ndof",      &dvNoVtx_ndof_);
    tree_->Branch("dvNoVtx_isValid",   &dvNoVtx_isValid_);
    tree_->Branch("dvNoVtx_xyCov",     &dvNoVtx_xyCov_);
    tree_->Branch("dvNoVtx_xzCov",     &dvNoVtx_xzCov_);
    tree_->Branch("dvNoVtx_yzCov",     &dvNoVtx_yzCov_);

    // ----- displaced vertices (Vtx) -----
    tree_->Branch("dvVtx_x",         &dvVtx_x_);
    tree_->Branch("dvVtx_y",         &dvVtx_y_);
    tree_->Branch("dvVtx_z",         &dvVtx_z_);
    tree_->Branch("dvVtx_zError",    &dvVtx_zError_);
    tree_->Branch("dvVtx_xError",    &dvVtx_xError_);
    tree_->Branch("dvVtx_yError",    &dvVtx_yError_);
    tree_->Branch("dvVtx_tracksSize",&dvVtx_tracksSize_);
    tree_->Branch("dvVtx_chi2",      &dvVtx_chi2_);
    tree_->Branch("dvVtx_ndof",      &dvVtx_ndof_);
    tree_->Branch("dvVtx_isValid",   &dvVtx_isValid_);
    tree_->Branch("dvVtx_xyCov",     &dvVtx_xyCov_);
    tree_->Branch("dvVtx_xzCov",     &dvVtx_xzCov_);
    tree_->Branch("dvVtx_yzCov",     &dvVtx_yzCov_);
  }

  if (storeMET_) {
    // ----- MET & rho -----
    tree_->Branch("met_pt",  &met_pt_);
    tree_->Branch("met_phi", &met_phi_);
    tree_->Branch("rho",     &rho_, "rho/D");
  }

  if (storeHLT_) {
    // ----- Trigger bits -----
    tree_->Branch("hlt_bits", &hlt_bits_);
  }
}


// ======================================================================
// analyze
// ======================================================================

void ScoutingNtuplizer::analyze(const edm::Event& ev, const edm::EventSetup&) {
  // clear per-event containers
  mu_pt_.clear();   mu_eta_.clear();   mu_phi_.clear();
  mu_charge_.clear();
  mu_normChi2_.clear(); mu_ecalIso_.clear(); mu_hcalIso_.clear(); mu_trackIso_.clear();
  mu_nPixLayMeas_.clear(); mu_nTrkLayMeas_.clear();
  mu_nStandAlone_.clear(); mu_nRecoStations_.clear();
  mu_m_.clear();
  mu_type_.clear();
  mu_nValidStandAloneHits_.clear();
  mu_nValidRecoHits_.clear();
  mu_nValidPixelHits_.clear();
  mu_nValidStripHits_.clear();
  mu_nRecoMuonChambers_.clear();
  mu_nRecoMuonChambersCSCorDT_.clear();
  mu_nRecoMuonMatches_.clear();
  mu_nRecoMuonExpectedMatchedStations_.clear();
  mu_recoMuonStationMask_.clear();
  mu_nRecoMuonMatchedRPCLayers_.clear();
  mu_recoMuonRPCLayerMask_.clear();

  tk_pt_.clear(); tk_eta_.clear(); tk_phi_.clear();
  tk_chi2_.clear(); tk_ndof_.clear();
  tk_charge_.clear();
  tk_dxy_.clear(); tk_dz_.clear();
  tk_nValidPixelHits_.clear(); tk_nTrkLayMeas_.clear();
  tk_nValidStripHits_.clear();
  tk_qoverp_.clear();
  tk_lambda_.clear();
  tk_dxyErr_.clear();
  tk_dzErr_.clear();
  tk_qoverpErr_.clear();
  tk_lambdaErr_.clear();
  tk_phiErr_.clear();
  tk_dsz_.clear();
  tk_dszErr_.clear();

  el_pt_.clear(); el_eta_.clear(); el_phi_.clear(); el_m_.clear();
  el_rawE_.clear(); el_preshowerE_.clear(); el_corrEcalEnergyError_.clear();
  el_dEtaIn_.clear(); el_dPhiIn_.clear(); el_sigmaIetaIeta_.clear(); el_hOverE_.clear(); el_ooEMOop_.clear();
  el_missingHits_.clear();
  el_ecalIso_.clear(); el_hcalIso_.clear(); el_trackIso_.clear();
  el_r9_.clear(); el_sMin_.clear(); el_sMaj_.clear();
  el_seedId_.clear(); el_nClusters_.clear(); el_nCrystals_.clear();

  ph_pt_.clear(); ph_eta_.clear(); ph_phi_.clear(); ph_m_.clear();
  ph_rawE_.clear(); ph_preshowerE_.clear(); ph_corrEcalEnergyError_.clear();
  ph_sigmaIetaIeta_.clear(); ph_hOverE_.clear();
  ph_ecalIso_.clear(); ph_hcalIso_.clear(); ph_trkIso_.clear();
  ph_r9_.clear(); ph_sMin_.clear(); ph_sMaj_.clear();
  ph_seedId_.clear(); ph_nClusters_.clear(); ph_nCrystals_.clear();

  jet_pt_.clear(); jet_eta_.clear(); jet_phi_.clear(); jet_m_.clear();
  jet_area_.clear();
  jet_chHadE_.clear(); jet_neuHadE_.clear(); jet_photonE_.clear();
  jet_electronE_.clear(); jet_muonE_.clear();
  jet_HFHadE_.clear(); jet_HFEME_.clear(); jet_HOEnergy_.clear();
  jet_chHadMult_.clear(); jet_neuHadMult_.clear();
  jet_photonMult_.clear(); jet_electronMult_.clear(); jet_muonMult_.clear();
  jet_HFHadMult_.clear(); jet_HFEMMult_.clear();
  jet_csv_.clear(); jet_mvaDisc_.clear();

  pfc_pt_.clear(); pfc_eta_.clear(); pfc_phi_.clear();
  pfc_pdgId_.clear();
  pfc_normchi2_.clear();
  pfc_dz_.clear(); pfc_dxy_.clear(); pfc_dzsig_.clear(); pfc_dxysig_.clear();
  pfc_lostInnerHits_.clear(); pfc_quality_.clear();

  pv_x_.clear(); pv_y_.clear(); pv_z_.clear();
  pv_zError_.clear(); pv_xError_.clear(); pv_yError_.clear();
  pv_tracksSize_.clear(); pv_chi2_.clear(); pv_ndof_.clear();
  pv_isValid_.clear();
  pv_xyCov_.clear(); pv_xzCov_.clear(); pv_yzCov_.clear();

  dvNoVtx_x_.clear(); dvNoVtx_y_.clear(); dvNoVtx_z_.clear();
  dvNoVtx_zError_.clear(); dvNoVtx_xError_.clear(); dvNoVtx_yError_.clear();
  dvNoVtx_tracksSize_.clear(); dvNoVtx_chi2_.clear(); dvNoVtx_ndof_.clear();
  dvNoVtx_isValid_.clear();
  dvNoVtx_xyCov_.clear(); dvNoVtx_xzCov_.clear(); dvNoVtx_yzCov_.clear();

  dvVtx_x_.clear(); dvVtx_y_.clear(); dvVtx_z_.clear();
  dvVtx_zError_.clear(); dvVtx_xError_.clear(); dvVtx_yError_.clear();
  dvVtx_tracksSize_.clear(); dvVtx_chi2_.clear(); dvVtx_ndof_.clear();
  dvVtx_isValid_.clear();
  dvVtx_xyCov_.clear(); dvVtx_xzCov_.clear(); dvVtx_yzCov_.clear();

  met_pt_.clear(); met_phi_.clear();
  rho_ = 0.0;
  hlt_bits_.clear();

  // event ids
  run_  = ev.id().run();
  lumi_ = ev.luminosityBlock();
  event_= ev.id().event();

  // read inputs
  edm::Handle<std::vector<Run3ScoutingMuon>>     hMu;
  edm::Handle<std::vector<Run3ScoutingTrack>>    hTk;
  edm::Handle<std::vector<Run3ScoutingElectron>> hEle;
  edm::Handle<std::vector<Run3ScoutingPhoton>>   hPho;
  edm::Handle<std::vector<Run3ScoutingPFJet>>    hJet;
  edm::Handle<std::vector<Run3ScoutingParticle>> hPFC;

  edm::Handle<std::vector<Run3ScoutingVertex>> hPV;
  edm::Handle<std::vector<Run3ScoutingVertex>> hDVNoVtx;
  edm::Handle<std::vector<Run3ScoutingVertex>> hDVVtx;

  edm::Handle<double> hMetPt, hMetPhi, hRho;
  edm::Handle<edm::TriggerResults> hTrig;

  ev.getByToken(muonToken_,  hMu);
  ev.getByToken(trackToken_, hTk);
  ev.getByToken(eleToken_,   hEle);
  ev.getByToken(phoToken_,   hPho);
  ev.getByToken(jetToken_,   hJet);
  ev.getByToken(pfcandToken_,hPFC);

  ev.getByToken(pvToken_,      hPV);
  ev.getByToken(dvNoVtxToken_, hDVNoVtx);
  ev.getByToken(dvVtxToken_,   hDVVtx);

  ev.getByToken(metPtToken_,  hMetPt);
  ev.getByToken(metPhiToken_, hMetPhi);
  ev.getByToken(rhoToken_,    hRho);

  ev.getByToken(trigToken_,   hTrig);

  // ----- muons -----
  if (hMu.isValid()) {
    for (const auto& mu : *hMu) {
      if (mu.pt() > 0) {
        mu_pt_.push_back(mu.pt());
        mu_eta_.push_back(mu.eta());
        mu_phi_.push_back(mu.phi());
        mu_charge_.push_back(mu.charge());
        mu_normChi2_.push_back(mu.normalizedChi2());
        mu_ecalIso_.push_back(mu.ecalIso());
        mu_hcalIso_.push_back(mu.hcalIso());
        mu_trackIso_.push_back(mu.trackIso());
        mu_nPixLayMeas_.push_back(mu.nPixelLayersWithMeasurement());
        mu_nTrkLayMeas_.push_back(mu.nTrackerLayersWithMeasurement());
        mu_nStandAlone_.push_back(mu.nStandAloneMuonMatchedStations());
        mu_nRecoStations_.push_back(mu.nRecoMuonMatchedStations());
        mu_m_.push_back(mu.m());
        mu_type_.push_back(mu.type());
        mu_nValidStandAloneHits_.push_back(mu.nValidStandAloneMuonHits());
        mu_nValidRecoHits_.push_back(mu.nValidRecoMuonHits());
        mu_nValidPixelHits_.push_back(mu.nValidPixelHits());
        mu_nValidStripHits_.push_back(mu.nValidStripHits());
        mu_nRecoMuonChambers_.push_back(mu.nRecoMuonChambers());
        mu_nRecoMuonChambersCSCorDT_.push_back(mu.nRecoMuonChambersCSCorDT());
        mu_nRecoMuonMatches_.push_back(mu.nRecoMuonMatches());
        mu_nRecoMuonExpectedMatchedStations_.push_back(mu.nRecoMuonExpectedMatchedStations());
        mu_recoMuonStationMask_.push_back(mu.recoMuonStationMask());
        mu_nRecoMuonMatchedRPCLayers_.push_back(mu.nRecoMuonMatchedRPCLayers());
        mu_recoMuonRPCLayerMask_.push_back(mu.recoMuonRPClayerMask());

      }
    }
  }

  // ----- tracks -----
  if (hTk.isValid()) {
    for (const auto& tk : *hTk) {
      if (tk.tk_pt() > 0) {
        tk_pt_.push_back(tk.tk_pt());
        tk_eta_.push_back(tk.tk_eta());
        tk_phi_.push_back(tk.tk_phi());
        tk_chi2_.push_back(tk.tk_chi2());
        tk_ndof_.push_back(tk.tk_ndof());
        tk_charge_.push_back(tk.tk_charge());
        tk_dxy_.push_back(tk.tk_dxy());
        tk_dz_.push_back(tk.tk_dz());
        tk_nValidPixelHits_.push_back(tk.tk_nValidPixelHits());
        tk_nTrkLayMeas_.push_back(tk.tk_nTrackerLayersWithMeasurement());
        tk_nValidStripHits_.push_back(tk.tk_nValidStripHits());
        tk_qoverp_.push_back(tk.tk_qoverp());
        tk_lambda_.push_back(tk.tk_lambda());
        tk_dxyErr_.push_back(tk.tk_dxy_Error());
        tk_dzErr_.push_back(tk.tk_dz_Error());
        tk_qoverpErr_.push_back(tk.tk_qoverp_Error());
        tk_lambdaErr_.push_back(tk.tk_lambda_Error());
        tk_phiErr_.push_back(tk.tk_phi_Error());
        tk_dsz_.push_back(tk.tk_dsz());
        tk_dszErr_.push_back(tk.tk_dsz_Error());
      }
    }
  }

  // ----- electrons -----
  if (hEle.isValid()) {
    for (const auto& el : *hEle) {
      if (el.pt() > 0) {
        el_pt_.push_back(el.pt());
        el_eta_.push_back(el.eta());
        el_phi_.push_back(el.phi());
        el_m_.push_back(el.m());

        el_rawE_.push_back(el.rawEnergy());
        el_preshowerE_.push_back(el.preshowerEnergy());
        el_corrEcalEnergyError_.push_back(el.corrEcalEnergyError());

        el_dEtaIn_.push_back(el.dEtaIn());
        el_dPhiIn_.push_back(el.dPhiIn());
        el_sigmaIetaIeta_.push_back(el.sigmaIetaIeta());
        el_hOverE_.push_back(el.hOverE());
        el_ooEMOop_.push_back(el.ooEMOop());

        el_missingHits_.push_back(el.missingHits());

        el_ecalIso_.push_back(el.ecalIso());
        el_hcalIso_.push_back(el.hcalIso());
        el_trackIso_.push_back(el.trackIso());

        el_r9_.push_back(el.r9());
        el_sMin_.push_back(el.sMin());
        el_sMaj_.push_back(el.sMaj());
        el_seedId_.push_back(el.seedId());
        el_nClusters_.push_back(el.nClusters());
        el_nCrystals_.push_back(el.nCrystals());
      }
    }
  }

  // ----- photons -----
  if (hPho.isValid()) {
    for (const auto& ph : *hPho) {
      if (ph.pt() > 0) {
        ph_pt_.push_back(ph.pt());
        ph_eta_.push_back(ph.eta());
        ph_phi_.push_back(ph.phi());
        ph_m_.push_back(ph.m());

        ph_rawE_.push_back(ph.rawEnergy());
        ph_preshowerE_.push_back(ph.preshowerEnergy());
        ph_corrEcalEnergyError_.push_back(ph.corrEcalEnergyError());

        ph_sigmaIetaIeta_.push_back(ph.sigmaIetaIeta());
        ph_hOverE_.push_back(ph.hOverE());
        ph_ecalIso_.push_back(ph.ecalIso());
        ph_hcalIso_.push_back(ph.hcalIso());
        ph_trkIso_.push_back(ph.trkIso());

        ph_r9_.push_back(ph.r9());
        ph_sMin_.push_back(ph.sMin());
        ph_sMaj_.push_back(ph.sMaj());
        ph_seedId_.push_back(ph.seedId());
        ph_nClusters_.push_back(ph.nClusters());
        ph_nCrystals_.push_back(ph.nCrystals());
      }
    }
  }

  // ----- PF jets -----
  if (hJet.isValid()) {
    for (const auto& jet : *hJet) {
      if (jet.pt() > 0) {
        jet_pt_.push_back(jet.pt());
        jet_eta_.push_back(jet.eta());
        jet_phi_.push_back(jet.phi());
        jet_m_.push_back(jet.m());
        jet_area_.push_back(jet.jetArea());

        jet_chHadE_.push_back(jet.chargedHadronEnergy());
        jet_neuHadE_.push_back(jet.neutralHadronEnergy());
        jet_photonE_.push_back(jet.photonEnergy());
        jet_electronE_.push_back(jet.electronEnergy());
        jet_muonE_.push_back(jet.muonEnergy());
        jet_HFHadE_.push_back(jet.HFHadronEnergy());
        jet_HFEME_.push_back(jet.HFEMEnergy());
        jet_HOEnergy_.push_back(jet.HOEnergy());

        jet_chHadMult_.push_back(jet.chargedHadronMultiplicity());
        jet_neuHadMult_.push_back(jet.neutralHadronMultiplicity());
        jet_photonMult_.push_back(jet.photonMultiplicity());
        jet_electronMult_.push_back(jet.electronMultiplicity());
        jet_muonMult_.push_back(jet.muonMultiplicity());
        jet_HFHadMult_.push_back(jet.HFHadronMultiplicity());
        jet_HFEMMult_.push_back(jet.HFEMMultiplicity());

        jet_csv_.push_back(jet.csv());
        jet_mvaDisc_.push_back(jet.mvaDiscriminator());
      }
    }
  }

  // ----- PF particles -----
  if (hPFC.isValid()) {
    for (const auto& p : *hPFC) {
      if (p.pt() > 0) {
        pfc_pt_.push_back(p.pt());
        pfc_eta_.push_back(p.eta());
        pfc_phi_.push_back(p.phi());
        pfc_pdgId_.push_back(p.pdgId());

        pfc_normchi2_.push_back(p.normchi2());
        pfc_dz_.push_back(p.dz());
        pfc_dxy_.push_back(p.dxy());
        pfc_dzsig_.push_back(p.dzsig());
        pfc_dxysig_.push_back(p.dxysig());
        pfc_lostInnerHits_.push_back(p.lostInnerHits());
        pfc_quality_.push_back(p.quality());
      }
    }
  }

  // ----- primary vertices -----
  if (hPV.isValid()) {
    for (const auto& v : *hPV) {
      pv_x_.push_back(v.x());
      pv_y_.push_back(v.y());
      pv_z_.push_back(v.z());
      pv_zError_.push_back(v.zError());
      pv_xError_.push_back(v.xError());
      pv_yError_.push_back(v.yError());
      pv_tracksSize_.push_back(v.tracksSize());
      pv_chi2_.push_back(v.chi2());
      pv_ndof_.push_back(v.ndof());
      pv_isValid_.push_back(v.isValidVtx());
      pv_xyCov_.push_back(v.xyCov());
      pv_xzCov_.push_back(v.xzCov());
      pv_yzCov_.push_back(v.yzCov());
    }
  }

  // ----- displaced vertices (NoVtx) -----
  if (hDVNoVtx.isValid()) {
    for (const auto& v : *hDVNoVtx) {
      dvNoVtx_x_.push_back(v.x());
      dvNoVtx_y_.push_back(v.y());
      dvNoVtx_z_.push_back(v.z());
      dvNoVtx_zError_.push_back(v.zError());
      dvNoVtx_xError_.push_back(v.xError());
      dvNoVtx_yError_.push_back(v.yError());
      dvNoVtx_tracksSize_.push_back(v.tracksSize());
      dvNoVtx_chi2_.push_back(v.chi2());
      dvNoVtx_ndof_.push_back(v.ndof());
      dvNoVtx_isValid_.push_back(v.isValidVtx());
      dvNoVtx_xyCov_.push_back(v.xyCov());
      dvNoVtx_xzCov_.push_back(v.xzCov());
      dvNoVtx_yzCov_.push_back(v.yzCov());
    }
  }

  // ----- displaced vertices (Vtx) -----
  if (hDVVtx.isValid()) {
    for (const auto& v : *hDVVtx) {
      dvVtx_x_.push_back(v.x());
      dvVtx_y_.push_back(v.y());
      dvVtx_z_.push_back(v.z());
      dvVtx_zError_.push_back(v.zError());
      dvVtx_xError_.push_back(v.xError());
      dvVtx_yError_.push_back(v.yError());
      dvVtx_tracksSize_.push_back(v.tracksSize());
      dvVtx_chi2_.push_back(v.chi2());
      dvVtx_ndof_.push_back(v.ndof());
      dvVtx_isValid_.push_back(v.isValidVtx());
      dvVtx_xyCov_.push_back(v.xyCov());
      dvVtx_xzCov_.push_back(v.xzCov());
      dvVtx_yzCov_.push_back(v.yzCov());
    }
  }

  // ----- MET & rho -----
  if (hMetPt.isValid() && hMetPhi.isValid()) {
    const double metPt  = *hMetPt;
    const double metPhi = *hMetPhi;
    met_pt_.push_back(metPt);
    met_phi_.push_back(metPhi);
  }

  if (hRho.isValid()) {
    rho_ = *hRho;
  }

  // ----- TriggerResults -----
  if (hTrig.isValid()) {
    const auto& tr = *hTrig;
    const size_t n = tr.size();
    hlt_bits_.resize(n);
    for (size_t i = 0; i < n; ++i) {
      hlt_bits_[i] = static_cast<unsigned char>(tr.accept(i));
    }
  }

  tree_->Fill();
}

DEFINE_FWK_MODULE(ScoutingNtuplizer);