// Extended version of ScoutingMuonNtuplizer with J/ψ, Z→μμ and W→μν transverse mass analysis

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h" // For MET
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"

// Trigger
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include <numeric>
#include <unordered_set>
#include <map>

class ScoutingMuonNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit ScoutingMuonNtuplizer(const edm::ParameterSet&);
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void beginJob() override;
    void endJob() override;

private:
    bool scoutingMuonID(const Run3ScoutingMuon&) const;
    void fillLxyHistograms(const Run3ScoutingVertex&, double avgX, double avgY, double mass, int lumisection);

    edm::InputTag algInputTag_;
    edm::InputTag l1tAlgBlkInputTag_;
    edm::InputTag l1tExtBlkInputTag_;
    std::vector<std::string> l1Seeds_;
    bool doL1_;
    
    edm::EDGetTokenT<std::vector<Run3ScoutingMuon>>   muonToken_;
    edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> vertexPrimaryToken_;
    edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> vertexMuonToken_;
    edm::EDGetTokenT<double>                          pfMetPhiToken_;
    edm::EDGetTokenT<double>                          pfMetPtToken_;
    edm::InputTag extInputTag_;
    edm::EDGetToken algToken_;

    edm::Service<TFileService> fs;
    TTree* tree_ = nullptr;
    TFileDirectory histoSubDir;

    std::vector<float> muon_pt_, muon_eta_, muon_phi_, muon_m_, absIso_, relIso_;
    std::vector<int> muon_charge_;
    
    std::map<std::string, TH1D*> m_1dhist_;
    std::map<std::string, TH2D*> m_2dhist_;

    std::map<int, std::vector<float>> jpsi_all_
                                    , jpsi_ss_
                                    , jpsi_os_
                                    , jpsi_soft_iso_
                                    , jpsi_hard_iso_
                                    , jpsi_barrel_
                                    , jpsi_endcap_
                                    , jpsi_passed_
                                    , jpsi_failed_
                                    , jpsi_prompt_
                                    , jpsi_non_prompt_;

    std::map<int, std::vector<int>> jpsi_cand_multiplicity_;
    std::map<int, int> triggers_;

    std::map<int, std::vector<float>> z_counts_, z_counts_failed_;
    std::map<int, std::vector<float>> wp_counts_, wn_counts_;

    std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
    std::vector<bool> l1Result_;
};

ScoutingMuonNtuplizer::ScoutingMuonNtuplizer(const edm::ParameterSet& iConfig)
  : histoSubDir(fs->mkdir("histograms")) {
    usesResource("TFileService");
    muonToken_ = consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muonCollection"));
    vertexPrimaryToken_ = consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("primaryVertexCollection"));
    vertexMuonToken_ = consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("muonVertexCollection"));
    pfMetPhiToken_ = consumes<double>(iConfig.getParameter<edm::InputTag>("pfMetPhiValue"));
    pfMetPtToken_ = consumes<double>(iConfig.getParameter<edm::InputTag>("pfMetPtValue"));

    doL1_               = iConfig.getParameter<bool>("doL1");

    algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
    extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
    algToken_    = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
    l1Seeds_     = iConfig.getParameter<std::vector<std::string>>("l1Seeds");
    l1GtUtils_   = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);

}

void ScoutingMuonNtuplizer::beginJob() {
    
    tree_ = fs->make<TTree>("Events", "Scouting muons");
    tree_->Branch("muon_pt", &muon_pt_);
    tree_->Branch("muon_eta", &muon_eta_);
    tree_->Branch("muon_phi", &muon_phi_);
    tree_->Branch("muon_m", &muon_m_);
    tree_->Branch("muon_charge", &muon_charge_);
    tree_->Branch("absIso", &absIso_);
    tree_->Branch("relIso", &relIso_);
    tree_->Branch("l1Result", "std::vector<bool>", &l1Result_);

    m_1dhist_["lxy"] = histoSubDir.make<TH1D>("lxy", "L_{xy}", 1000, -1, 1);
    m_1dhist_["lxySig"] = histoSubDir.make<TH1D>("lxySig", "L_{xy}/#sigma", 1000, -10, 100);

    m_1dhist_["allJpsi"] = histoSubDir.make<TH1D>("J/psi", "J/psi", 50, 2.6, 3.6);
    m_1dhist_["allZ"] = histoSubDir.make<TH1D>("Z count", "Z count", 60, 60, 120);
    m_1dhist_["allW"] = histoSubDir.make<TH1D>("W count", "W count", 60, 20, 140);

    // 2D histograms
    // Mass vs pt of probe
    m_2dhist_["mass_vs_pt"] = histoSubDir.make<TH2D>( "mass_vs_pt", "Mass of Jpsi with Tag > 20 GeV"
                                                    , 50, 2.6, 3.6
                                                    , 200, 0, 100 
                                                    );

    m_2dhist_["mass_vs_pt_failed"] = histoSubDir.make<TH2D>( "mass_vs_pt_failed", "Mass of Jpsi with Tag > 20 GeV"
                                                           , 50, 2.6, 3.6
                                                           , 200, 0, 100 
                                                           );

    // Mass vs soft iso
    m_2dhist_["mass_vs_iso1"] = histoSubDir.make<TH2D>( "mass_vs_iso1", "Mass of Jpsi with mu1 reliso"
                                                      , 50, 2.6, 3.6
                                                      , 100, -5, 5
                                                      );

    // Mass vs hard iso
    m_2dhist_["mass_vs_iso2"] = histoSubDir.make<TH2D>( "mass_vs_iso2", "Mass of Jpsi with mu2 reliso"
                                                      , 50, 2.6, 3.6
                                                      , 100, -5, 5
                                                      );
                                                      
}


bool ScoutingMuonNtuplizer::scoutingMuonID(const Run3ScoutingMuon& mu) const {
    
    math::PtEtaPhiMLorentzVector particle(mu.pt(), mu.eta(), mu.phi(), 0.10566);
    double normchisq_threshold = 3.0;
    double pt_threshold = 3.0;
    double eta_threshold = 2.4;
    int layer_threshold = 4;

    if (mu.pt() > pt_threshold && fabs(mu.eta()) < eta_threshold && mu.normalizedChi2() < normchisq_threshold &&
        mu.isGlobalMuon() && mu.nTrackerLayersWithMeasurement() > layer_threshold) {
        return true;
    }
    return false;
}

void ScoutingMuonNtuplizer::fillLxyHistograms(const Run3ScoutingVertex& vtx, double avgX, double avgY, double mass, int ls) {
    
    double dx = vtx.x() - avgX;
    double dy = vtx.y() - avgY;
    double Lxy = std::sqrt(dx * dx + dy * dy);
    double LxyErr = Lxy > 0 ? std::sqrt(dx * dx * vtx.xError() * vtx.xError() + dy * dy * vtx.yError() * vtx.yError()) / Lxy : 0;
    double LxySig = LxyErr > 0 ? Lxy / LxyErr : -1;

    m_1dhist_["lxy"]->Fill(Lxy);
    m_1dhist_["lxySig"]->Fill(LxySig);

    auto& target = (LxySig < 3 ? jpsi_prompt_[ls] : jpsi_non_prompt_[ls]);
    target.push_back(mass);
}

void ScoutingMuonNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
    muon_pt_.clear(); muon_eta_.clear(); muon_phi_.clear(); muon_m_.clear(); muon_charge_.clear(); absIso_.clear(); relIso_.clear();

    if (doL1_) {
        l1Result_.clear();
        l1GtUtils_->retrieveL1(iEvent, iSetup, algToken_);
        for (unsigned iseed = 0; iseed < l1Seeds_.size(); ++iseed) {
            bool l1htbit = false;

            l1GtUtils_->getFinalDecisionByName(l1Seeds_[iseed], l1htbit);
            l1Result_.push_back(l1htbit);

            if (l1htbit) {
                std::cout << l1Seeds_[iseed] << std::endl;
            }
            
        }

        auto const& algResults = l1GtUtils_->decisionsInitial(); // vector<pair<string_view, bool>>
        for (const auto& [algoname, fired] : algResults) {
            if (fired) {
                std::cout << "L1 trigger fired: " << algoname << std::endl;
                // You can also: std::string myname = std::string(algoname);
                // ... or push_back to your own vector<string>
            }
        }
    }

    edm::Handle<std::vector<Run3ScoutingMuon>> muons;
    edm::Handle<std::vector<Run3ScoutingVertex>> primaryVertices;
    edm::Handle<std::vector<Run3ScoutingVertex>> muonVertices;
    
    edm::Handle<double> pfMetPhiH;
    edm::Handle<double> pfMetPtH;

    iEvent.getByToken(muonToken_, muons);
    iEvent.getByToken(vertexPrimaryToken_, primaryVertices);
    iEvent.getByToken(vertexMuonToken_, muonVertices);
    iEvent.getByToken(pfMetPhiToken_, pfMetPhiH);
    iEvent.getByToken(pfMetPtToken_, pfMetPtH);
    
    double metPhi = *pfMetPhiH;
    double metPt = *pfMetPtH;

    int ls = iEvent.luminosityBlock();
    triggers_[ls] += 1;

    if (muons->empty() || primaryVertices->empty() || muonVertices->empty()) return;

    double avgX = std::accumulate(primaryVertices->begin(), primaryVertices->end(), 0.0,
        [](double acc, const auto& vtx){ return acc + vtx.x(); }) / primaryVertices->size();
    double avgY = std::accumulate(primaryVertices->begin(), primaryVertices->end(), 0.0,
        [](double acc, const auto& vtx){ return acc + vtx.y(); }) / primaryVertices->size();


    int candidate_per_event = 0;
    const double MUON_MASS = 0.105658;

    // --- Z and J/ψ analysis (unchanged logic) ---
    for (size_t i = 0; i < muons->size(); ++i) {
        const auto& mu1 = muons->at(i);
        
        if (!scoutingMuonID(mu1)) continue;
        
        TLorentzVector v1;
        v1.SetPtEtaPhiM(mu1.pt(), mu1.eta(), mu1.phi(), MUON_MASS);

        for (size_t j = i + 1; j < muons->size(); ++j) {
            
            const auto& mu2 = muons->at(j);

            TLorentzVector v2;
            v2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), MUON_MASS);

            TLorentzVector pair = v1 + v2;
            double mass = pair.M();

            if (mass > 60 && mass < 120) {
                if (mu1.charge() + mu2.charge() == 0) {
                    auto& target = scoutingMuonID(mu2) ? z_counts_[ls] : z_counts_failed_[ls];
                    target.push_back(mass);

                    m_1dhist_["allZ"]->Fill(mass);
                }
            }

            // Jpsi mass requirement
            if (mass < 2.6 || mass > 3.6) continue;

            // All
            jpsi_all_[ls].push_back(mass);

            if (mu1.charge() + mu2.charge() != 0) {
                jpsi_ss_[ls].push_back(mass);
                continue;
            }

            m_1dhist_["allJpsi"]->Fill(mass);

            candidate_per_event += 1;

            // Opposite Sign
            jpsi_os_[ls].push_back(mass);

            // Isolations
            float absIso1 = mu1.trackIso() + mu1.ecalIso() + mu1.hcalIso();
            float relIso1 = absIso1 / mu1.pt();

            float absIso2 = mu2.trackIso() + mu2.ecalIso() + mu2.hcalIso();
            float relIso2 = absIso2 / mu2.pt();

            // "Hard isolation"
            if ((relIso1 < 0.2) && (relIso2 < 0.2)) {
                jpsi_hard_iso_[ls].push_back(mass);
                
            }

            // "Soft isolation"
            if ((relIso1 < 0.3) || (relIso2 < 0.3)) {
                jpsi_soft_iso_[ls].push_back(mass);
            }

            m_2dhist_["mass_vs_iso1"]->Fill(mass, relIso1);
            m_2dhist_["mass_vs_iso2"]->Fill(mass, relIso2);

            // BARREL
            if ((fabs(mu1.eta()) <= 0.9) && (fabs(mu2.eta()) <= 0.9)) {
                jpsi_barrel_[ls].push_back(mass);
            }

            // ENDCAP
            if ((fabs(mu1.eta()) > 0.9) && (fabs(mu2.eta() > 0.9))) {
                jpsi_endcap_[ls].push_back(mass);
            }        

            // tag and probe!
            auto& target = scoutingMuonID(mu2) ? jpsi_passed_[ls] : jpsi_failed_[ls];
            target.push_back(mass);

            // prompt, non-prompt
            std::unordered_set<int> v1set(mu1.vtxIndx().begin(), mu1.vtxIndx().end());
            for (int vidx : mu2.vtxIndx()) {
                if (v1set.count(vidx) && static_cast<size_t>(vidx) < muonVertices->size()) {
                    fillLxyHistograms(muonVertices->at(vidx), avgX, avgY, mass, ls);
                }
            }

            // Do 1D histograms
            // Vs pt 2 Tag. Relatively hard cut for now
            if (mu1.pt() >= 20) {
                if (scoutingMuonID(mu2)) {
                    m_2dhist_["mass_vs_pt"]->Fill(mass, mu2.pt());
                } else {
                    m_2dhist_["mass_vs_pt_failed"]->Fill(mass, mu2.pt());
                }
            }   
        }
    }

    jpsi_cand_multiplicity_[ls].push_back(candidate_per_event);

    // --- W analysis: one muon + MET ---
    for (const auto& mu : *muons) {
        if (mu.pt() > 5) { // && scoutingMuonID(mu)) {
            
            if (metPt < 20.0) continue;

            double dphi = std::abs(reco::deltaPhi(mu.phi(), metPhi));
            double mt = std::sqrt(2 * mu.pt() * metPt * (1 - std::cos(dphi)));
            
            if (mu.charge() > 0) {
                wp_counts_[ls].push_back(mt);
            } else {
                wn_counts_[ls].push_back(mt);
            }

            m_1dhist_["allW"]->Fill(mt);
        }
    }

    tree_->Fill();
}

void ScoutingMuonNtuplizer::endJob() {
    
    // fill histograms
    int minLS = std::numeric_limits<int>::max();
    int maxLS = std::numeric_limits<int>::min();

    for (const auto& [key, _] : z_counts_) {
        if (key < minLS) minLS = key;
        if (key > maxLS) maxLS = key;
    }

    const int nLumiBins = maxLS - minLS + 1;

    auto make2D = [&](const std::string& name, int dataBin, double minM, double maxM) {
        return histoSubDir.make<TH2D>(name.c_str(), name.c_str(), nLumiBins, minLS - 0.5, maxLS + 0.5, dataBin, minM, maxM);
    };

    // Make trigger histogram
    m_1dhist_["trigger"] = histoSubDir.make<TH1D>("trigger", "Number of triggers per LS", nLumiBins, minLS - 0.5, maxLS + 0.5);
    for (const auto& [ls, val] : triggers_) m_1dhist_["trigger"]->Fill(ls, val);

    m_2dhist_["multiplicity"] = make2D("multiplicity", 10, 0, 10);

    m_2dhist_["all"] = make2D("all", 50, 2.6, 3.6);
    m_2dhist_["ss"] = make2D("ss", 50, 2.6, 3.6);
    m_2dhist_["os"] = make2D("os", 50, 2.6, 3.6);
    
    m_2dhist_["softiso"] = make2D("softiso", 50, 2.6, 3.6);
    m_2dhist_["hardiso"] = make2D("hardiso", 50, 2.6, 3.6);
    
    m_2dhist_["barrel"] = make2D("barrel", 50, 2.6, 3.6);
    m_2dhist_["endcap"] = make2D("endcap", 50, 2.6, 3.6);
    
    m_2dhist_["prompt"] = make2D("prompt", 50, 2.6, 3.6);
    m_2dhist_["non_prompt"] = make2D("non_prompt", 50, 2.6, 3.6);
    
    m_2dhist_["passed"] = make2D("passed", 50, 2.6, 3.6);
    m_2dhist_["failed"] = make2D("failed", 50, 2.6, 3.6);

    m_2dhist_["z_passed"] = make2D("z_passed", 60, 60, 120);
    m_2dhist_["z_failed"] = make2D("z_failed", 60, 60, 120);

    m_2dhist_["w_p"] = make2D("w_p", 60, 20, 140);
    m_2dhist_["w_n"] = make2D("w_n", 60, 20, 140);

    for (const auto& [ls, vec] : jpsi_cand_multiplicity_) for (auto m : vec) m_2dhist_["multiplicity"]->Fill(ls, m);

    for (const auto& [ls, vec] : jpsi_all_) for (auto m : vec) m_2dhist_["all"]->Fill(ls, m);
    
    for (const auto& [ls, vec] : jpsi_ss_) for (auto m : vec) m_2dhist_["ss"]->Fill(ls, m);
    for (const auto& [ls, vec] : jpsi_os_) for (auto m : vec) m_2dhist_["os"]->Fill(ls, m);

    for (const auto& [ls, vec] : jpsi_soft_iso_) for (auto m : vec) m_2dhist_["softiso"]->Fill(ls, m);
    for (const auto& [ls, vec] : jpsi_hard_iso_) for (auto m : vec) m_2dhist_["hardiso"]->Fill(ls, m);

    for (const auto& [ls, vec] : jpsi_barrel_) for (auto m : vec) m_2dhist_["barrel"]->Fill(ls, m);
    for (const auto& [ls, vec] : jpsi_endcap_) for (auto m : vec) m_2dhist_["endcap"]->Fill(ls, m);
    
    for (const auto& [ls, vec] : jpsi_prompt_) for (auto m : vec) m_2dhist_["prompt"]->Fill(ls, m);
    for (const auto& [ls, vec] : jpsi_non_prompt_) for (auto m : vec) m_2dhist_["non_prompt"]->Fill(ls, m);
    
    for (const auto& [ls, vec] : jpsi_passed_) for (auto m : vec) m_2dhist_["passed"]->Fill(ls, m);
    for (const auto& [ls, vec] : jpsi_failed_) for (auto m : vec) m_2dhist_["failed"]->Fill(ls, m);
    
    for (const auto& [ls, vec] : z_counts_) for (auto m : vec) m_2dhist_["z_passed"]->Fill(ls, m);
    for (const auto& [ls, vec] : z_counts_failed_) for (auto m : vec) m_2dhist_["z_failed"]->Fill(ls, m);

    for (const auto& [ls, vec] : wp_counts_) for (auto m : vec) m_2dhist_["w_p"]->Fill(ls, m);
    for (const auto& [ls, vec] : wn_counts_) for (auto m : vec) m_2dhist_["w_n"]->Fill(ls, m);
}

DEFINE_FWK_MODULE(ScoutingMuonNtuplizer);