#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/Muon.h"  // <-- для pat::Muon
#include "DataFormats/MuonReco/interface/Muon.h"       // <-- для reco::Muon

#include "TH1D.h"
#include "TH2D.h"

#include "TLorentzVector.h"

#include <iostream>
#include <cmath>

class ScoutingMCJpsiAnalyzer : public edm::one::EDAnalyzer<> {
public:

    void beginJob() override {
        edm::Service<TFileService> fs;
        m_1dhist_["jpsi_mass_mc"] = fs->make<TH1D>("jpsi_mass_mc", "J/psi invariant mass;m_{#mu#mu} [GeV];Events", 500, 0, 10);
        m_1dhist_["jpsi_mass_ss"] = fs->make<TH1D>("jpsi_mass_ss", "J/psi invariant mass;m_{#mu#mu} [GeV];Events", 500, 0, 10);
        m_1dhist_["jpsi_mass_os"] = fs->make<TH1D>("jpsi_mass_os", "J/psi invariant mass;m_{#mu#mu} [GeV];Events", 500, 0, 10);
        m_1dhist_["jpsi_mass_pat"] = fs->make<TH1D>("jpsi_mass_pat", "J/psi mass (PAT);m_{#mu#mu} [GeV];Events", 500, 0, 10);
        m_1dhist_["jpsi_mass_reco"] = fs->make<TH1D>("jpsi_mass_reco", "J/psi mass (RECO);m_{#mu#mu} [GeV];Events", 500, 0, 10);

        m_1dhist_["pileup_true"] = fs->make<TH1D>("pileup_true", "pileup_true", 15, 0, 150);
        m_1dhist_["pileup_reco"] = fs->make<TH1D>("pileup_reco", "pileup_reco", 15, 0, 150);
    
        // Gen muons
        m_1dhist_["gen_muons_deltaR"] = fs->make<TH1D>("gen_muons_deltaR", "gen_muons_deltaR", 1000, 0, 10);

        m_2dhist_["gen_muons_pt1_pt2"]  = fs->make<TH2D>("gen_muons_pt1_pt2",  "GEN muons; p_{T}^{1}; p_{T}^{2}", 100, 0, 50, 100, 0, 50);
        m_2dhist_["gen_muons_eta1_eta2"] = fs->make<TH2D>("gen_muons_eta1_eta2", "GEN muons; #eta_{1}; #eta_{2}", 80, -2.5, 2.5, 80, -2.5, 2.5);
        m_2dhist_["gen_muons_phi1_phi2"] = fs->make<TH2D>("gen_muons_phi1_phi2", "GEN muons; #phi_{1}; #phi_{2}", 64, -3.2, 3.2, 64, -3.2, 3.2);


        m_1dhist_["hlt_muons_deltaR"] = fs->make<TH1D>("hlt_muons_deltaR", "hlt_muons_deltaR", 1000, 0, 10);

        m_2dhist_["hlt_muons_pt1_pt2"]  = fs->make<TH2D>("hlt_muons_pt1_pt2",  "HLT muons; p_{T}^{1}; p_{T}^{2}", 100, 0, 50, 100, 0, 50);
        m_2dhist_["hlt_muons_eta1_eta2"] = fs->make<TH2D>("hlt_muons_eta1_eta2", "HLT muons; #eta_{1}; #eta_{2}", 80, -2.5, 2.5, 80, -2.5, 2.5);
        m_2dhist_["hlt_muons_phi1_phi2"] = fs->make<TH2D>("hlt_muons_phi1_phi2", "HLT muons; #phi_{1}; #phi_{2}", 64, -3.2, 3.2, 64, -3.2, 3.2);

        m_2dhist_["reconstucted_jpsi_hlt"] = fs->make<TH2D>("reconstucted_jpsi_hlt", "Jpsi vs True Pile-up", 150, 0, 150, 30, 0, 1.5);
        m_2dhist_["reconstucted_jpsi_pat"] = fs->make<TH2D>("reconstucted_jpsi_pat", "Jpsi vs True Pile-up", 150, 0, 150, 30, 0, 1.5);

    }

    void endJob() override {
        
        for (const auto& [pile, true_count] : jpsi_count_true) {
            
            int reco_count = 0;

            auto it = jpsi_count_reco.find(pile);
            if (it != jpsi_count_reco.end()) {
                reco_count = it->second;
            }

            double ratio = 0.0;
            if (true_count > 0) {
                ratio = static_cast<double>(reco_count) / true_count;
            }

            int pat_reco_count = 0;

            auto ir = jpsi_count_pat.find(pile);
            if (ir != jpsi_count_pat.end()) {
                pat_reco_count = ir->second;
            }

            double pat_ratio = 0.0;
            if (true_count > 0) {
                pat_ratio = static_cast<double>(pat_reco_count) / true_count;
            }

            std::cout << "Pileup: " << pile << ", "
                      << "True: " << true_count << ", "
                      << "PAT: " << pat_reco_count << ", "
                      << "HLT Reco: " << reco_count << ", "
                      << "PAT Ratio: " << pat_ratio << ", "
                      << "HLT Ratio: " << ratio << std::endl;

            m_2dhist_["reconstucted_jpsi_hlt"]->Fill(pile, ratio);
            m_2dhist_["reconstucted_jpsi_pat"]->Fill(pile, pat_ratio);
        }

    }

    explicit ScoutingMCJpsiAnalyzer(const edm::ParameterSet& iConfig)
        : muonsToken_(consumes<std::vector<Run3ScoutingMuon>>(edm::InputTag("hltScoutingMuonPacker"))),
        vtxToken_(consumes<std::vector<Run3ScoutingVertex>>(edm::InputTag("hltScoutingPrimaryVertexPacker","primaryVtx"))),
        genToken_(consumes<std::vector<reco::GenParticle>>(edm::InputTag("prunedGenParticles"))), 
        patMuonsToken_(consumes<std::vector<pat::Muon>>(edm::InputTag("slimmedMuons"))),
        recoMuonsToken_(consumes<std::vector<reco::Muon>>(edm::InputTag("muons"))),
        pileupToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))) {}

    void analyze(const edm::Event& iEvent, const edm::EventSetup&) override {
        edm::Handle<std::vector<Run3ScoutingMuon>> muons;
        edm::Handle<std::vector<Run3ScoutingVertex>> vertices;
        edm::Handle<std::vector<reco::GenParticle>> gens;
        edm::Handle<std::vector<reco::Muon>> recoMuons;
        
        iEvent.getByToken(muonsToken_, muons);
        iEvent.getByToken(vtxToken_, vertices);
        iEvent.getByToken(genToken_, gens);
        iEvent.getByToken(recoMuonsToken_, recoMuons);

        edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
        iEvent.getByToken(pileupToken_, puInfo);

        float pile_up_true = 0;

        if (puInfo.isValid()) {
            for (const auto& pu : *puInfo) {
                if (pu.getBunchCrossing() == 0) { // Только in-time pileup
                    float ntrue = pu.getTrueNumInteractions(); // Истинное число PU взаимодействий (сгенерированных)
                    int npu = pu.getPU_NumInteractions();      // Число PU вершин, найденных реконструкцией
                    //std::cout << "TruePU: " << ntrue << "   PU: " << npu << std::endl;
                    m_1dhist_["pileup_true"]->Fill(ntrue);
                    m_1dhist_["pileup_reco"]->Fill(npu);

                    int rounded = static_cast<int>(std::round(ntrue / 2.0)) * 2;
                    pile_up_true = rounded;
                }
            }
        }

        // 1. J/ψ in MC (GEN)
        if (gens.isValid()) {
            std::vector<const reco::GenParticle*> muons;
            const reco::GenParticle* jpsi;

            for (const auto& gp : *gens) {
                if (abs(gp.pdgId()) == 13 && gp.status() == 1) {
                    muons.push_back(&gp);
                }

                if (std::abs(gp.pdgId()) == 443) { // 443 = J/psi
                    jpsi = &gp;
                    jpsi_count_true[pile_up_true] += 1;
                }

            }
            for (size_t i = 0; i < muons.size(); ++i) {
                for (size_t j = i + 1; j < muons.size(); ++j) {

                    if ((abs(muons[i]->eta()) < 2.4) && (abs(muons[j]->eta()) < 2.4)) {

                        m_2dhist_["gen_muons_pt1_pt2"]->Fill(muons[i]->pt(), muons[j]->pt());
                        m_2dhist_["gen_muons_eta1_eta2"]->Fill(muons[i]->eta(), muons[j]->eta());
                        m_2dhist_["gen_muons_phi1_phi2"]->Fill(muons[i]->phi(), muons[j]->phi());

                        double deltaR = reco::deltaR(muons[i]->eta(), muons[i]->phi(), muons[j]->eta(), muons[j]->phi());
                        m_1dhist_["gen_muons_deltaR"]->Fill(deltaR);
                        m_1dhist_["jpsi_mass_mc"]->Fill(jpsi->mass());
                    }
                }
            }
        }


        // 2. Кандидаты J/ψ из scouting-мюонов (reco)
        if (muons.isValid()) {
            for (size_t i = 0; i < muons->size(); ++i) {
                for (size_t j = i+1; j < muons->size(); ++j) {
                    const auto& mu1 = (*muons)[i];
                    const auto& mu2 = (*muons)[j];

                    TLorentzVector v1;
                    v1.SetPtEtaPhiM(mu1.pt(), mu1.eta(), mu1.phi(), MUON_MASS);

                    TLorentzVector v2;
                    v2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), MUON_MASS);

                    TLorentzVector pair = v1 + v2;
                    double mass = pair.M();

                    if ((abs(mu1.eta()) < 2.4) && (abs(mu2.eta()) < 2.4)) {

                        if (mu1.charge() + mu2.charge() != 0) {
                            //std::cout << "[RECO] J/psi candidate ss: mass=" << mass << " pt1=" << mu1.pt() << " pt2=" << mu2.pt() << std::endl;
                            m_1dhist_["jpsi_mass_ss"]->Fill(mass);
                        } else {
                            //std::cout << "[RECO] J/psi candidate os: mass=" << mass << " pt1=" << mu1.pt() << " pt2=" << mu2.pt() << std::endl;
                            m_1dhist_["jpsi_mass_os"]->Fill(mass);
                        }

                        if (mass > 2.9 && mass < 3.3) {
                            jpsi_count_reco[pile_up_true] += 1;
                        }

                        m_2dhist_["hlt_muons_pt1_pt2"]->Fill(mu1.pt(), mu2.pt());
                        m_2dhist_["hlt_muons_eta1_eta2"]->Fill(mu1.eta(), mu2.eta());
                        m_2dhist_["hlt_muons_phi1_phi2"]->Fill(mu1.phi(), mu2.phi());

                        double deltaR = reco::deltaR(mu1.eta(), mu1.phi(), mu2.eta(), mu2.phi());
                        m_1dhist_["hlt_muons_deltaR"]->Fill(deltaR);
                    }
    
                }
            }
        }

        // PAT Muons
        edm::Handle<std::vector<pat::Muon>> patMuons;
        iEvent.getByToken(patMuonsToken_, patMuons);

        if (patMuons.isValid()) {
            for (size_t i = 0; i < patMuons->size(); ++i) {
                for (size_t j = i+1; j < patMuons->size(); ++j) {
                    const auto& mu1 = (*patMuons)[i];
                    const auto& mu2 = (*patMuons)[j];
                    TLorentzVector v1, v2;
                    v1.SetPtEtaPhiM(mu1.pt(), mu1.eta(), mu1.phi(), MUON_MASS);
                    v2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), MUON_MASS);
                    TLorentzVector pair = v1 + v2;
                    double mass = pair.M();
                    if ((abs(mu1.eta()) < 2.4) && (abs(mu2.eta()) < 2.4)) {
                        if (mu1.charge() + mu2.charge() == 0) {
                            m_1dhist_["jpsi_mass_pat"]->Fill(mass);
                            if (mass > 2.9 && mass < 3.3) {
                                jpsi_count_pat[pile_up_true] += 1;
                            }
                        }
                            
                    }
                }
            }
        }

        // RECO
        if (recoMuons.isValid()) {
            for (size_t i = 0; i < recoMuons->size(); ++i) {
                for (size_t j = i+1; j < recoMuons->size(); ++j) {
                    const auto& mu1 = (*recoMuons)[i];
                    const auto& mu2 = (*recoMuons)[j];
                    TLorentzVector v1, v2;
                    v1.SetPtEtaPhiM(mu1.pt(), mu1.eta(), mu1.phi(), MUON_MASS);
                    v2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), MUON_MASS);
                    TLorentzVector pair = v1 + v2;
                    double mass = pair.M();
                    if (mu1.charge() + mu2.charge() == 0)
                        m_1dhist_["jpsi_mass_reco"]->Fill(mass);
                }
            }
        }
    }
private:
    edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muonsToken_;
    edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> vtxToken_;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genToken_;
    edm::EDGetTokenT<std::vector<pat::Muon>> patMuonsToken_;
    edm::EDGetTokenT<std::vector<reco::Muon>> recoMuonsToken_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupToken_;

    std::map<std::string, TH1D*> m_1dhist_;
    std::map<std::string, TH2D*> m_2dhist_;

    const double MUON_MASS = 0.105658;

    // pile up vs number of reconstructed jpsi
    std::map<int, int> jpsi_count_true;
    std::map<int, int> jpsi_count_reco;
    std::map<int, int> jpsi_count_pat;
};

DEFINE_FWK_MODULE(ScoutingMCJpsiAnalyzer);