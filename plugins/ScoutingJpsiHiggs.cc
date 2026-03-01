// ScoutingJpsiHiggs — H→J/ψJ/ψ (4μ) и H→J/ψγ (μμγ) для HLT-Scouting
// Чистая реализация: димю-пары в окне J/ψ, лучшая разбивка, фотон-канал с FSR veto, контрольки (SB/SS/prompt/NP).

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"

#include "TLorentzVector.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include <cmath>
#include <numeric>
#include <algorithm>
#include <unordered_set>
#include <vector>
#include <string>
#include <map>

class ScoutingJpsiHiggs : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingJpsiHiggs(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override {}

private:
  // ===== helpers =====
  struct DiMu {
    int i{-1}, j{-1};
    TLorentzVector v;
    double m{0.0};
    double dR{0.0};
    bool os{false};
    bool jpsiLike{false};
    bool zLike{false};
  };

  static inline double deltaR(double e1, double p1, double e2, double p2) {
    double dphi = std::fabs(p1 - p2); if (dphi > M_PI) dphi = 2*M_PI - dphi;
    double deta = e1 - e2;
    return std::sqrt(deta*deta + dphi*dphi);
  }
  inline bool disjoint(const DiMu& a, const DiMu& b) const {
    return a.i!=b.i && a.i!=b.j && a.j!=b.i && a.j!=b.j;
  }

  struct Pairing { const DiMu* A{nullptr}; const DiMu* B{nullptr}; double chi2{1e99}; };
  template <typename PredA, typename PredB>
  Pairing pickBest(const std::vector<DiMu>& v, PredA acceptA, PredB acceptB,
                   double m0A, double sA, double m0B, double sB) const {
    Pairing best;
    for (size_t i = 0; i < v.size(); ++i) if (acceptA(v[i]))
    for (size_t j = i+1; j < v.size(); ++j) if (acceptB(v[j])) {
      if (!disjoint(v[i], v[j])) continue;
      const double chi2 = std::pow((v[i].m - m0A)/sA, 2) + std::pow((v[j].m - m0B)/sB, 2);
      if (chi2 < best.chi2) best = { &v[i], &v[j], chi2 };
    }
    return best;
  }

  bool scoutingMuonID(const Run3ScoutingMuon&) const;
  bool scoutingPhotonID(const Run3ScoutingPhoton&) const;

  double computeLxySigIfCommonVertex(const Run3ScoutingMuon& mu1,
                                     const Run3ScoutingMuon& mu2,
                                     const std::vector<Run3ScoutingVertex>& muonVertices,
                                     double avgX, double avgY) const;

  void sumw2(TH1* h) { if (h) h->Sumw2(); }
  void sumw2(TH2* h) { if (h) h->Sumw2(); }

  void book1D(const std::string& name, const std::string& title, int n, double lo, double hi) {
    m_1dhist_[name] = histoSubDir.make<TH1D>(name.c_str(), title.c_str(), n, lo, hi);
    sumw2(m_1dhist_[name]);
  }
  void book2D(const std::string& name, const std::string& title, int nx, double xlo, double xhi, int ny, double ylo, double yhi) {
    m_2dhist_[name] = histoSubDir.make<TH2D>(name.c_str(), title.c_str(), nx, xlo, xhi, ny, ylo, yhi);
    sumw2(m_2dhist_[name]);
  }

  // ===== cfg / tokens =====
  edm::EDGetTokenT<std::vector<Run3ScoutingMuon>>   muonToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> vertexPrimaryToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> vertexMuonToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingPhoton>> photonToken_;

  // muon/dimuon windows & thresholds
  double muonMass_{0.105658};
  double jMin_{2.95},  jMax_{3.15};     // J/psi signal window
  double hMin_{115.0}, hMax_{135.0};    // blind window for H

  double minPtForJpsiMuon_{4.0};
  double maxAbsEta_{2.4};
  double minDeltaR_{0.02};              // general ΔR(μ,μ) cut
  double minDeltaRLowM_{0.02};          // not critical for J/psi; оставлен параметром
  double sigmaJpsi_{0.5};             // ~45 MeV (tune из данных)

  // photons
  bool   doPhotons_{true};
  double minPtPhoton_{25.0};
  double maxAbsEtaPhoton_{2.5};
  double minDRMuGamma_{0.2};            // FSR veto

  // outputs
  edm::Service<TFileService> fs;
  TFileDirectory histoSubDir;
  TTree* tree_{nullptr};

  // flat branches (минимум для быстрой проверки)
  std::vector<float> muon_pt_, muon_eta_, muon_phi_;
  std::vector<int>   muon_charge_;

  std::map<std::string, TH1D*> m_1dhist_;
  std::map<std::string, TH2D*> m_2dhist_;
};

// ===== impl =====

ScoutingJpsiHiggs::ScoutingJpsiHiggs(const edm::ParameterSet& iConfig)
  : histoSubDir(fs->mkdir("histograms"))
{
  usesResource("TFileService");

  muonToken_          = consumes<std::vector<Run3ScoutingMuon>>(   iConfig.getParameter<edm::InputTag>("muonCollection"));
  vertexPrimaryToken_ = consumes<std::vector<Run3ScoutingVertex>>( iConfig.getParameter<edm::InputTag>("primaryVertexCollection"));
  vertexMuonToken_    = consumes<std::vector<Run3ScoutingVertex>>( iConfig.getParameter<edm::InputTag>("muonVertexCollection"));

  doPhotons_          = iConfig.getUntrackedParameter<bool>("doPhotons", true);
  if (doPhotons_ && iConfig.existsAs<edm::InputTag>("photonCollection"))
    photonToken_      = consumes<std::vector<Run3ScoutingPhoton>>(iConfig.getParameter<edm::InputTag>("photonCollection"));

  auto get = [&](const char* k, double& v){ if (iConfig.existsAs<double>(k)) v = iConfig.getParameter<double>(k); };
  get("jMin", jMin_); get("jMax", jMax_);
  get("hMin", hMin_); get("hMax", hMax_);
  get("minPtForJpsiMuon", minPtForJpsiMuon_);
  get("maxAbsEta", maxAbsEta_);
  get("minDeltaR", minDeltaR_);
  get("minDeltaRLowM", minDeltaRLowM_);
  get("sigmaJpsi", sigmaJpsi_);
  get("minPtPhoton", minPtPhoton_);
  get("maxAbsEtaPhoton", maxAbsEtaPhoton_);
  get("minDRMuGamma", minDRMuGamma_);
}

void ScoutingJpsiHiggs::beginJob() {
  // flat event tree (для быстрых sanity-чеков)
  tree_ = fs->make<TTree>("Events", "Scouting Jpsi Higgs");
  tree_->Branch("muon_pt",   &muon_pt_);
  tree_->Branch("muon_eta",  &muon_eta_);
  tree_->Branch("muon_phi",  &muon_phi_);
  tree_->Branch("muon_charge", &muon_charge_);

  // ---- Histograms ----
  // Inclusive dimuon & J/psi QA
  book1D("h_dimu_all", "m_{#mu#mu} [GeV]", 240, 0, 12);
  book1D("h_Jpsi_os",  "OS J/#psi-like m_{#mu#mu} [GeV]", 60, jMin_-0.15, jMax_+0.15);
  book1D("h_Z_os",  "OS Z-like m_{#mu#mu} [GeV]", 80, 60, 140);
  book1D("h_Jpsi_ss",  "SS J/#psi-like m_{#mu#mu} [GeV]", 60, jMin_-0.15, jMax_+0.15);
  book1D("h_LxySig", "L_{xy}/#sigma (J/#psi cand)", 500, -5, 45);
  book1D("h_dR_mumu_J", "#Delta R_{#mu#mu}(J/#psi)", 500, 0, 5);
  book1D("h_ptJ", "p_{T}(J/#psi) [GeV]", 100, 0, 200);
  book1D("h_ptJ_prompt", "p_{T}(J/#psi) prompt [GeV]", 100, 0, 200);
  book1D("h_ptJ_nonprompt", "p_{T}(J/#psi) non-prompt [GeV]", 100, 0, 200);

  // H → J/ψ J/ψ
  book1D("JJ_all",          "m_{J/#psi+J/#psi} [GeV]", 200, 0, 200);
  book1D("h_mH_JJ",          "m_{J/#psi+J/#psi} [GeV]", 80, 60, 220);
  book1D("h_mH_JJ_prompt",   "m_{J/#psi+J/#psi} prompt [GeV]", 80, 60, 220);
  book1D("h_mH_JJ_nonprompt","m_{J/#psi+J/#psi} non-prompt [GeV]", 80, 60, 220);
  book2D("h2_mJ1_vs_mJ2", "m_{J/#psi1} vs m_{J/#psi2}; m_{1} [GeV]; m_{2} [GeV]", 50, jMin_-0.15, jMax_+0.15, 50, jMin_-0.15, jMax_+0.15);
  book2D("h2_ptJ1_vs_ptJ2", "pt_{J/#psi1} vs pt_{J/#psi2}; pt_{1} [GeV]; pt_{2} [GeV]", 100, 0, 200, 100, 0, 200);

  // H → J/ψ γ
  if (doPhotons_) {
    book1D("h_ptGamma", "p_{T}(#gamma) [GeV]", 80, 0, 160);
    book1D("h_dr_mu_gamma_min", "min #DeltaR(#mu,#gamma)", 60, 0, 0.6);
    book1D("h_dphi_Vgamma", "#Delta#phi(J/#psi,#gamma)", 64, 0, M_PI);
    book1D("h_m_Jpsigamma", "m_{J/#psi+#gamma} [GeV]", 80, 60, 220);
    book1D("h_m_Jpsigamma_high", "m_{J/#psi+#gamma} [GeV]", 40, 110, 130);
    book1D("h_m_Jpsigamma_prompt", "m_{J/#psi+#gamma} prompt [GeV]", 80, 60, 220);
    book1D("h_m_Jpsigamma_nonprompt", "m_{J/#psi+#gamma} non-prompt [GeV]", 80, 60, 220);
    book2D("h2_mJ_vs_mH", "m_{J/#psi} vs m_{J/#psi#gamma}; m_{J/#psi} [GeV]; m_{J/#psi#gamma} [GeV]", 60, jMin_-0.15, jMax_+0.15, 120, 60, 160);

    // H → Z γ
    book1D("h_m_Zgamma", "m_{Z+#gamma} [GeV]", 80, 60, 220);
    book1D("h_m_Zgamma_high", "m_{Z+#gamma} [GeV]", 40, 110, 130);
    book2D("h2_mZ_vs_mH", "m_{Z} vs m_{Z#gamma}; m_{Z} [GeV]; m_{J/#psi#gamma} [GeV]", 80, 60, 140, 120, 60, 160);
  }
}

bool ScoutingJpsiHiggs::scoutingMuonID(const Run3ScoutingMuon& mu) const {
  const double normchisq_threshold = 3.0;
  const int    layer_threshold = 4;
  if (mu.pt() < minPtForJpsiMuon_) return false;
  if (std::abs(mu.eta()) >= maxAbsEta_) return false;
  if (!mu.isGlobalMuon()) return false;
  if (mu.normalizedChi2() >= normchisq_threshold) return false;
  if (mu.nTrackerLayersWithMeasurement() <= layer_threshold) return false;
  return true;
}

bool ScoutingJpsiHiggs::scoutingPhotonID(const Run3ScoutingPhoton& g) const {
  if (g.pt() < minPtPhoton_) return false;
  if (std::abs(g.eta()) > maxAbsEtaPhoton_) return false;
  return true;
}

double ScoutingJpsiHiggs::computeLxySigIfCommonVertex(
    const Run3ScoutingMuon& mu1,
    const Run3ScoutingMuon& mu2,
    const std::vector<Run3ScoutingVertex>& muonVertices,
    double avgX, double avgY) const {

  std::unordered_set<int> idx(mu1.vtxIndx().begin(), mu1.vtxIndx().end());
  for (int vidx : mu2.vtxIndx()) {
    if (!idx.count(vidx) || vidx < 0 || static_cast<size_t>(vidx) >= muonVertices.size()) continue;
    const auto& v = muonVertices[vidx];
    const double dx = v.x() - avgX, dy = v.y() - avgY;
    const double Lxy = std::sqrt(dx*dx + dy*dy);
    const double err2 = dx*dx*v.xError()*v.xError() + dy*dy*v.yError()*v.yError();
    if (Lxy <= 0.0 || err2 <= 0.0) return -1.0;
    return (Lxy*Lxy) / std::sqrt(err2);
  }
  return -1.0;
}

void ScoutingJpsiHiggs::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // reset flat
  muon_pt_.clear(); muon_eta_.clear(); muon_phi_.clear(); muon_charge_.clear();

  // inputs
  edm::Handle<std::vector<Run3ScoutingMuon>>   muons;
  edm::Handle<std::vector<Run3ScoutingVertex>> primaryVertices, muonVertices;
  edm::Handle<std::vector<Run3ScoutingPhoton>> photonsH;

  iEvent.getByToken(muonToken_,          muons);
  iEvent.getByToken(vertexPrimaryToken_, primaryVertices);
  iEvent.getByToken(vertexMuonToken_,    muonVertices);
  if (doPhotons_ && photonToken_.isUninitialized()==false)
    iEvent.getByToken(photonToken_,      photonsH);

  if (!muons || muons->empty() || !primaryVertices || primaryVertices->empty() || !muonVertices) return;

  // PV averages (для LxySig)
  const double avgX = std::accumulate(primaryVertices->begin(), primaryVertices->end(), 0.0,
                        [](double a, const auto& v){ return a + v.x(); }) / primaryVertices->size();
  const double avgY = std::accumulate(primaryVertices->begin(), primaryVertices->end(), 0.0,
                        [](double a, const auto& v){ return a + v.y(); }) / primaryVertices->size();

  // flat muon arrays
  muon_pt_.reserve(muons->size()); muon_eta_.reserve(muons->size());
  muon_phi_.reserve(muons->size()); muon_charge_.reserve(muons->size());
  for (const auto& mu : *muons) {
    muon_pt_.push_back(mu.pt());
    muon_eta_.push_back(mu.eta());
    muon_phi_.push_back(mu.phi());
    muon_charge_.push_back(mu.charge());
  }

  // ===== Build dimuons, tag J/psi-like =====
  std::vector<DiMu> dimu; dimu.reserve(muons->size()*2);

  for (size_t i = 0; i < muons->size(); ++i) {
    const auto& mu1 = muons->at(i);
    if (!scoutingMuonID(mu1)) continue;

    TLorentzVector v1; v1.SetPtEtaPhiM(mu1.pt(), mu1.eta(), mu1.phi(), muonMass_);

    for (size_t j = i+1; j < muons->size(); ++j) {
      const auto& mu2 = muons->at(j);
      if (!scoutingMuonID(mu2)) continue;

      const double dR = deltaR(mu1.eta(), mu1.phi(), mu2.eta(), mu2.phi());
      const double dRcut = (/*low-mass logic*/ false ? minDeltaRLowM_ : minDeltaR_);
      if (dR < dRcut) continue;

      TLorentzVector v2; v2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), muonMass_);
      const TLorentzVector vp = v1 + v2;
      const double m  = vp.M();
      const bool   os = (mu1.charge() + mu2.charge() == 0);

      m_1dhist_["h_dimu_all"]->Fill(m);

      // SS monitoring
      if (!os && m>jMin_ && m<jMax_)    m_1dhist_["h_Jpsi_ss"]->Fill(m);

      const bool jLike = os && (m>jMin_ && m<jMax_);
      const bool zLike = os && (m>60 && m<140);
      if (jLike) {
        m_1dhist_["h_Jpsi_os"]->Fill(m);
        m_1dhist_["h_dR_mumu_J"]->Fill(dR);
        m_1dhist_["h_ptJ"]->Fill(vp.Pt());
        // LxySig for prompt/non-prompt study
        double lxySig = computeLxySigIfCommonVertex(mu1, mu2, *muonVertices, avgX, avgY);
        if (std::isfinite(lxySig)) {
          m_1dhist_["h_LxySig"]->Fill(lxySig);
          if (lxySig <= 3.0)  m_1dhist_["h_ptJ_prompt"]->Fill(vp.Pt());
          if (lxySig > 3.0)  m_1dhist_["h_ptJ_nonprompt"]->Fill(vp.Pt());
        }
      }
      
      if (zLike) {
        m_1dhist_["h_Z_os"]->Fill(m);
      }

      DiMu d; d.i = int(i); d.j = int(j); d.v = vp; d.m = m; d.os = os; d.dR = dR; d.jpsiLike = jLike; d.zLike = zLike;
      dimu.push_back(std::move(d));
    }
  }

  // ===== H → J/ψ J/ψ (best pairing) =====
  {
    auto accJ = [&](const DiMu& d){ return d.jpsiLike; };
    const auto best = pickBest(dimu, accJ, accJ, 3.0969, sigmaJpsi_, 3.0969, sigmaJpsi_);
    if (best.A && best.B) {
      TLorentzVector vH = best.A->v + best.B->v;

      const double mH = vH.M();
      m_1dhist_["JJ_all"]->Fill(mH);
      m_1dhist_["h_mH_JJ"]->Fill(mH);
      m_2dhist_["h2_mJ1_vs_mJ2"]->Fill(best.A->m, best.B->m);
      m_2dhist_["h2_ptJ1_vs_ptJ2"]->Fill(best.A->v.Pt(), best.B->v.Pt());

      // Prompt/non-prompt categorization по двум парам
      // NB: считаем повторно LxySig (дёшево)
      const auto &muA1 = muons->at(best.A->i), &muA2 = muons->at(best.A->j);
      const auto &muB1 = muons->at(best.B->i), &muB2 = muons->at(best.B->j);
      const double lxyA = computeLxySigIfCommonVertex(muA1, muA2, *muonVertices, avgX, avgY);
      const double lxyB = computeLxySigIfCommonVertex(muB1, muB2, *muonVertices, avgX, avgY);
      const bool isPromptPair = (lxyA<=3.0 && lxyB<=3.0);
      const bool isNonPrompt  = (lxyA>3.0 || lxyB>3.0);
      if (isPromptPair) m_1dhist_["h_mH_JJ_prompt"]->Fill(mH);
      if (isNonPrompt)  m_1dhist_["h_mH_JJ_nonprompt"]->Fill(mH);
    }
  }

  // ===== H → J/ψ γ =====
  if (doPhotons_ && photonsH && !photonsH->empty()) {
    // собрать годные фотоны
    std::vector<Run3ScoutingPhoton> photons; photons.reserve(photonsH->size());
    for (const auto& g : *photonsH) if (scoutingPhotonID(g)) photons.push_back(g);

    if (!photons.empty()) {
      for (const auto& J : dimu) {
        if (!J.jpsiLike && !J.zLike) continue;
        const TLorentzVector vJ = J.v;

        for (const auto& g : photons) {
          // FSR veto: min ΔR(γ, μ из J/ψ)
          const auto& mu1 = muons->at(J.i);
          const auto& mu2 = muons->at(J.j);
          const double dR1 = deltaR(mu1.eta(), mu1.phi(), g.eta(), g.phi());
          const double dR2 = deltaR(mu2.eta(), mu2.phi(), g.eta(), g.phi());
          const double dRmin = std::min(dR1, dR2);
          m_1dhist_["h_dr_mu_gamma_min"]->Fill(dRmin);
          if (dRmin < minDRMuGamma_) continue;

          TLorentzVector vg; vg.SetPtEtaPhiM(g.pt(), g.eta(), g.phi(), 0.0);
          const TLorentzVector vH = vJ + vg;
          const double mH = vH.M();

          if (J.jpsiLike) {
            m_1dhist_["h_ptGamma"]->Fill(vg.Pt());
            const double dphi = std::fabs(std::atan2(std::sin(vJ.Phi()-vg.Phi()), std::cos(vJ.Phi()-vg.Phi())));
            m_1dhist_["h_dphi_Vgamma"]->Fill(dphi);
            m_1dhist_["h_m_Jpsigamma_high"]->Fill(mH);
            m_1dhist_["h_m_Jpsigamma"]->Fill(mH);
            m_2dhist_["h2_mJ_vs_mH"]->Fill(J.m, mH);

            // Prompt/non-prompt split по LxySig(J/ψ)
            const double lxy = computeLxySigIfCommonVertex(mu1, mu2, *muonVertices, avgX, avgY);
            if (std::isfinite(lxy)) {
              if (lxy < 3.0) m_1dhist_["h_m_Jpsigamma_prompt"]->Fill(mH);
              else if (lxy > 5.0) m_1dhist_["h_m_Jpsigamma_nonprompt"]->Fill(mH);
            }
          }

          if (J.zLike) {
            m_1dhist_["h_m_Zgamma_high"]->Fill(mH);
            m_1dhist_["h_m_Zgamma"]->Fill(mH);
            m_2dhist_["h2_mZ_vs_mH"]->Fill(J.m, mH);
          }
        }
      }
    }
  }

  //tree_->Fill();
}

DEFINE_FWK_MODULE(ScoutingJpsiHiggs);