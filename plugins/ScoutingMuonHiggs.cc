// ScoutingMuonHiggs — streamlined selection & histos for:
//   H→ZZ→4μ, H→Z(μμ)+J/ψ(μμ), H→J/ψ+J/ψ  (+ control regions)

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
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// L1T
//#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
//#include "DataFormats/L1Trigger/interface/BXVector.h"
//#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

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

class ScoutingMuonHiggs : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingMuonHiggs(const edm::ParameterSet&);
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override {}

private:
  struct DiMu {
    int i{-1}, j{-1};
    TLorentzVector v;
    double m{0.0};
    double dR{0.0};
    bool os{false};
    bool zLike{false};
    bool jpsiLike{false};
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
  Pairing pickBest(const std::vector<DiMu>& v, PredA acceptA, PredB acceptB, double m0A, double sA, double m0B, double sB) const {
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

  // cfg / tokens
  bool doL1_{false}, doW_{false}, doWrongPairing_{true};
  //edm::InputTag algInputTag_, extInputTag_;
  //std::vector<std::string> l1Seeds_;
  //edm::EDGetTokenT<l1t::BXVector<l1t::GlobalAlgBlk>> algToken_;

  edm::EDGetTokenT<std::vector<Run3ScoutingMuon>>   muonToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> vertexPrimaryToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> vertexMuonToken_;
  edm::EDGetTokenT<double>                          pfMetPhiToken_;
  edm::EDGetTokenT<double>                          pfMetPtToken_;

  // windows / thresholds
  double muonMass_{0.105658};
  double zMin_{60.0}, zMax_{120.0};
  double zCoreMin_{80.0}, zCoreMax_{100.0};
  double jMin_{2.6},  jMax_{3.6};
  double jSbLoMax_{2.90}, jSbHiMin_{3.30};
  double zSbLoMin_{76.0}, zSbLoMax_{86.0}, zSbHiMin_{96.0}, zSbHiMax_{106.0};
  double hMin_{115.0}, hMax_{135.0};

  double minPtForJpsiMuon_{3.0};
  double minPtForZMuon_{10.0};
  double maxAbsEta_{2.4};
  double minDeltaR_{0.02};

  double sigmaZ_{2.5};
  double sigmaJpsi_{0.045};

  double promptLxySigMax_{3.0};
  double nonpromptLxySigMin_{5.0};

  // outputs
  edm::Service<TFileService> fs;
  TFileDirectory histoSubDir;
  TTree* tree_{nullptr};
  std::vector<float> muon_pt_, muon_eta_, muon_phi_, muon_m_, absIso_, relIso_;
  std::vector<int>   muon_charge_;
  std::vector<bool>  l1Result_;

  std::map<std::string, TH1D*> m_1dhist_;
  std::map<std::string, TH2D*> m_2dhist_;

  //std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
};

ScoutingMuonHiggs::ScoutingMuonHiggs(const edm::ParameterSet& iConfig)
  : histoSubDir(fs->mkdir("histograms"))
{
  usesResource("TFileService");

  muonToken_          = consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muonCollection"));
  vertexPrimaryToken_ = consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("primaryVertexCollection"));
  vertexMuonToken_    = consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("muonVertexCollection"));
  pfMetPhiToken_      = consumes<double>(iConfig.getParameter<edm::InputTag>("pfMetPhiValue"));
  pfMetPtToken_       = consumes<double>(iConfig.getParameter<edm::InputTag>("pfMetPtValue"));

  doL1_        = iConfig.getParameter<bool>("doL1");
  doW_         = iConfig.getUntrackedParameter<bool>("doW", false);
  doWrongPairing_ = iConfig.getUntrackedParameter<bool>("doWrongPairing", true);

  //algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
  //extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
  //algToken_    = consumes<l1t::BXVector<l1t::GlobalAlgBlk>>(algInputTag_);
  //l1Seeds_     = iConfig.getParameter<std::vector<std::string>>("l1Seeds");
  //l1GtUtils_   = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);

  auto get = [&](const char* k, double& v){ if (iConfig.existsAs<double>(k)) v = iConfig.getParameter<double>(k); };
  get("zMin", zMin_); get("zMax", zMax_); get("zCoreMin", zCoreMin_); get("zCoreMax", zCoreMax_);
  get("jMin", jMin_); get("jMax", jMax_); get("jSbLoMax", jSbLoMax_); get("jSbHiMin", jSbHiMin_);
  get("zSbLoMin", zSbLoMin_); get("zSbLoMax", zSbLoMax_); get("zSbHiMin", zSbHiMin_); get("zSbHiMax", zSbHiMax_);
  get("hMin", hMin_); get("hMax", hMax_);
  get("minPtForJpsiMuon", minPtForJpsiMuon_); get("minPtForZMuon", minPtForZMuon_);
  get("maxAbsEta", maxAbsEta_); get("minDeltaR", minDeltaR_);
  get("sigmaZ", sigmaZ_); get("sigmaJpsi", sigmaJpsi_);
  get("promptLxySigMax", promptLxySigMax_); get("nonpromptLxySigMin", nonpromptLxySigMin_);
}

void ScoutingMuonHiggs::beginJob() {
  // flat event tree
  tree_ = fs->make<TTree>("Events", "Scouting muons (flat)");
  tree_->Branch("muon_pt",     &muon_pt_);
  tree_->Branch("muon_eta",    &muon_eta_);
  tree_->Branch("muon_phi",    &muon_phi_);
  tree_->Branch("muon_m",      &muon_m_);
  tree_->Branch("muon_charge", &muon_charge_);
  tree_->Branch("absIso",      &absIso_);
  tree_->Branch("relIso",      &relIso_);
  tree_->Branch("l1Result",    "std::vector<bool>", &l1Result_);

  // histos
  book1D("h_dimu_all", "m_{#mu#mu} [GeV]", 20000, 0, 200);
  book1D("h_Z_os",     "OS Z-like m_{#mu#mu} [GeV]", 60, zMin_, zMax_);
  book1D("h_Jpsi_os",  "OS J/#psi-like m_{#mu#mu} [GeV]", 50, jMin_, jMax_);
  book1D("h_Jpsi_ss",  "SS J/#psi-like m_{#mu#mu} [GeV]", 50, jMin_, jMax_);
  book1D("h_Jpsi_sb_lo", "J/#psi sideband low; m_{#mu#mu} [GeV]", 60, 2.6, jSbLoMax_);
  book1D("h_Jpsi_sb_hi", "J/#psi sideband high; m_{#mu#mu} [GeV]", 60, jSbHiMin_, 3.6);
  book1D("h_Z_sb_lo", "Z sideband low; m_{#mu#mu} [GeV]", 40, zSbLoMin_, zSbLoMax_);
  book1D("h_Z_sb_hi", "Z sideband high; m_{#mu#mu} [GeV]", 40, zSbHiMin_, zSbHiMax_);
  book1D("h_LxySig", "L_{xy}/#sigma (J/#psi cand)", 80, -10, 30);
  book1D("h_m4mu",   "m_{4#mu} [GeV]", 160, 60, 220);
  book1D("h_mH_ZJ",  "m_{Z+J/#psi} [GeV]", 160, 60, 220);
  book1D("h_mH_JJ",  "m_{J/#psi+J/#psi} [GeV]", 160, 60, 220);
  if (doWrongPairing_) book1D("h_m4mu_wrongPair", "m_{4#mu} (wrong pairing) [GeV]", 160, 60, 220);
  book1D("h_pt4", "p_{T}(4-body) [GeV]", 100, 0, 200);
  book1D("h_ptZ", "p_{T}(Z) [GeV]",      100, 0, 200);
  book1D("h_ptJ", "p_{T}(J/#psi) [GeV]", 100, 0, 100);
  book1D("h_dR_mumu", "#Delta R(#mu,#mu)", 50, 0, 5);
  book2D("h2_mZ_vs_mJ",     "m_{Z} vs m_{J/#psi}; m_{Z} [GeV]; m_{J/#psi} [GeV]", 60, zMin_, zMax_, 50, jMin_, jMax_);
  book2D("h2_mJ1_vs_mJ2",   "m_{J/#psi1} vs m_{J/#psi2}; m_{1} [GeV]; m_{2} [GeV]", 50, jMin_, jMax_, 50, jMin_, jMax_);
  book2D("h2_LxySig_vs_mJ", "L_{xy}/#sigma vs m_{J/#psi}; m_{J/#psi} [GeV]; L_{xy}/#sigma", 50, jMin_, jMax_, 60, -5, 25);
}

bool ScoutingMuonHiggs::scoutingMuonID(const Run3ScoutingMuon& mu) const {
  const double normchisq_threshold = 3.0;
  const int    layer_threshold = 4;
  if (mu.pt() < 3.0) return false;
  if (std::abs(mu.eta()) >= maxAbsEta_) return false;
  if (!mu.isGlobalMuon()) return false;
  if (mu.normalizedChi2() >= normchisq_threshold) return false;
  if (mu.nTrackerLayersWithMeasurement() <= layer_threshold) return false;
  return true;
}

double ScoutingMuonHiggs::computeLxySigIfCommonVertex(const Run3ScoutingMuon& mu1,
                                                      const Run3ScoutingMuon& mu2,
                                                      const std::vector<Run3ScoutingVertex>& muonVertices,
                                                      double avgX, double avgY) const {
  std::unordered_set<int> set1(mu1.vtxIndx().begin(), mu1.vtxIndx().end());
  for (int vidx : mu2.vtxIndx()) {
    if (set1.count(vidx) && static_cast<size_t>(vidx) < muonVertices.size()) {
      const auto& v = muonVertices[vidx];
      const double dx = v.x() - avgX;
      const double dy = v.y() - avgY;
      const double Lxy = std::sqrt(dx*dx + dy*dy);
      const double err2 = dx*dx*v.xError()*v.xError() + dy*dy*v.yError()*v.yError();
      if (Lxy > 0.0 && err2 > 0.0) return Lxy / std::sqrt(err2);
      return -1.0;
    }
  }
  return -1.0;
}

void ScoutingMuonHiggs::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  muon_pt_.clear(); muon_eta_.clear(); muon_phi_.clear(); muon_m_.clear();
  muon_charge_.clear(); absIso_.clear(); relIso_.clear();

  if (doL1_) {
    /*
    l1Result_.clear();
    l1GtUtils_->retrieveL1(iEvent, iSetup, algToken_);
    for (const auto& name : l1Seeds_) {
      bool fired = false;
      l1GtUtils_->getFinalDecisionByName(name, fired);
      l1Result_.push_back(fired);
      if (fired) edm::LogInfo("ScoutingMuonHiggs") << "L1 fired: " << name;
    }
    */
  }

  edm::Handle<std::vector<Run3ScoutingMuon>> muons;
  edm::Handle<std::vector<Run3ScoutingVertex>> primaryVertices, muonVertices;
  edm::Handle<double> pfMetPhiH, pfMetPtH;

  iEvent.getByToken(muonToken_,          muons);
  iEvent.getByToken(vertexPrimaryToken_, primaryVertices);
  iEvent.getByToken(vertexMuonToken_,    muonVertices);
  iEvent.getByToken(pfMetPhiToken_,      pfMetPhiH);
  iEvent.getByToken(pfMetPtToken_,       pfMetPtH);

  if (!muons || muons->empty() || !primaryVertices || primaryVertices->empty() || !muonVertices) return;

  const double metPhi = (pfMetPhiH.isValid() ? *pfMetPhiH : 0.0);
  const double metPt  = (pfMetPtH.isValid()  ? *pfMetPtH  : 0.0);

  const double avgX = std::accumulate(primaryVertices->begin(), primaryVertices->end(), 0.0,
                        [](double a, const auto& v){ return a + v.x(); }) / primaryVertices->size();
  const double avgY = std::accumulate(primaryVertices->begin(), primaryVertices->end(), 0.0,
                        [](double a, const auto& v){ return a + v.y(); }) / primaryVertices->size();

  muon_pt_.reserve(muons->size()); muon_eta_.reserve(muons->size());
  muon_phi_.reserve(muons->size()); muon_m_.reserve(muons->size());
  muon_charge_.reserve(muons->size()); absIso_.reserve(muons->size()); relIso_.reserve(muons->size());
  for (const auto& mu : *muons) {
    const float aIso = mu.trackIso() + mu.ecalIso() + mu.hcalIso();
    muon_pt_.push_back(mu.pt());
    muon_eta_.push_back(mu.eta());
    muon_phi_.push_back(mu.phi());
    muon_m_.push_back(muonMass_);
    muon_charge_.push_back(mu.charge());
    absIso_.push_back(aIso);
    relIso_.push_back( (mu.pt()>0) ? aIso/mu.pt() : 999.f );
  }

  // build dimu list
  std::vector<DiMu> dimu; dimu.reserve(muons->size()*2);

  for (size_t i = 0; i < muons->size(); ++i) {
    const auto& mu1 = muons->at(i);
    if (!scoutingMuonID(mu1)) continue;

    TLorentzVector v1; v1.SetPtEtaPhiM(mu1.pt(), mu1.eta(), mu1.phi(), muonMass_);

    for (size_t j = i+1; j < muons->size(); ++j) {
      const auto& mu2 = muons->at(j);
      if (!scoutingMuonID(mu2)) continue;

      const double dR = deltaR(mu1.eta(), mu1.phi(), mu2.eta(), mu2.phi());
      if (dR < minDeltaR_) continue;

      TLorentzVector v2; v2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), muonMass_);
      const TLorentzVector vp = v1 + v2;
      const double m  = vp.M();
      const bool   os = (mu1.charge() + mu2.charge() == 0);

      m_1dhist_["h_dimu_all"]->Fill(m);
      m_1dhist_["h_dR_mumu"]->Fill(dR);

      const bool passJpsiKine = (mu1.pt()>=minPtForJpsiMuon_ && mu2.pt()>=minPtForJpsiMuon_);
      const bool passZKine    = (mu1.pt()>=minPtForZMuon_    && mu2.pt()>=minPtForZMuon_);

      const bool jLike = os && passJpsiKine && (m > jMin_ && m < jMax_);
      const bool zLike = os && passZKine    && (m > zMin_ && m < zMax_);

      if (jLike) m_1dhist_["h_Jpsi_os"]->Fill(m);
      if (!os && (m > jMin_ && m < jMax_)) m_1dhist_["h_Jpsi_ss"]->Fill(m);
      if (zLike) m_1dhist_["h_Z_os"]->Fill(m);

      if (m >= 2.6 && m <= jSbLoMax_) m_1dhist_["h_Jpsi_sb_lo"]->Fill(m);
      if (m >= jSbHiMin_ && m <= 3.6) m_1dhist_["h_Jpsi_sb_hi"]->Fill(m);
      if (m >= zSbLoMin_ && m <= zSbLoMax_) m_1dhist_["h_Z_sb_lo"]->Fill(m);
      if (m >= zSbHiMin_ && m <= zSbHiMax_) m_1dhist_["h_Z_sb_hi"]->Fill(m);

      if (jLike) {
        double lxySig = computeLxySigIfCommonVertex(mu1, mu2, *muonVertices, avgX, avgY);
        if (std::isfinite(lxySig)) m_1dhist_["h_LxySig"]->Fill(lxySig);
        m_2dhist_["h2_LxySig_vs_mJ"]->Fill(m, std::max(std::min(lxySig, 25.0), -5.0));
      }

      DiMu d; d.i = int(i); d.j = int(j); d.v = vp; d.m = m; d.os = os; d.dR = dR;
      d.zLike = zLike; d.jpsiLike = jLike;
      dimu.push_back(std::move(d));
    }
  }

  // best pairing per channel
  {
    auto accZ = [&](const DiMu& d){ return d.zLike; };
    const auto best = pickBest(dimu, accZ, accZ, 91.1876, sigmaZ_, 91.1876, sigmaZ_);
    if (best.A && best.B) {
      TLorentzVector v4 = best.A->v + best.B->v;
      m_1dhist_["h_m4mu"]->Fill(v4.M());
      m_1dhist_["h_pt4"]->Fill(v4.Pt());
      m_1dhist_["h_ptZ"]->Fill(best.A->v.Pt());
      m_1dhist_["h_ptZ"]->Fill(best.B->v.Pt());
      if (doWrongPairing_) {
        for (const auto& d : dimu) {
          if (&d==best.A || &d==best.B) continue;
          if (!d.zLike) continue;
          if (!disjoint(*best.A, d)) continue;
          TLorentzVector v4w = best.A->v + d.v;
          m_1dhist_["h_m4mu_wrongPair"]->Fill(v4w.M());
          break;
        }
      }
    }
  }

  {
    auto accZ = [&](const DiMu& d){ return d.zLike; };
    auto accJ = [&](const DiMu& d){ return d.jpsiLike; };
    const auto best = pickBest(dimu, accZ, accJ, 91.1876, sigmaZ_, 3.0969, sigmaJpsi_);
    if (best.A && best.B) {
      TLorentzVector vH = best.A->v + best.B->v;
      m_1dhist_["h_mH_ZJ"]->Fill(vH.M());
      m_1dhist_["h_pt4"]->Fill(vH.Pt());
      m_1dhist_["h_ptZ"]->Fill(best.A->v.Pt());
      m_1dhist_["h_ptJ"]->Fill(best.B->v.Pt());
      m_2dhist_["h2_mZ_vs_mJ"]->Fill(best.A->m, best.B->m);
    }
  }

  {
    auto accJ = [&](const DiMu& d){ return d.jpsiLike; };
    const auto best = pickBest(dimu, accJ, accJ, 3.0969, sigmaJpsi_, 3.0969, sigmaJpsi_);
    if (best.A && best.B) {
      TLorentzVector vH = best.A->v + best.B->v;
      m_1dhist_["h_mH_JJ"]->Fill(vH.M());
      m_1dhist_["h_pt4"]->Fill(vH.Pt());
      m_1dhist_["h_ptJ"]->Fill(best.A->v.Pt());
      m_1dhist_["h_ptJ"]->Fill(best.B->v.Pt());
      m_2dhist_["h2_mJ1_vs_mJ2"]->Fill(best.A->m, best.B->m);
    }
  }

  if (doW_ && muons->size()==1) {
    const auto& mu = muons->at(0);
    if (mu.pt()>20.0 && metPt>20.0) {
      const double dphi = std::fabs(reco::deltaPhi(mu.phi(), metPhi));
      (void)dphi; // not used; keep quiet
    }
  }

  //tree_->Fill();
}

DEFINE_FWK_MODULE(ScoutingMuonHiggs);