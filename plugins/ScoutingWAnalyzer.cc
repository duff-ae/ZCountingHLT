// ScoutingWAnalyzer.cc
// Build and validate PU-mitigated MET variants in Run-3 Scouting.
// Output: histograms vs true PU (nTrueInt), rho, nPV; compare to genMET.
//
// MET variants:
//  - met_hlt: scouting PF MET from event-level values
//  - met_chs: charged-only PF candidates associated to best PV
//  - met_puppilite: charged PV (w=1) + neutrals with PUPPI-lite weight from local charged PV vs charged PU activity

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "TH1F.h"
#include "TH2F.h"

#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <utility>

namespace {
  inline float dR2(float eta1, float phi1, float eta2, float phi2) {
    float dphi = reco::deltaPhi(phi1, phi2);
    float deta = eta1 - eta2;
    return dphi*dphi + deta*deta;
  }

  struct PtVec2 {
    float px{0.f}, py{0.f};
    inline void addPtPhi(float pt, float phi, float w=1.f) {
      px += w * pt * std::cos(phi);
      py += w * pt * std::sin(phi);
    }
    inline float pt()  const { return std::sqrt(px*px + py*py); }
    inline float phi() const { return std::atan2(py, px); }
  };

  // Very small helper: interpret PF-candidate "chargedness" from pdgId.
  // Adjust if your Run3ScoutingParticle pdgId convention differs.
  inline bool isChargedFromPdgId(int pdgId) {
    int id = std::abs(pdgId);
    // e, mu, charged hadrons, proton
    return (id == 11 || id == 13 || id == 211 || id == 321 || id == 2212);
  }
  inline bool isNeutralFromPdgId(int pdgId) {
    return !isChargedFromPdgId(pdgId);
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
    float mu_pt{0.f}, mu_eta{0.f}, mu_phi{0.f};
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
    out.mu_pt = mu->pt(); out.mu_eta = mu->eta(); out.mu_phi = mu->phi();
    out.nu_pt = nu->pt(); out.nu_eta = nu->eta(); out.nu_phi = nu->phi();
    return out;
  }

  inline float clamp01(float x) {
    if (x < 0.f) return 0.f;
    if (x > 1.f) return 1.f;
    return x;
  }
} // namespace

class ScoutingWAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingWAnalyzer(const edm::ParameterSet&);
  ~ScoutingWAnalyzer() override = default;

  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  // Config
  float fid_mu_pt_min_{25.f};
  float fid_mu_abs_eta_max_{2.4f};

  // PV selection cuts (optional)
  int   pv_min_ndof_{4};
  float pv_max_abs_z_{24.f};
  float pv_max_rho_{2.f};

  // Track/PF quality cuts (simple, tune later)
  float trk_max_normchi2_{10.f};
  float trk_max_abs_dz_{0.5f};
  float trk_max_abs_dxy_{0.2f};
  uint8_t trk_max_lostInnerHits_{2};

  // PUPPI-lite local weighting parameters
  float puppi_drmax_{0.4f};
  float puppi_eps_{0.02f};
  int   puppi_max_neighbors_{64};
  float puppi_beta_pow_{1.0f}; // use pT/(dR+eps)^pow

  // Tokens
  edm::EDGetTokenT<std::vector<Run3ScoutingMuon>>      scMuToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingVertex>>   scVtxToken_;
  edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> scPfcToken_;

  edm::EDGetTokenT<double> scMetPtToken_;
  edm::EDGetTokenT<double> scMetPhiToken_;
  edm::EDGetTokenT<double> rhoToken_;

  edm::EDGetTokenT<std::vector<PileupSummaryInfo>>    puToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>>    genToken_;
  edm::EDGetTokenT<std::vector<reco::GenMET>>         genMetToken_;

  edm::Service<TFileService> fs_;

  // Histograms
  TH1F* h_nTrueInt_{nullptr};
  TH2F* h_met_res_vs_nTrue_hlt_{nullptr};
  TH2F* h_met_res_vs_nTrue_chs_{nullptr};
  TH2F* h_met_res_vs_nTrue_puppi_{nullptr};

  TH2F* h_dphi_vs_nTrue_hlt_{nullptr};
  TH2F* h_dphi_vs_nTrue_chs_{nullptr};
  TH2F* h_dphi_vs_nTrue_puppi_{nullptr};

  TH2F* h_resp_vs_nTrue_hlt_{nullptr};
  TH2F* h_resp_vs_nTrue_chs_{nullptr};
  TH2F* h_resp_vs_nTrue_puppi_{nullptr};

  TH2F* h_met_pt_vs_nTrue_hlt_{nullptr};
  TH2F* h_met_pt_vs_nTrue_chs_{nullptr};
  TH2F* h_met_pt_vs_nTrue_puppi_{nullptr};

  TH2F* h_chgPV_vs_chgPU_{nullptr};
  TH2F* h_fracChgPU_vs_nTrue_{nullptr};
  TH1I* h_chg_vertexIdx_{nullptr};
  TH1I* h_pvIndex_{nullptr};

  TH1I* h_pfc_pdgId_{nullptr};
  TH1I* h_pfc_absPdgId_{nullptr};
  TH1F* h_pfc_pt_{nullptr};

  TH1I* h_cnt_{nullptr}; // 6-bin counter

  // Helpers
  int pickBestPV_(const std::vector<Run3ScoutingVertex>& vtxs) const;

  bool passTrkQuality_(const Run3ScoutingParticle& p) const {
    // Uses fields you showed: normchi2, dz, dxy, lostInnerHits, quality.
    if (p.normchi2() > trk_max_normchi2_) return false;
    if (std::abs(p.dz()) > trk_max_abs_dz_) return false;
    if (std::abs(p.dxy()) > trk_max_abs_dxy_) return false;
    if (p.lostInnerHits() > trk_max_lostInnerHits_) return false;
    // quality() is uint8_t; do not over-cut here until you confirm meaning.
    return true;
  }

  std::pair<float,float> buildMetCHS_(
      const std::vector<Run3ScoutingParticle>& pfc,
      int pvIndex) const;

  std::pair<float,float> buildMetPuppiLite_(
      const std::vector<Run3ScoutingParticle>& pfc,
      int pvIndex) const;
};

ScoutingWAnalyzer::ScoutingWAnalyzer(const edm::ParameterSet& cfg) {
  usesResource("TFileService");

  fid_mu_pt_min_      = cfg.getParameter<double>("fid_mu_pt_min");
  fid_mu_abs_eta_max_ = cfg.getParameter<double>("fid_mu_abs_eta_max");

  pv_min_ndof_     = cfg.getParameter<int>("pv_min_ndof");
  pv_max_abs_z_    = cfg.getParameter<double>("pv_max_abs_z");
  pv_max_rho_      = cfg.getParameter<double>("pv_max_rho");

  trk_max_normchi2_      = cfg.getParameter<double>("trk_max_normchi2");
  trk_max_abs_dz_        = cfg.getParameter<double>("trk_max_abs_dz");
  trk_max_abs_dxy_       = cfg.getParameter<double>("trk_max_abs_dxy");
  trk_max_lostInnerHits_ = static_cast<uint8_t>(cfg.getParameter<int>("trk_max_lostInnerHits"));

  puppi_drmax_          = cfg.getParameter<double>("puppi_drmax");
  puppi_eps_            = cfg.getParameter<double>("puppi_eps");
  puppi_max_neighbors_  = cfg.getParameter<int>("puppi_max_neighbors");
  puppi_beta_pow_       = cfg.getParameter<double>("puppi_beta_pow");

  scMuToken_  = consumes<std::vector<Run3ScoutingMuon>>(cfg.getParameter<edm::InputTag>("scoutingMuons"));
  scVtxToken_ = consumes<std::vector<Run3ScoutingVertex>>(cfg.getParameter<edm::InputTag>("scoutingVertices"));
  scPfcToken_ = consumes<std::vector<Run3ScoutingParticle>>(cfg.getParameter<edm::InputTag>("scoutingParticles"));

  scMetPtToken_  = consumes<double>(cfg.getParameter<edm::InputTag>("pfMetPtValue"));
  scMetPhiToken_ = consumes<double>(cfg.getParameter<edm::InputTag>("pfMetPhiValue"));
  rhoToken_      = consumes<double>(cfg.getParameter<edm::InputTag>("rhoValue"));

  puToken_      = consumes<std::vector<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("pileup"));
  genToken_     = consumes<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"));
  genMetToken_  = consumes<std::vector<reco::GenMET>>(cfg.getParameter<edm::InputTag>("genMetTrue"));
}

void ScoutingWAnalyzer::beginJob() {
  h_nTrueInt_ = fs_->make<TH1F>("nTrueInt", "true interactions; nTrueInt; events", 80, 0, 80);

  // MET resolution: (MET - genMET) [GeV] vs nTrueInt
  h_met_res_vs_nTrue_hlt_   = fs_->make<TH2F>("metRes_vs_nTrue_hlt",   "HLT MET - genMET; nTrueInt; MET-genMET [GeV]", 80, 0, 80, 200, -200, 200);
  h_met_res_vs_nTrue_chs_   = fs_->make<TH2F>("metRes_vs_nTrue_chs",   "CHS MET - genMET; nTrueInt; MET-genMET [GeV]", 80, 0, 80, 200, -200, 200);
  h_met_res_vs_nTrue_puppi_ = fs_->make<TH2F>("metRes_vs_nTrue_puppi", "PUPPI-lite MET - genMET; nTrueInt; MET-genMET [GeV]", 80, 0, 80, 200, -200, 200);

  // dphi(MET, genMET) vs nTrueInt
  h_dphi_vs_nTrue_hlt_   = fs_->make<TH2F>("dphi_vs_nTrue_hlt",   "dphi(HLT MET, genMET); nTrueInt; #Delta#phi", 80, 0, 80, 64, 0, 3.2);
  h_dphi_vs_nTrue_chs_   = fs_->make<TH2F>("dphi_vs_nTrue_chs",   "dphi(CHS MET, genMET); nTrueInt; #Delta#phi", 80, 0, 80, 64, 0, 3.2);
  h_dphi_vs_nTrue_puppi_ = fs_->make<TH2F>("dphi_vs_nTrue_puppi", "dphi(PUPPI-lite MET, genMET); nTrueInt; #Delta#phi", 80, 0, 80, 64, 0, 3.2);

  // response MET/genMET vs nTrueInt
  h_resp_vs_nTrue_hlt_   = fs_->make<TH2F>("resp_vs_nTrue_hlt",   "HLT MET / genMET; nTrueInt; MET/genMET", 80, 0, 80, 120, 0, 3.0);
  h_resp_vs_nTrue_chs_   = fs_->make<TH2F>("resp_vs_nTrue_chs",   "CHS MET / genMET; nTrueInt; MET/genMET", 80, 0, 80, 120, 0, 3.0);
  h_resp_vs_nTrue_puppi_ = fs_->make<TH2F>("resp_vs_nTrue_puppi", "PUPPI-lite MET / genMET; nTrueInt; MET/genMET", 80, 0, 80, 120, 0, 3.0);

  // MET pt vs nTrueInt (to see broadening)
  h_met_pt_vs_nTrue_hlt_   = fs_->make<TH2F>("metPt_vs_nTrue_hlt",   "HLT MET pT; nTrueInt; MET pT [GeV]", 80, 0, 80, 200, 0, 400);
  h_met_pt_vs_nTrue_chs_   = fs_->make<TH2F>("metPt_vs_nTrue_chs",   "CHS MET pT; nTrueInt; MET pT [GeV]", 80, 0, 80, 200, 0, 400);
  h_met_pt_vs_nTrue_puppi_ = fs_->make<TH2F>("metPt_vs_nTrue_puppi", "PUPPI-lite MET pT; nTrueInt; MET pT [GeV]", 80, 0, 80, 200, 0, 400);

  // Test
  h_chgPV_vs_chgPU_ = fs_->make<TH2F>(
    "chgPV_vs_chgPU",
    "charged counts; N(chg PV); N(chg PU)",
    200, 0, 2000, 200, 0, 2000
  );

  h_fracChgPU_vs_nTrue_ = fs_->make<TH2F>(
    "fracChgPU_vs_nTrue",
    "charged PU fraction; nTrueInt; N(chg PU)/(N(chg PV)+N(chg PU))",
    80, 0, 80, 110, 0.0, 1.1
  );

  h_chg_vertexIdx_ = fs_->make<TH1I>(
    "chg_vertexIdx",
    "Run3ScoutingParticle::vertex() for charged; vertex index; candidates",
    41, -1, 40
  );

  h_pvIndex_ = fs_->make<TH1I>(
    "pvIndex",
    "picked PV index; pvIndex; events",
    41, -1, 40
  );

  h_pfc_pdgId_ = fs_->make<TH1I>("pfc_pdgId", "Run3ScoutingParticle pdgId; pdgId; candidates", 401, -200, 200);
  h_pfc_absPdgId_ = fs_->make<TH1I>("pfc_absPdgId", "|pdgId|; |pdgId|; candidates", 201, 0, 200);
  h_pfc_pt_ = fs_->make<TH1F>("pfc_pt", "PF cand p_{T}; p_{T} [GeV]; candidates", 200, 0, 200);

  h_cnt_ = fs_->make<TH1I>("debug_counts",
    "counts; step; entries",
    6, 0, 6);
  h_cnt_->GetXaxis()->SetBinLabel(1, "events");
  h_cnt_->GetXaxis()->SetBinLabel(2, "pfc_total");
  h_cnt_->GetXaxis()->SetBinLabel(3, "pdg!=0");
  h_cnt_->GetXaxis()->SetBinLabel(4, "chargedByMap");
  h_cnt_->GetXaxis()->SetBinLabel(5, "passTrkQuality");
  h_cnt_->GetXaxis()->SetBinLabel(6, "filled_chg_vertexIdx");

}

int ScoutingWAnalyzer::pickBestPV_(const std::vector<Run3ScoutingVertex>& vtxs) const {
  int best = -1;
  int best_ntrk = -1;
  float best_score = -1e9f;

  for (int i = 0; i < (int)vtxs.size(); ++i) {
    const auto& v = vtxs[i];
    if (!v.isValidVtx()) continue;
    if (v.ndof() < pv_min_ndof_) continue;
    if (std::abs(v.z()) > pv_max_abs_z_) continue;
    float rho = std::sqrt(v.x()*v.x() + v.y()*v.y());
    if (rho > pv_max_rho_) continue;

    // Score: prefer more tracks; break ties with lower chi2/ndof
    float chi2ndof = (v.ndof() > 0) ? (v.chi2() / v.ndof()) : 999.f;
    float score = (float)v.tracksSize() - 0.1f * chi2ndof;

    if (v.tracksSize() > best_ntrk || score > best_score) {
      best = i;
      best_ntrk = v.tracksSize();
      best_score = score;
    }
  }
  return best;
}

std::pair<float,float> ScoutingWAnalyzer::buildMetCHS_(
    const std::vector<Run3ScoutingParticle>& pfc,
    int pvIndex) const {
  PtVec2 sum;
  for (const auto& p : pfc) {
    if (!isChargedFromPdgId(p.pdgId())) continue;
    if (p.vertex() != pvIndex) continue;
    if (!passTrkQuality_(p)) continue;
    sum.addPtPhi(p.pt(), p.phi(), 1.f);
  }
  // MET vector is negative of visible sum
  PtVec2 met;
  met.px = -sum.px; met.py = -sum.py;
  return {met.pt(), met.phi()};
}

std::pair<float,float> ScoutingWAnalyzer::buildMetPuppiLite_(
    const std::vector<Run3ScoutingParticle>& pfc,
    int pvIndex) const {

  // Pre-collect charged PV and charged PU lists (with quality cuts)
  struct Chg {
    float pt, eta, phi;
    int vtx;
  };
  std::vector<Chg> chgPV;
  std::vector<Chg> chgPU;
  chgPV.reserve(pfc.size());
  chgPU.reserve(pfc.size());

  for (const auto& p : pfc) {
    if (!isChargedFromPdgId(p.pdgId())) continue;
    if (!passTrkQuality_(p)) continue;
    Chg c{p.pt(), p.eta(), p.phi(), p.vertex()};
    if (p.vertex() == pvIndex) chgPV.push_back(c);
    else chgPU.push_back(c);
  }

  // Charged PV contribution: weight 1
  PtVec2 vis;
  for (const auto& c : chgPV) vis.addPtPhi(c.pt, c.phi, 1.f);

  // Neutrals: weight from local charged PV vs PU activity
  const float dr2max = puppi_drmax_ * puppi_drmax_;
  const float eps = puppi_eps_;

  for (const auto& p : pfc) {
    if (!isNeutralFromPdgId(p.pdgId())) continue;

    // compute alpha (PV) and beta (PU) from nearby charged tracks
    float alpha = 0.f, beta = 0.f;

    // PV neighbors
    int ncount = 0;
    for (const auto& c : chgPV) {
      if (ncount >= puppi_max_neighbors_) break;
      float dr2 = dR2(p.eta(), p.phi(), c.eta, c.phi);
      if (dr2 > dr2max) continue;
      float dr = std::sqrt(dr2);
      float denom = std::pow(dr + eps, puppi_beta_pow_);
      alpha += c.pt / denom;
      ++ncount;
    }

    // PU neighbors
    ncount = 0;
    for (const auto& c : chgPU) {
      if (ncount >= puppi_max_neighbors_) break;
      float dr2 = dR2(p.eta(), p.phi(), c.eta, c.phi);
      if (dr2 > dr2max) continue;
      float dr = std::sqrt(dr2);
      float denom = std::pow(dr + eps, puppi_beta_pow_);
      beta += c.pt / denom;
      ++ncount;
    }

    float w = 0.f;
    if (alpha + beta > 0.f) w = alpha / (alpha + beta);
    w = clamp01(w);

    vis.addPtPhi(p.pt(), p.phi(), w);
  }

  PtVec2 met;
  met.px = -vis.px; met.py = -vis.py;
  return {met.pt(), met.phi()};
}

void ScoutingWAnalyzer::analyze(const edm::Event& ev, const edm::EventSetup&) {
  // nTrueInt
  float nTrueInt = -1.f;
  edm::Handle<std::vector<PileupSummaryInfo>> hPU;
  ev.getByToken(puToken_, hPU);
  if (hPU.isValid()) {
    for (const auto& pu : *hPU) {
      if (pu.getBunchCrossing() == 0) { nTrueInt = pu.getTrueNumInteractions(); break; }
    }
  }
  if (nTrueInt < 0.f) return;
  h_nTrueInt_->Fill(nTrueInt);

  // GEN
  edm::Handle<std::vector<reco::GenParticle>> hGen;
  ev.getByToken(genToken_, hGen);
  if (!hGen.isValid()) return;

  GenWmunu g = selectGenWmunu(*hGen);
  if (!g.ok) return;
  if (!(g.mu_pt > fid_mu_pt_min_ && std::abs(g.mu_eta) < fid_mu_abs_eta_max_)) return;

  // genMET true
  float gen_met_pt = -1.f, gen_met_phi = 0.f;
  edm::Handle<std::vector<reco::GenMET>> hGenMet;
  ev.getByToken(genMetToken_, hGenMet);
  if (!hGenMet.isValid() || hGenMet->empty()) return;
  gen_met_pt = hGenMet->at(0).pt();
  gen_met_phi = hGenMet->at(0).phi();
  if (gen_met_pt <= 0.f) return;

  // Scouting PV and PF candidates
  edm::Handle<std::vector<Run3ScoutingVertex>> hVtx;
  ev.getByToken(scVtxToken_, hVtx);
  if (!hVtx.isValid() || hVtx->empty()) return;
  int pvIndex = pickBestPV_(*hVtx);
  if (pvIndex < 0) return;

  edm::Handle<std::vector<Run3ScoutingParticle>> hPfc;
  ev.getByToken(scPfcToken_, hPfc);
  if (!hPfc.isValid()) return;

  // TEST ====================================================================================
  h_cnt_->Fill(0); // events

  int n_total = 0, n_pdg_nonzero = 0, n_chg_map = 0, n_pass_q = 0, n_filled = 0;

  for (const auto& p : *hPfc) {
    ++n_total;

    int pdg = p.pdgId();
    h_pfc_pdgId_->Fill(pdg);
    h_pfc_absPdgId_->Fill(std::abs(pdg));
    h_pfc_pt_->Fill(p.pt());

    if (pdg != 0) ++n_pdg_nonzero;

    bool chg = isChargedFromPdgId(pdg);
    if (chg) ++n_chg_map;

    if (!chg) continue;

    if (!passTrkQuality_(p)) continue;
    ++n_pass_q;

    h_chg_vertexIdx_->Fill(p.vertex());
    ++n_filled;
  }

  h_cnt_->Fill(1, n_total);
  h_cnt_->Fill(2, n_pdg_nonzero);
  h_cnt_->Fill(3, n_chg_map);
  h_cnt_->Fill(4, n_pass_q);
  h_cnt_->Fill(5, n_filled);

  h_pvIndex_->Fill(pvIndex);

  int nChgPV = 0;
  int nChgPU = 0;

  for (const auto& p : *hPfc) {
    if (!isChargedFromPdgId(p.pdgId())) continue;
    if (!passTrkQuality_(p)) continue;

    h_chg_vertexIdx_->Fill(p.vertex());

    if (p.vertex() == pvIndex) ++nChgPV;
    else ++nChgPU;
  }

  h_chgPV_vs_chgPU_->Fill(nChgPV, nChgPU);

  float denom = float(nChgPV + nChgPU);
  float frac = (denom > 0.f) ? (float(nChgPU) / denom) : 0.f;
  h_fracChgPU_vs_nTrue_->Fill(nTrueInt, frac);
  // TEST ====================================================================================

  // HLT scouting MET
  float sc_met_pt = -1.f, sc_met_phi = 0.f;
  edm::Handle<double> hMetPt, hMetPhi;
  ev.getByToken(scMetPtToken_, hMetPt);
  ev.getByToken(scMetPhiToken_, hMetPhi);
  if (hMetPt.isValid())  sc_met_pt  = (float)(*hMetPt);
  if (hMetPhi.isValid()) sc_met_phi = (float)(*hMetPhi);
  if (sc_met_pt < 0.f) return;

  // Build mitigated METs
  auto [chs_pt, chs_phi] = buildMetCHS_(*hPfc, pvIndex);
  auto [pup_pt, pup_phi] = buildMetPuppiLite_(*hPfc, pvIndex);

  // Fill pt vs PU
  h_met_pt_vs_nTrue_hlt_->Fill(nTrueInt, sc_met_pt);
  h_met_pt_vs_nTrue_chs_->Fill(nTrueInt, chs_pt);
  h_met_pt_vs_nTrue_puppi_->Fill(nTrueInt, pup_pt);

  // Resolution vs gen
  auto fill_one = [&](float met_pt, float met_phi,
                      TH2F* hRes, TH2F* hDphi, TH2F* hResp) {
    float res = met_pt - gen_met_pt;
    float dphi = std::abs(reco::deltaPhi(met_phi, gen_met_phi));
    float resp = met_pt / gen_met_pt;
    hRes->Fill(nTrueInt, res);
    hDphi->Fill(nTrueInt, dphi);
    hResp->Fill(nTrueInt, resp);
  };

  fill_one(sc_met_pt, sc_met_phi, h_met_res_vs_nTrue_hlt_,   h_dphi_vs_nTrue_hlt_,   h_resp_vs_nTrue_hlt_);
  fill_one(chs_pt,    chs_phi,    h_met_res_vs_nTrue_chs_,   h_dphi_vs_nTrue_chs_,   h_resp_vs_nTrue_chs_);
  fill_one(pup_pt,    pup_phi,    h_met_res_vs_nTrue_puppi_, h_dphi_vs_nTrue_puppi_, h_resp_vs_nTrue_puppi_);
}

DEFINE_FWK_MODULE(ScoutingWAnalyzer);