// ScoutingHiggsRho.cc
// Clean refactor (MC flag + strict GEN fiducial + scouting uses ONLY Run3ScoutingParticle)
// - evt: keep only nTrueInt
// - flags: removed
// - gen: keep all + fid; fid is STRICT (event enters fid only if gamma + both pions in fid and hcand exists)
// - gen jets: removed completely
// - scouting: no tracks; rho built from PFCands (OS pion pairs); inclusive + truthMatched
// - matching: GEN<->SCOUT using photons + PFCands
// - scouting jets: kept (inclusive diagnostics + optional "jet-as-rho" toy)

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Jet.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <algorithm>

namespace {
  constexpr float M_PI_CH = 0.13957039f;
  constexpr float M_RHO   = 0.77526f;
  constexpr float M_H     = 125.0f;

  inline float toF(double x) { return static_cast<float>(x); }
  inline float toF(float  x) { return x; }
  inline float toF(int    x) { return static_cast<float>(x); }

  inline void setP4_massless(TLorentzVector& v, float pt, float eta, float phi) { v.SetPtEtaPhiM(pt, eta, phi, 0.0); }
  inline void setP4_pion(TLorentzVector& v, float pt, float eta, float phi)     { v.SetPtEtaPhiM(pt, eta, phi, M_PI_CH); }

  inline bool passAbsEta(float eta, float maxAbsEta) { return (maxAbsEta <= 0.f) ? true : (std::abs(eta) < maxAbsEta); }

  inline float deltaR2(float eta1, float phi1, float eta2, float phi2) {
    const float dEta = eta1 - eta2;
    float dPhi = std::fabs(phi1 - phi2);
    if (dPhi > float(M_PI)) dPhi = float(2.0*M_PI) - dPhi;
    return dEta*dEta + dPhi*dPhi;
  }

  inline float absDeltaPhi(float phi1, float phi2) {
    float dphi = std::fabs(phi1 - phi2);
    if (dphi > float(M_PI)) dphi = float(2.0*M_PI) - dphi;
    return dphi;
  }

  struct SimpleCand {
    float pt{0.f}, eta{0.f}, phi{0.f}, mass{0.f};
    int charge{0}; // 0 for neutral
  };

  struct RhoCand {
    SimpleCand piPlus;
    SimpleCand piMinus;
    TLorentzVector p4;
    float mass{-1.f};
    float pt()  const { return toF(p4.Pt()); }
    float eta() const { return toF(p4.Eta()); }
    float phi() const { return toF(p4.Phi()); }
  };

  struct HCand {
    SimpleCand gamma;
    RhoCand rho;
    TLorentzVector p4;
    float mass{-1.f};
    float pt()  const { return toF(p4.Pt()); }
    float eta() const { return toF(p4.Eta()); }
    float phi() const { return toF(p4.Phi()); }
  };

  // Bounded DFS on mother pointers (works for packed/pruned).
  inline bool hasAncestorPdgIdAbs(const reco::Candidate* start, int absPdgId, int guardMax = 300) {
    if (!start) return false;
    std::vector<const reco::Candidate*> stack;
    stack.reserve(32);
    const size_t nm0 = start->numberOfMothers();
    for (size_t i = 0; i < nm0; ++i) {
      const reco::Candidate* m = start->mother(i);
      if (m) stack.push_back(m);
    }
    int guard = 0;
    while (!stack.empty() && guard++ < guardMax) {
      const reco::Candidate* cur = stack.back();
      stack.pop_back();
      if (!cur) continue;
      if (std::abs(cur->pdgId()) == absPdgId) return true;
      const size_t nm = cur->numberOfMothers();
      for (size_t i = 0; i < nm; ++i) {
        const reco::Candidate* mm = cur->mother(i);
        if (mm) stack.push_back(mm);
      }
    }
    return false;
  }

  // ----------------------------
  // PF kinematics policy (Run3ScoutingParticle)
  // ----------------------------
  inline bool pf_canUseTrkKinematics(const Run3ScoutingParticle& p, bool wantUseTrkVars) {
    if (!wantUseTrkVars) return false;
    // In NanoAOD-land you'd always have separate arrays; here scouting has a flag.
    // "relative_trk_vars()==true" means track vars are *relative* and not usable as absolute.
    return !p.relative_trk_vars();
  }
  inline float pf_pt (const Run3ScoutingParticle& p, bool wantUseTrkVars) { return pf_canUseTrkKinematics(p, wantUseTrkVars) ? toF(p.trk_pt())  : toF(p.pt()); }
  inline float pf_eta(const Run3ScoutingParticle& p, bool wantUseTrkVars) { return pf_canUseTrkKinematics(p, wantUseTrkVars) ? toF(p.trk_eta()) : toF(p.eta()); }
  inline float pf_phi(const Run3ScoutingParticle& p, bool wantUseTrkVars) { return pf_canUseTrkKinematics(p, wantUseTrkVars) ? toF(p.trk_phi()) : toF(p.phi()); }

  inline int chargeFromPdgId(int pdgId) { return (pdgId > 0) ? +1 : (pdgId < 0 ? -1 : 0); }

  // Choose best OS pion pair closest to rho mass (from a vector of SimpleCand).
  inline bool bestPionPairToRho(const std::vector<SimpleCand>& pions, SimpleCand& piPlus, SimpleCand& piMinus, RhoCand& rhoOut) {
    float bestDiff = std::numeric_limits<float>::infinity();
    bool found = false;
    TLorentzVector best; float bestMass = -1.f;
    SimpleCand bestP{}, bestM{};
    for (size_t i = 0; i < pions.size(); ++i) {
      for (size_t j = i + 1; j < pions.size(); ++j) {
        const auto& a = pions[i];
        const auto& b = pions[j];
        if (a.charge * b.charge >= 0) continue;
        TLorentzVector va, vb;
        setP4_pion(va, a.pt, a.eta, a.phi);
        setP4_pion(vb, b.pt, b.eta, b.phi);
        TLorentzVector v = va + vb;
        const float m = toF(v.M());
        const float diff = std::abs(m - M_RHO);
        if (diff < bestDiff) {
          bestDiff = diff; best = v; bestMass = m; found = true;
          if (a.charge > 0) { bestP = a; bestM = b; } else { bestP = b; bestM = a; }
        }
      }
    }
    if (!found) return false;
    piPlus = bestP; piMinus = bestM;
    rhoOut.piPlus = bestP; rhoOut.piMinus = bestM; rhoOut.p4 = best; rhoOut.mass = bestMass;
    return true;
  }

  // ----------------------------
  // Histogram packs (minimal + consistent)
  // ----------------------------
  struct GenPack {
    // Higgs
    TH1F* h_H_pt{nullptr};
    TH1F* h_H_eta{nullptr};

    // final-state objects
    TH1F* h_g_pt{nullptr};   TH1F* h_g_eta{nullptr};   TH1F* h_g_phi{nullptr};
    TH1F* h_pi_pt{nullptr};  TH1F* h_pi_eta{nullptr};  TH1F* h_pi_phi{nullptr};
    TH1F* h_rho_m{nullptr};  TH1F* h_rho_pt{nullptr};  TH1F* h_rho_eta{nullptr}; TH1F* h_rho_phi{nullptr};
    TH1F* h_h_m{nullptr};    TH1F* h_h_pt{nullptr};    TH1F* h_h_eta{nullptr};   TH1F* h_h_phi{nullptr};

    // correlations
    TH1F* h_dr_pi_pi{nullptr};
    TH1F* h_dr_g_rho{nullptr};
    TH1F* h_dphi_g_rho{nullptr};
    TH1F* h_ptbal_rho_g{nullptr};

    void book(TFileDirectory& dir, const std::string& label) {
      (void)label;
      h_H_pt  = dir.make<TH1F>("H_pt",  "GEN Higgs p_{T};p_{T} [GeV];events", 200, 0., 400.);
      h_H_eta = dir.make<TH1F>("H_eta", "GEN Higgs #eta;#eta;events",         120, -6., 6.);

      h_g_pt  = dir.make<TH1F>("g_pt",  "GEN #gamma p_{T};p_{T} [GeV];events", 200, 0., 400.);
      h_g_eta = dir.make<TH1F>("g_eta", "GEN #gamma #eta;#eta;events",         120, -6., 6.);
      h_g_phi = dir.make<TH1F>("g_phi", "GEN #gamma #phi;#phi;events",         128, -3.2, 3.2);

      h_pi_pt  = dir.make<TH1F>("pi_pt",  "GEN #pi^{#pm} p_{T};p_{T} [GeV];cands", 200, 0., 100.);
      h_pi_eta = dir.make<TH1F>("pi_eta", "GEN #pi^{#pm} #eta;#eta;cands",         120, -6., 6.);
      h_pi_phi = dir.make<TH1F>("pi_phi", "GEN #pi^{#pm} #phi;#phi;cands",         128, -3.2, 3.2);

      h_rho_m   = dir.make<TH1F>("rho_m",   "GEN m(#pi#pi) (best OS);m [GeV];events", 200, 0., 2.0);
      h_rho_pt  = dir.make<TH1F>("rho_pt",  "GEN #rho p_{T};p_{T} [GeV];events",      200, 0., 250.);
      h_rho_eta = dir.make<TH1F>("rho_eta", "GEN #rho #eta;#eta;events",              120, -6., 6.);
      h_rho_phi = dir.make<TH1F>("rho_phi", "GEN #rho #phi;#phi;events",              128, -3.2, 3.2);

      h_h_m   = dir.make<TH1F>("hcand_m",   "GEN m(#pi#pi#gamma);m [GeV];events", 240, 0., 240.);
      h_h_pt  = dir.make<TH1F>("hcand_pt",  "GEN p_{T}(#pi#pi#gamma);p_{T} [GeV];events", 200, 0., 400.);
      h_h_eta = dir.make<TH1F>("hcand_eta", "GEN #eta(#pi#pi#gamma);#eta;events",        120, -6., 6.);
      h_h_phi = dir.make<TH1F>("hcand_phi", "GEN #phi(#pi#pi#gamma);#phi;events",        128, -3.2, 3.2);

      h_dr_pi_pi   = dir.make<TH1F>("dr_pi_pi",   "GEN #DeltaR(#pi,#pi);#DeltaR;events",           120, 0., 6.0);
      h_dr_g_rho   = dir.make<TH1F>("dr_g_rho",   "GEN #DeltaR(#gamma,#rho);#DeltaR;events",       120, 0., 6.0);
      h_dphi_g_rho = dir.make<TH1F>("dphi_g_rho", "GEN |#Delta#phi(#gamma,#rho)|;|#Delta#phi|;events", 128, 0., 3.2);
      h_ptbal_rho_g= dir.make<TH1F>("ptbal_rho_g","GEN p_{T}(#rho)/p_{T}(#gamma);ratio;events",    120, 0., 3.0);
    }

    void fillH(const reco::GenParticle& H) {
      if (h_H_pt)  h_H_pt->Fill(toF(H.pt()));
      if (h_H_eta) h_H_eta->Fill(toF(H.eta()));
    }

    void fillPions(const std::vector<SimpleCand>& pions) {
      for (const auto& pi : pions) {
        if (h_pi_pt)  h_pi_pt->Fill(pi.pt);
        if (h_pi_eta) h_pi_eta->Fill(pi.eta);
        if (h_pi_phi) h_pi_phi->Fill(pi.phi);
      }
    }

    void fillGamma(const SimpleCand& g) {
      if (h_g_pt)  h_g_pt->Fill(g.pt);
      if (h_g_eta) h_g_eta->Fill(g.eta);
      if (h_g_phi) h_g_phi->Fill(g.phi);
    }

    void fillRho(const RhoCand& r) {
      if (h_rho_m)   h_rho_m->Fill(r.mass);
      if (h_rho_pt)  h_rho_pt->Fill(r.pt());
      if (h_rho_eta) h_rho_eta->Fill(r.eta());
      if (h_rho_phi) h_rho_phi->Fill(r.phi());
      if (h_dr_pi_pi) {
        const float dr = std::sqrt(deltaR2(r.piPlus.eta, r.piPlus.phi, r.piMinus.eta, r.piMinus.phi));
        h_dr_pi_pi->Fill(dr);
      }
    }

    void fillHcand(const HCand& h) {
      if (h_h_m)   h_h_m->Fill(h.mass);
      if (h_h_pt)  h_h_pt->Fill(h.pt());
      if (h_h_eta) h_h_eta->Fill(h.eta());
      if (h_h_phi) h_h_phi->Fill(h.phi());
      if (h_dr_g_rho)   h_dr_g_rho->Fill(std::sqrt(deltaR2(h.gamma.eta, h.gamma.phi, h.rho.eta(), h.rho.phi())));
      if (h_dphi_g_rho) h_dphi_g_rho->Fill(absDeltaPhi(h.gamma.phi, h.rho.phi()));
      if (h_ptbal_rho_g && h.gamma.pt > 0.f) h_ptbal_rho_g->Fill(h.rho.pt() / h.gamma.pt);
    }
  };

  struct ScoutPack {
    // object-level
    TH1F* h_g_pt{nullptr};   TH1F* h_g_eta{nullptr};   TH1F* h_g_phi{nullptr};
    TH1F* h_pi_pt{nullptr};  TH1F* h_pi_eta{nullptr};  TH1F* h_pi_phi{nullptr};
    TH1F* h_rho_m{nullptr};  TH1F* h_rho_pt{nullptr};  TH1F* h_rho_eta{nullptr}; TH1F* h_rho_phi{nullptr};
    TH1F* h_h_m{nullptr};    TH1F* h_h_pt{nullptr};    TH1F* h_h_eta{nullptr};   TH1F* h_h_phi{nullptr};

    // event-candidate correlations
    TH1F* h_dr_pi_pi{nullptr};
    TH1F* h_dr_g_rho{nullptr};
    TH1F* h_dphi_g_rho{nullptr};
    TH1F* h_ptbal_rho_g{nullptr};

    void book(TFileDirectory& dir) {
      h_g_pt  = dir.make<TH1F>("g_pt",  "SCOUT #gamma p_{T};p_{T} [GeV];events", 200, 0., 400.);
      h_g_eta = dir.make<TH1F>("g_eta", "SCOUT #gamma #eta;#eta;events",         120, -6., 6.);
      h_g_phi = dir.make<TH1F>("g_phi", "SCOUT #gamma #phi;#phi;events",         128, -3.2, 3.2);

      h_pi_pt  = dir.make<TH1F>("pi_pt",  "SCOUT #pi^{#pm} p_{T};p_{T} [GeV];cands", 200, 0., 100.);
      h_pi_eta = dir.make<TH1F>("pi_eta", "SCOUT #pi^{#pm} #eta;#eta;cands",         120, -6., 6.);
      h_pi_phi = dir.make<TH1F>("pi_phi", "SCOUT #pi^{#pm} #phi;#phi;cands",         128, -3.2, 3.2);

      h_rho_m   = dir.make<TH1F>("rho_m",   "SCOUT m(#pi#pi) (best OS);m [GeV];events", 200, 0., 2.0);
      h_rho_pt  = dir.make<TH1F>("rho_pt",  "SCOUT #rho p_{T};p_{T} [GeV];events",      200, 0., 250.);
      h_rho_eta = dir.make<TH1F>("rho_eta", "SCOUT #rho #eta;#eta;events",              120, -6., 6.);
      h_rho_phi = dir.make<TH1F>("rho_phi", "SCOUT #rho #phi;#phi;events",              128, -3.2, 3.2);

      h_h_m   = dir.make<TH1F>("hcand_m",   "SCOUT m(#pi#pi#gamma);m [GeV];events", 240, 0., 240.);
      h_h_pt  = dir.make<TH1F>("hcand_pt",  "SCOUT p_{T}(#pi#pi#gamma);p_{T} [GeV];events", 200, 0., 400.);
      h_h_eta = dir.make<TH1F>("hcand_eta", "SCOUT #eta(#pi#pi#gamma);#eta;events",        120, -6., 6.);
      h_h_phi = dir.make<TH1F>("hcand_phi", "SCOUT #phi(#pi#pi#gamma);#phi;events",        128, -3.2, 3.2);

      h_dr_pi_pi   = dir.make<TH1F>("dr_pi_pi",   "SCOUT #DeltaR(#pi,#pi);#DeltaR;events",           120, 0., 6.0);
      h_dr_g_rho   = dir.make<TH1F>("dr_g_rho",   "SCOUT #DeltaR(#gamma,#rho);#DeltaR;events",       120, 0., 6.0);
      h_dphi_g_rho = dir.make<TH1F>("dphi_g_rho", "SCOUT |#Delta#phi(#gamma,#rho)|;|#Delta#phi|;events", 128, 0., 3.2);
      h_ptbal_rho_g= dir.make<TH1F>("ptbal_rho_g","SCOUT p_{T}(#rho)/p_{T}(#gamma);ratio;events",    120, 0., 3.0);
    }

    void fillGamma(const SimpleCand& g) {
      if (h_g_pt)  h_g_pt->Fill(g.pt);
      if (h_g_eta) h_g_eta->Fill(g.eta);
      if (h_g_phi) h_g_phi->Fill(g.phi);
    }

    void fillPions(const std::vector<SimpleCand>& pions) {
      for (const auto& pi : pions) {
        if (h_pi_pt)  h_pi_pt->Fill(pi.pt);
        if (h_pi_eta) h_pi_eta->Fill(pi.eta);
        if (h_pi_phi) h_pi_phi->Fill(pi.phi);
      }
    }

    void fillRho(const RhoCand& r) {
      if (h_rho_m)   h_rho_m->Fill(r.mass);
      if (h_rho_pt)  h_rho_pt->Fill(r.pt());
      if (h_rho_eta) h_rho_eta->Fill(r.eta());
      if (h_rho_phi) h_rho_phi->Fill(r.phi());
      if (h_dr_pi_pi) {
        const float dr = std::sqrt(deltaR2(r.piPlus.eta, r.piPlus.phi, r.piMinus.eta, r.piMinus.phi));
        h_dr_pi_pi->Fill(dr);
      }
    }

    void fillHcand(const HCand& h) {
      if (h_h_m)   h_h_m->Fill(h.mass);
      if (h_h_pt)  h_h_pt->Fill(h.pt());
      if (h_h_eta) h_h_eta->Fill(h.eta());
      if (h_h_phi) h_h_phi->Fill(h.phi());
      if (h_dr_g_rho)   h_dr_g_rho->Fill(std::sqrt(deltaR2(h.gamma.eta, h.gamma.phi, h.rho.eta(), h.rho.phi())));
      if (h_dphi_g_rho) h_dphi_g_rho->Fill(absDeltaPhi(h.gamma.phi, h.rho.phi()));
      if (h_ptbal_rho_g && h.gamma.pt > 0.f) h_ptbal_rho_g->Fill(h.rho.pt() / h.gamma.pt);
    }
  };

  struct BestCandPack {
    // “best-per-event” after cuts (keep it simple; no massive zoo)
    TH1F* h_rho_m{nullptr}; TH1F* h_rho_pt{nullptr};
    TH1F* h_h_m{nullptr};   TH1F* h_h_pt{nullptr};
    TH1F* h_dr_pi_pi{nullptr};
    TH1F* h_dr_g_rho{nullptr};
    TH1F* h_dphi_g_rho{nullptr};
    TH1F* h_ptbal_rho_g{nullptr};

    void book(TFileDirectory& dir) {
      h_rho_m  = dir.make<TH1F>("rho_m_best",  "BEST: m(#pi#pi);m [GeV];events", 200, 0., 2.0);
      h_rho_pt = dir.make<TH1F>("rho_pt_best", "BEST: p_{T}(#rho);p_{T} [GeV];events", 200, 0., 250.);
      h_h_m    = dir.make<TH1F>("hcand_m_best","BEST: m(#pi#pi#gamma);m [GeV];events", 240, 0., 240.);
      h_h_pt   = dir.make<TH1F>("hcand_pt_best","BEST: p_{T}(Hcand);p_{T} [GeV];events", 200, 0., 400.);
      h_dr_pi_pi  = dir.make<TH1F>("dr_pi_pi_best", "BEST: #DeltaR(#pi,#pi);#DeltaR;events", 120, 0., 2.0);
      h_dr_g_rho  = dir.make<TH1F>("dr_g_rho_best", "BEST: #DeltaR(#gamma,#rho);#DeltaR;events", 120, 0., 6.0);
      h_dphi_g_rho= dir.make<TH1F>("dphi_g_rho_best","BEST: |#Delta#phi(#gamma,#rho)|;|#Delta#phi|;events", 128, 0., 3.2);
      h_ptbal_rho_g = dir.make<TH1F>("ptbal_rho_g_best","BEST: p_{T}(#rho)/p_{T}(#gamma);ratio;events", 120, 0., 3.0);
    }

    void fill(const HCand& h) {
      const float dr_pp = std::sqrt(deltaR2(h.rho.piPlus.eta, h.rho.piPlus.phi, h.rho.piMinus.eta, h.rho.piMinus.phi));
      const float dr_gr = std::sqrt(deltaR2(h.gamma.eta, h.gamma.phi, h.rho.eta(), h.rho.phi()));
      const float dphi  = absDeltaPhi(h.gamma.phi, h.rho.phi());
      const float ptbal = (h.gamma.pt > 0.f) ? (h.rho.pt() / h.gamma.pt) : 0.f;

      if (h_rho_m)  h_rho_m->Fill(h.rho.mass);
      if (h_rho_pt) h_rho_pt->Fill(h.rho.pt());
      if (h_h_m)    h_h_m->Fill(h.mass);
      if (h_h_pt)   h_h_pt->Fill(h.pt());
      if (h_dr_pi_pi)   h_dr_pi_pi->Fill(dr_pp);
      if (h_dr_g_rho)   h_dr_g_rho->Fill(dr_gr);
      if (h_dphi_g_rho) h_dphi_g_rho->Fill(dphi);
      if (h_ptbal_rho_g) h_ptbal_rho_g->Fill(ptbal);
    }
  };

  struct MatchPack {
    TH1F* h_dr2_g{nullptr};
    TH1F* h_dr2_piP{nullptr};
    TH1F* h_dr2_piM{nullptr};
    TH1F* h_resp_g_pt{nullptr};
    TH1F* h_resp_rho_m{nullptr};
    TH1F* h_resp_h_m{nullptr};

    TH1F* h_eff_g_den{nullptr}; TH1F* h_eff_g_num{nullptr};
    TH1F* h_eff_rho_den{nullptr}; TH1F* h_eff_rho_num{nullptr};
    TH1F* h_eff_g_den_fid{nullptr}; TH1F* h_eff_g_num_fid{nullptr};
    TH1F* h_eff_rho_den_fid{nullptr}; TH1F* h_eff_rho_num_fid{nullptr};

    void book(TFileDirectory& dir) {
      h_dr2_g   = dir.make<TH1F>("dr2_g",   "#DeltaR^{2}(GEN #gamma, SCOUT #gamma);#DeltaR^{2};events", 120, 0., 0.5);
      h_dr2_piP = dir.make<TH1F>("dr2_piP", "#DeltaR^{2}(GEN #pi^{+}, SCOUT pfCand);#DeltaR^{2};events", 120, 0., 0.5);
      h_dr2_piM = dir.make<TH1F>("dr2_piM", "#DeltaR^{2}(GEN #pi^{-}, SCOUT pfCand);#DeltaR^{2};events", 120, 0., 0.5);

      h_resp_g_pt  = dir.make<TH1F>("resp_g_pt",  "p_{T} response #gamma; p_{T}^{scout}/p_{T}^{gen};events", 120, 0., 3.0);
      h_resp_rho_m = dir.make<TH1F>("resp_rho_m", "mass response #rho; m^{scout}-m^{gen} [GeV];events", 160, -0.4, 0.4);
      h_resp_h_m   = dir.make<TH1F>("resp_h_m",   "mass response Hcand; m^{scout}-m^{gen} [GeV];events", 200, -10., 10.);

      h_eff_g_den   = dir.make<TH1F>("eff_g_den",   "DEN: GEN #gamma p_{T};p_{T}^{gen} [GeV];events", 200, 0., 400.);
      h_eff_g_num   = dir.make<TH1F>("eff_g_num",   "NUM: matched SCOUT #gamma;p_{T}^{gen} [GeV];events", 200, 0., 400.);
      h_eff_rho_den = dir.make<TH1F>("eff_rho_den", "DEN: GEN #rho p_{T};p_{T}^{gen} [GeV];events", 200, 0., 250.);
      h_eff_rho_num = dir.make<TH1F>("eff_rho_num", "NUM: matched SCOUT #rho;p_{T}^{gen} [GeV];events", 200, 0., 250.);

      h_eff_g_den_fid   = dir.make<TH1F>("eff_g_den_fid",   "DEN(FID): GEN #gamma p_{T};p_{T}^{gen} [GeV];events", 200, 0., 400.);
      h_eff_g_num_fid   = dir.make<TH1F>("eff_g_num_fid",   "NUM(FID): matched SCOUT #gamma;p_{T}^{gen} [GeV];events", 200, 0., 400.);
      h_eff_rho_den_fid = dir.make<TH1F>("eff_rho_den_fid", "DEN(FID): GEN #rho p_{T};p_{T}^{gen} [GeV];events", 200, 0., 250.);
      h_eff_rho_num_fid = dir.make<TH1F>("eff_rho_num_fid", "NUM(FID): matched SCOUT #rho;p_{T}^{gen} [GeV];events", 200, 0., 250.);
    }
  };

  struct ScoutJetPack {
    TH1F* h_njet{nullptr};
    TH1F* h_ht{nullptr};
    TH1F* h_j1_pt{nullptr}; TH1F* h_j1_eta{nullptr}; TH1F* h_j1_phi{nullptr}; TH1F* h_j1_m{nullptr}; TH1F* h_j1_area{nullptr};
    TH1F* h_min_dr2_g_j{nullptr};
    TH1F* h_min_dr2_rho_j{nullptr};

    // toy: jet-as-rho
    TH1F* h_h_m_jetAsRho{nullptr};
    TH1F* h_h_pt_jetAsRho{nullptr};

    void book(TFileDirectory& dir) {
      h_njet = dir.make<TH1F>("njet", "N(jets) (selection);N;events", 30, -0.5, 29.5);
      h_ht   = dir.make<TH1F>("ht",   "H_{T} (sum jet p_{T});H_{T} [GeV];events", 400, 0., 1200.);

      h_j1_pt   = dir.make<TH1F>("j1_pt",   "leading jet p_{T};p_{T} [GeV];events", 300, 0., 600.);
      h_j1_eta  = dir.make<TH1F>("j1_eta",  "leading jet #eta;#eta;events",         160, -5., 5.);
      h_j1_phi  = dir.make<TH1F>("j1_phi",  "leading jet #phi;#phi;events",         128, -3.2, 3.2);
      h_j1_m    = dir.make<TH1F>("j1_m",    "leading jet mass;m [GeV];events",      200, 0., 200.);
      h_j1_area = dir.make<TH1F>("j1_area", "leading jet area;area;events",         120, 0., 2.0);

      h_min_dr2_g_j   = dir.make<TH1F>("min_dr2_g_j",   "min #DeltaR^{2}(#gamma,jet);#DeltaR^{2};events", 200, 0., 25.);
      h_min_dr2_rho_j = dir.make<TH1F>("min_dr2_rho_j", "min #DeltaR^{2}(#rho,jet);#DeltaR^{2};events",   200, 0., 25.);

      h_h_m_jetAsRho  = dir.make<TH1F>("hcand_m_jetAsRho",  "TOY: m(jet+#gamma);m [GeV];events", 300, 0., 600.);
      h_h_pt_jetAsRho = dir.make<TH1F>("hcand_pt_jetAsRho", "TOY: p_{T}(jet+#gamma);p_{T} [GeV];events", 300, 0., 600.);
    }
  };

  struct L1Pack {
    // object-level
    TH1F* h_neg{nullptr};
    TH1F* h_njet{nullptr};

    TH1F* h_eg_pt{nullptr};  TH1F* h_eg_eta{nullptr};  TH1F* h_eg_phi{nullptr};
    TH1F* h_j_pt{nullptr};   TH1F* h_j_eta{nullptr};   TH1F* h_j_phi{nullptr};

    // best-pair per event (bx=0)
    TH1F* h_m_egj{nullptr};
    TH1F* h_pt_egj{nullptr};
    TH1F* h_dr_eg_j{nullptr};
    TH1F* h_dphi_eg_j{nullptr};
    TH1F* h_ptbal_j_eg{nullptr};

    // optional GEN matching diagnostics (only if isMC_ && hasG)
    TH1F* h_dr2_genG_l1EG{nullptr};
    TH1F* h_resp_eg_pt{nullptr};

    void book(TFileDirectory& dir) {
      h_neg  = dir.make<TH1F>("nEG",  "N(L1 EGamma) BX0;N;events", 30, -0.5, 29.5);
      h_njet = dir.make<TH1F>("nJet", "N(L1 Jet) BX0;N;events",   60, -0.5, 59.5);

      h_eg_pt  = dir.make<TH1F>("eg_pt",  "L1 EG p_{T};p_{T} [GeV];cands", 300, 0., 600.);
      h_eg_eta = dir.make<TH1F>("eg_eta", "L1 EG #eta;#eta;cands",         160, -5., 5.);
      h_eg_phi = dir.make<TH1F>("eg_phi", "L1 EG #phi;#phi;cands",         128, -3.2, 3.2);

      h_j_pt   = dir.make<TH1F>("jet_pt",  "L1 Jet p_{T};p_{T} [GeV];cands", 400, 0., 800.);
      h_j_eta  = dir.make<TH1F>("jet_eta", "L1 Jet #eta;#eta;cands",         160, -5., 5.);
      h_j_phi  = dir.make<TH1F>("jet_phi", "L1 Jet #phi;#phi;cands",         128, -3.2, 3.2);

      h_m_egj     = dir.make<TH1F>("m_egjet",   "L1 m(EG+Jet);m [GeV];events", 400, 0., 800.);
      h_pt_egj    = dir.make<TH1F>("pt_egjet",  "L1 p_{T}(EG+Jet);p_{T} [GeV];events", 400, 0., 800.);
      h_dr_eg_j   = dir.make<TH1F>("dr_eg_jet", "L1 #DeltaR(EG,Jet);#DeltaR;events", 200, 0., 6.0);
      h_dphi_eg_j = dir.make<TH1F>("dphi_eg_jet","L1 |#Delta#phi(EG,Jet)|;|#Delta#phi|;events", 128, 0., 3.2);
      h_ptbal_j_eg= dir.make<TH1F>("ptbal_j_eg","L1 p_{T}(Jet)/p_{T}(EG);ratio;events", 160, 0., 4.0);

      h_dr2_genG_l1EG = dir.make<TH1F>("dr2_genG_l1EG", "#DeltaR^{2}(GEN #gamma, L1 EG);#DeltaR^{2};events", 120, 0., 0.5);
      h_resp_eg_pt    = dir.make<TH1F>("resp_eg_pt", "L1 EG p_{T} response; p_{T}^{L1}/p_{T}^{gen};events", 160, 0., 4.0);
    }
  };

} // namespace

class ScoutingHiggsRho : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ScoutingHiggsRho(const edm::ParameterSet&);
  ~ScoutingHiggsRho() override = default;

  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  // ---------------- Config
  bool isMC_{true};

  edm::InputTag puTag_, prunedGenTag_, packedGenTag_;
  edm::InputTag phoTag_, pfcTag_;
  edm::InputTag jetTag_;
  edm::InputTag l1EGTag_, l1JetTag_;

  bool storeHistograms_{true};

  // GEN fiducial
  float fidPhoMaxAbsEta_{2.5f};
  float fidPiMaxAbsEta_{2.5f};
  float genPhoMinPt_{0.f};
  float genPiMinPt_{0.f};

  // SCOUT photon selection
  double phoMinPt_{0.0};
  double phoMaxAbsEta_{0.0};

  // SCOUT PFC quality
  double pfcMinPt_{1.0};
  double pfcMaxAbsEta_{2.4};
  double pfcMaxAbsDz_{0.2};
  double pfcMaxAbsDxy_{0.1};
  bool   pfPionOnly_{true};
  bool   pfUseTrkVars_{false};

  // best-cand cuts (inclusive)
  bool   pfRequireSameVtx_{true};
  double pfMaxAbsDzDiff_{0.05};
  double pfMaxAbsDxyDiff_{0.05};
  double pfMaxDR2PiPi_{0.25};          // e.g. (0.5)^2
  double pfRhoMinPt_{15.0};
  double pfMinAbsDPhiGRho_{2.7};
  double pfMinPtBal_{0.5};
  double pfMaxPtBal_{2.0};
  double pfRhoMassMin_{0.0};
  double pfRhoMassMax_{0.0};

  // matching
  bool  doTruthMatching_{true};
  float matchPhoDR2_{0.01f};
  float matchPionDR2_{0.01f};

  // jets (scouting only)
  double jetMinPt_{30.0};
  double jetMaxAbsEta_{2.4};
  double htJetMinPt_{30.0};

  // ---------------- L1 (Stage-2 RECO): EGamma + Jet (BXVector)
  bool   doL1_{true};
  bool   l1UseBX0Only_{true};

  double l1EGMinPt_{0.0};
  double l1JetMinPt_{0.0};
  double l1EGMaxAbsEta_{0.0};   // 0 => no cut
  double l1JetMaxAbsEta_{0.0};  // 0 => no cut

  // optional GEN match for EG
  float  matchL1EG_DR2_{0.04f}; // (0.2)^2

  // ---------------- Tokens
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>>      puTok_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>>      prunedGenTok_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> packedGenTok_;

  edm::EDGetTokenT<std::vector<Run3ScoutingPhoton>>     scoutPhoTok_;
  edm::EDGetTokenT<std::vector<Run3ScoutingParticle>>   scoutPFCTok_;
  edm::EDGetTokenT<std::vector<Run3ScoutingPFJet>>      scoutJetTok_;

  edm::EDGetTokenT<l1t::EGammaBxCollection> l1EGTok_;
  edm::EDGetTokenT<l1t::JetBxCollection>    l1JetTok_;

  // ---------------- Output
  edm::Service<TFileService> fs_;

  // evt
  TH1F* h_evt_nTrueInt_{nullptr};

  // gen
  GenPack genAll_;
  GenPack genFid_;

  // scouting
  ScoutPack scoutIncl_;
  BestCandPack scoutBest_;
  ScoutPack scoutTruth_; // objects reconstructed from matched pieces

  // matching
  MatchPack match_;

  // jets
  ScoutJetPack jets_;

  // L1
  L1Pack l1_;

  // ---------------- helpers
  bool passScoutPFQuality_(const Run3ScoutingParticle& p) const {
    const float pt  = pf_pt(p, pfUseTrkVars_);
    const float eta = pf_eta(p, pfUseTrkVars_);
    if (pt < pfcMinPt_) return false;
    if (std::abs(eta) > pfcMaxAbsEta_) return false;
    if (std::abs(toF(p.dz()))  > pfcMaxAbsDz_)  return false;
    if (std::abs(toF(p.dxy())) > pfcMaxAbsDxy_) return false;
    return true;
  }

  bool selectBestScoutPhoton_(const std::vector<Run3ScoutingPhoton>& phos, SimpleCand& out) const {
    int best = -1; float bestPt = -1.f;
    for (size_t i = 0; i < phos.size(); ++i) {
      const auto& p = phos[i];
      const float pt = toF(p.pt());
      const float eta = toF(p.eta());
      if (pt < toF(phoMinPt_)) continue;
      if (phoMaxAbsEta_ > 0.0 && std::abs(eta) > toF(phoMaxAbsEta_)) continue;
      if (pt > bestPt) { bestPt = pt; best = int(i); }
    }
    if (best < 0) return false;
    const auto& p = phos[best];
    out.pt = toF(p.pt()); out.eta = toF(p.eta()); out.phi = toF(p.phi()); out.mass = 0.f; out.charge = 0;
    return true;
  }

  // If rho known, prefer photon that makes m(rho+g) close to 125 and back-to-back-ish.
  bool selectBestScoutPhotonGivenRho_(const std::vector<Run3ScoutingPhoton>& phos, const TLorentzVector& rhoP4,
                                      SimpleCand& out) const {
    int best = -1; float bestScore = 1e9f;
    for (size_t i = 0; i < phos.size(); ++i) {
      const auto& p = phos[i];
      const float pt  = toF(p.pt());
      const float eta = toF(p.eta());
      const float phi = toF(p.phi());
      if (pt < toF(phoMinPt_)) continue;
      if (phoMaxAbsEta_ > 0.0 && std::abs(eta) > toF(phoMaxAbsEta_)) continue;

      TLorentzVector gP4; setP4_massless(gP4, pt, eta, phi);
      const TLorentzVector hP4 = rhoP4 + gP4;
      const float mHcand = toF(hP4.M());

      const float dphi = absDeltaPhi(phi, toF(rhoP4.Phi()));
      const float pull_m = (mHcand - M_H) / 10.0f;
      const float pull_top = (dphi - float(M_PI)) / 0.6f;
      const float pull_pt = (pt - 60.0f) / 25.0f;

      const float score = pull_m*pull_m + 0.25f*pull_top*pull_top + 0.20f*pull_pt*pull_pt;
      if (score < bestScore) { bestScore = score; best = int(i); }
    }
    if (best < 0) return false;
    const auto& p = phos[best];
    out.pt = toF(p.pt()); out.eta = toF(p.eta()); out.phi = toF(p.phi()); out.mass = 0.f; out.charge = 0;
    return true;
  }

  // Build best rho from PFCands (OS pion-like), return also the two indices for further vertex/dxy/dz checks.
  bool bestPFCandPairToRho_(const std::vector<Run3ScoutingParticle>& pfc,
                           int& idxP, int& idxM, RhoCand& rhoOut,
                           std::vector<SimpleCand>* outAllPions = nullptr) const {
    idxP = -1; idxM = -1; rhoOut = RhoCand{};
    float bestDiff = std::numeric_limits<float>::infinity();
    bool found = false;

    std::vector<int> piIdx;
    piIdx.reserve(pfc.size());
    for (int i = 0; i < (int)pfc.size(); ++i) {
      const auto& p = pfc[i];
      if (pfPionOnly_ && std::abs(p.pdgId()) != 211) continue;
      if (!passScoutPFQuality_(p)) continue;
      const int q = chargeFromPdgId(p.pdgId());
      if (q == 0) continue;
      piIdx.push_back(i);

      if (outAllPions) {
        SimpleCand c;
        c.pt = pf_pt(p, pfUseTrkVars_);
        c.eta = pf_eta(p, pfUseTrkVars_);
        c.phi = pf_phi(p, pfUseTrkVars_);
        c.mass = M_PI_CH;
        c.charge = q;
        outAllPions->push_back(c);
      }
    }

    for (int a = 0; a < (int)piIdx.size(); ++a) {
      const auto& pa = pfc[piIdx[a]];
      const int qa = chargeFromPdgId(pa.pdgId());
      const float pta = pf_pt(pa, pfUseTrkVars_), etaa = pf_eta(pa, pfUseTrkVars_), phia = pf_phi(pa, pfUseTrkVars_);
      for (int b = a + 1; b < (int)piIdx.size(); ++b) {
        const auto& pb = pfc[piIdx[b]];
        const int qb = chargeFromPdgId(pb.pdgId());
        if (qa * qb >= 0) continue;

        const float ptb = pf_pt(pb, pfUseTrkVars_), etab = pf_eta(pb, pfUseTrkVars_), phib = pf_phi(pb, pfUseTrkVars_);

        TLorentzVector vA, vB;
        setP4_pion(vA, pta, etaa, phia);
        setP4_pion(vB, ptb, etab, phib);
        const TLorentzVector vR = vA + vB;
        const float mR = toF(vR.M());
        const float diff = std::abs(mR - M_RHO);
        if (diff < bestDiff) {
          bestDiff = diff;
          found = true;

          const int iA = piIdx[a];
          const int iB = piIdx[b];
          if (qa > 0) { idxP = iA; idxM = iB; }
          else        { idxP = iB; idxM = iA; }

          SimpleCand piP{ pf_pt(pfc[idxP], pfUseTrkVars_), pf_eta(pfc[idxP], pfUseTrkVars_), pf_phi(pfc[idxP], pfUseTrkVars_), M_PI_CH, +1 };
          SimpleCand piM{ pf_pt(pfc[idxM], pfUseTrkVars_), pf_eta(pfc[idxM], pfUseTrkVars_), pf_phi(pfc[idxM], pfUseTrkVars_), M_PI_CH, -1 };

          rhoOut.piPlus = piP;
          rhoOut.piMinus = piM;
          rhoOut.p4 = vR;
          rhoOut.mass = mR;
        }
      }
    }
    return found;
  }

  // ---------------- GEN helpers
  const reco::GenParticle* findLeadingHiggs_(const std::vector<reco::GenParticle>& pruned) const {
    const reco::GenParticle* higgs = nullptr;
    for (const auto& gp : pruned) {
      if (std::abs(gp.pdgId()) != 25) continue;
      if (!higgs || gp.pt() > higgs->pt()) higgs = &gp;
    }
    return higgs;
  }

  const reco::GenParticle* findGammaFromPrunedDesc_(const reco::GenParticle* higgs) const {
    if (!higgs) return nullptr;
    const reco::GenParticle* best = nullptr;
    float bestPt = -1.f;
    std::vector<const reco::Candidate*> stack;
    stack.reserve(64);
    stack.push_back(higgs);
    int guard = 0;
    while (!stack.empty() && guard++ < 2000) {
      const reco::Candidate* cur = stack.back();
      stack.pop_back();
      if (!cur) continue;
      if (cur->pdgId() == 22) {
        const float pt = toF(cur->pt());
        if (pt > bestPt) {
          bestPt = pt;
          best = dynamic_cast<const reco::GenParticle*>(cur); // may fail; ok
        }
      }
      const size_t nd = cur->numberOfDaughters();
      for (size_t i = 0; i < nd; ++i) {
        const reco::Candidate* d = cur->daughter(i);
        if (d) stack.push_back(d);
      }
    }
    return best;
  }

  void collectPackedFromHiggs_(const std::vector<pat::PackedGenParticle>& packed,
                               std::vector<SimpleCand>& pionsOut,
                               const pat::PackedGenParticle*& bestGammaOut) const {
    pionsOut.clear();
    pionsOut.reserve(128);
    bestGammaOut = nullptr;
    float bestGamPt = -1.f;

    for (const auto& pg : packed) {
      const int id = pg.pdgId();
      if (!hasAncestorPdgIdAbs(&pg, 25)) continue;

      if (std::abs(id) == 211) {
        const float pt = toF(pg.pt());
        if (pt < genPiMinPt_) continue;
        SimpleCand c;
        c.pt  = pt;
        c.eta = toF(pg.eta());
        c.phi = toF(pg.phi());
        c.mass = M_PI_CH;
        c.charge = (id == 211 ? +1 : -1);
        pionsOut.push_back(c);
        continue;
      }

      if (id == 22) {
        const float pt = toF(pg.pt());
        if (pt < genPhoMinPt_) continue;
        if (pt > bestGamPt) { bestGamPt = pt; bestGammaOut = &pg; }
      }
    }
  }

  // Build GEN gamma + best rho + hcand (fidOnly uses fid cuts for gamma and both pions)
  bool buildGenCandidates_(const edm::Event& ev, bool fidOnly,
                          bool& hasH, bool& hasG, bool& hasPiPi, bool& hasHcand,
                          const reco::GenParticle*& higgsOut,
                          SimpleCand& gammaOut, RhoCand& rhoOut, HCand& hcOut,
                          std::vector<SimpleCand>* outAllPions = nullptr) const {
    hasH = hasG = hasPiPi = hasHcand = false;
    higgsOut = nullptr;
    gammaOut = SimpleCand{};
    rhoOut = RhoCand{};
    hcOut = HCand{};
    hcOut.p4.SetPxPyPzE(0,0,0,0);
    if (outAllPions) outAllPions->clear();

    edm::Handle<std::vector<reco::GenParticle>> hPruned;
    ev.getByToken(prunedGenTok_, hPruned);
    if (!hPruned.isValid()) return false;

    const reco::GenParticle* higgs = findLeadingHiggs_(*hPruned);
    if (!higgs) return true;
    hasH = true;
    higgsOut = higgs;

    const reco::GenParticle* gamPruned = findGammaFromPrunedDesc_(higgs);

    edm::Handle<std::vector<pat::PackedGenParticle>> hPacked;
    ev.getByToken(packedGenTok_, hPacked);

    std::vector<SimpleCand> pionsPacked;
    const pat::PackedGenParticle* gamPacked = nullptr;
    if (hPacked.isValid()) collectPackedFromHiggs_(*hPacked, pionsPacked, gamPacked);

    auto acceptGammaFid = [&](float eta) { return passAbsEta(eta, fidPhoMaxAbsEta_); };
    auto acceptPionFid  = [&](float eta) { return passAbsEta(eta, fidPiMaxAbsEta_); };

    // gamma: prefer pruned, fallback packed
    if (gamPruned) {
      const float pt  = toF(gamPruned->pt());
      const float eta = toF(gamPruned->eta());
      const float phi = toF(gamPruned->phi());
      if (pt >= genPhoMinPt_ && (!fidOnly || acceptGammaFid(eta))) {
        hasG = true;
        gammaOut = SimpleCand{pt, eta, phi, 0.f, 0};
      }
    }
    if (!hasG && gamPacked) {
      const float pt  = toF(gamPacked->pt());
      const float eta = toF(gamPacked->eta());
      const float phi = toF(gamPacked->phi());
      if (pt >= genPhoMinPt_ && (!fidOnly || acceptGammaFid(eta))) {
        hasG = true;
        gammaOut = SimpleCand{pt, eta, phi, 0.f, 0};
      }
    }

    // pions: packed first; if none, fallback to pruned
    std::vector<SimpleCand> pionCands;
    pionCands.reserve(pionsPacked.size());
    for (const auto& pi : pionsPacked) {
      if (pi.pt < genPiMinPt_) continue;
      if (fidOnly && !acceptPionFid(pi.eta)) continue;
      pionCands.push_back(pi);
    }

    if (pionCands.empty()) {
      for (const auto& gp : *hPruned) {
        if (std::abs(gp.pdgId()) != 211) continue;
        if (!hasAncestorPdgIdAbs(&gp, 25)) continue;
        const float pt  = toF(gp.pt());
        const float eta = toF(gp.eta());
        const float phi = toF(gp.phi());
        if (pt < genPiMinPt_) continue;
        if (fidOnly && !acceptPionFid(eta)) continue;
        SimpleCand c{pt, eta, phi, M_PI_CH, (gp.pdgId()==211?+1:-1)};
        pionCands.push_back(c);
      }
    }

    if (outAllPions) *outAllPions = pionCands;

    SimpleCand piP, piM;
    RhoCand rho;
    if (!pionCands.empty() && bestPionPairToRho(pionCands, piP, piM, rho)) {
      hasPiPi = true;
      rhoOut = rho;
    }

    if (hasG && hasPiPi) {
      TLorentzVector vg, vh;
      setP4_massless(vg, gammaOut.pt, gammaOut.eta, gammaOut.phi);
      vh = vg + rhoOut.p4;
      hasHcand = true;
      hcOut.gamma = gammaOut;
      hcOut.rho = rhoOut;
      hcOut.p4 = vh;
      hcOut.mass = toF(vh.M());
    }

    // STRICT fidOnly: if any product missing -> treat as “not in fid category”
    if (fidOnly) {
      if (!(hasG && hasPiPi && hasHcand)) {
        // keep hasH to reflect that Higgs exists, but caller must not fill fid pack unless hasHcand==true
        hasG = false; hasPiPi = false; hasHcand = false;
      }
    }

    return true;
  }
};

ScoutingHiggsRho::ScoutingHiggsRho(const edm::ParameterSet& cfg) {
  usesResource("TFileService");

  storeHistograms_ = cfg.getUntrackedParameter<bool>("storeHistograms", true);
  isMC_            = cfg.getUntrackedParameter<bool>("isMC", true);

  puTag_        = cfg.getUntrackedParameter<edm::InputTag>("pileupSummary",      edm::InputTag("slimmedAddPileupInfo"));
  prunedGenTag_ = cfg.getUntrackedParameter<edm::InputTag>("prunedGenParticles", edm::InputTag("prunedGenParticles"));
  packedGenTag_ = cfg.getUntrackedParameter<edm::InputTag>("packedGenParticles", edm::InputTag("packedGenParticles"));

  fidPhoMaxAbsEta_ = toF(cfg.getUntrackedParameter<double>("fidPhoMaxAbsEta", 2.5));
  fidPiMaxAbsEta_  = toF(cfg.getUntrackedParameter<double>("fidPiMaxAbsEta",  2.5));
  genPhoMinPt_     = toF(cfg.getUntrackedParameter<double>("genPhoMinPt", 0.0));
  genPiMinPt_      = toF(cfg.getUntrackedParameter<double>("genPiMinPt",  0.0));

  phoTag_ = cfg.getParameter<edm::InputTag>("photonCollection");
  pfcTag_ = cfg.getParameter<edm::InputTag>("pfCandCollection");
  jetTag_ = cfg.getParameter<edm::InputTag>("pfJetCollection");

  phoMinPt_     = toF(cfg.getUntrackedParameter<double>("phoMinPt", 0.0));
  phoMaxAbsEta_ = toF(cfg.getUntrackedParameter<double>("phoMaxAbsEta", 0.0));

  pfcMinPt_     = toF(cfg.getUntrackedParameter<double>("pfcMinPt", 1.0));
  pfcMaxAbsEta_ = toF(cfg.getUntrackedParameter<double>("pfcMaxAbsEta", 2.4));
  pfcMaxAbsDz_  = toF(cfg.getUntrackedParameter<double>("pfcMaxAbsDz", 0.2));
  pfcMaxAbsDxy_ = toF(cfg.getUntrackedParameter<double>("pfcMaxAbsDxy", 0.1));
  pfPionOnly_   = cfg.getUntrackedParameter<bool>("pfPionOnly", true);
  pfUseTrkVars_ = cfg.getUntrackedParameter<bool>("pfUseTrkVars", false);

  pfRequireSameVtx_   = cfg.getUntrackedParameter<bool>("pfRequireSameVtx", true);
  pfMaxAbsDzDiff_     = toF(cfg.getUntrackedParameter<double>("pfMaxAbsDzDiff",  0.05));
  pfMaxAbsDxyDiff_    = toF(cfg.getUntrackedParameter<double>("pfMaxAbsDxyDiff", 0.05));
  pfMaxDR2PiPi_       = toF(cfg.getUntrackedParameter<double>("pfMaxDR2PiPi", 0.25));
  pfRhoMinPt_         = toF(cfg.getUntrackedParameter<double>("pfRhoMinPt", 15.0));
  pfMinAbsDPhiGRho_   = toF(cfg.getUntrackedParameter<double>("pfMinAbsDPhiGRho", 2.7));
  pfMinPtBal_         = toF(cfg.getUntrackedParameter<double>("pfMinPtBal", 0.5));
  pfMaxPtBal_         = toF(cfg.getUntrackedParameter<double>("pfMaxPtBal", 2.0));
  pfRhoMassMin_       = toF(cfg.getUntrackedParameter<double>("pfRhoMassMin", 0.0));
  pfRhoMassMax_       = toF(cfg.getUntrackedParameter<double>("pfRhoMassMax", 0.0));

  doTruthMatching_ = cfg.getUntrackedParameter<bool>("doTruthMatching", true);
  matchPhoDR2_     = toF(cfg.getUntrackedParameter<double>("matchPhoDR2", 0.01));
  matchPionDR2_    = toF(cfg.getUntrackedParameter<double>("matchPionDR2", 0.01));

  jetMinPt_     = toF(cfg.getUntrackedParameter<double>("jetMinPt", 30.0));
  htJetMinPt_   = toF(cfg.getUntrackedParameter<double>("htJetMinPt", 30.0));
  jetMaxAbsEta_ = toF(cfg.getUntrackedParameter<double>("jetMaxAbsEta", 4.7));

  // L1 config (defaults match your available products: caloStage2Digis, gmtStage2Digis)
  doL1_          = cfg.getUntrackedParameter<bool>("doL1", true);
  l1UseBX0Only_  = cfg.getUntrackedParameter<bool>("l1UseBX0Only", true);

  l1EGTag_  = cfg.getUntrackedParameter<edm::InputTag>("l1EGamma", edm::InputTag("caloStage2Digis", "EGamma"));
  l1JetTag_ = cfg.getUntrackedParameter<edm::InputTag>("l1Jet",    edm::InputTag("caloStage2Digis", "Jet"));

  l1EGMinPt_      = toF(cfg.getUntrackedParameter<double>("l1EGMinPt", 0.0));
  l1JetMinPt_     = toF(cfg.getUntrackedParameter<double>("l1JetMinPt", 0.0));
  l1EGMaxAbsEta_  = toF(cfg.getUntrackedParameter<double>("l1EGMaxAbsEta", 0.0));
  l1JetMaxAbsEta_ = toF(cfg.getUntrackedParameter<double>("l1JetMaxAbsEta", 0.0));

  matchL1EG_DR2_  = toF(cfg.getUntrackedParameter<double>("matchL1EG_DR2", 0.04));

  // consumes
  puTok_        = consumes<std::vector<PileupSummaryInfo>>(puTag_);
  prunedGenTok_ = consumes<std::vector<reco::GenParticle>>(prunedGenTag_);
  packedGenTok_ = consumes<std::vector<pat::PackedGenParticle>>(packedGenTag_);

  scoutPhoTok_  = consumes<std::vector<Run3ScoutingPhoton>>(phoTag_);
  scoutPFCTok_  = consumes<std::vector<Run3ScoutingParticle>>(pfcTag_);
  scoutJetTok_  = consumes<std::vector<Run3ScoutingPFJet>>(jetTag_);

  l1EGTok_      = consumes<l1t::EGammaBxCollection>(l1EGTag_);
  l1JetTok_     = consumes<l1t::JetBxCollection>(l1JetTag_);
}

void ScoutingHiggsRho::beginJob() {
  if (!storeHistograms_) return;

  // evt
  {
    TFileDirectory dEvt = fs_->mkdir("evt");
    h_evt_nTrueInt_ = dEvt.make<TH1F>("nTrueInt", "nTrueInt;nTrueInt;events", 120, 0., 120.);
  }

  // gen
  {
    TFileDirectory dGen = fs_->mkdir("gen");
    TFileDirectory dAll = dGen.mkdir("all");
    TFileDirectory dFid = dGen.mkdir("fid_strict");
    genAll_.book(dAll, "all");
    genFid_.book(dFid, "fid");
  }

  // scouting
  {
    TFileDirectory dScout = fs_->mkdir("scout");
    TFileDirectory dIncl  = dScout.mkdir("inclusive");
    TFileDirectory dBest  = dScout.mkdir("bestCand");
    TFileDirectory dTM    = dScout.mkdir("truthMatched");
    scoutIncl_.book(dIncl);
    scoutBest_.book(dBest);
    scoutTruth_.book(dTM);
  }

  // matching
  {
    TFileDirectory dMatch = fs_->mkdir("matching");
    match_.book(dMatch);
  }

  // jets
  {
    TFileDirectory dJets = fs_->mkdir("jets");
    jets_.book(dJets);
  }

  // L1
  {
    TFileDirectory dL1 = fs_->mkdir("l1");
    l1_.book(dL1);
  }
}

void ScoutingHiggsRho::analyze(const edm::Event& ev, const edm::EventSetup&) {
  // ----------------------------
  // 0) event-level (safe-ish: data will just have invalid PU handle)
  // ----------------------------
  {
    edm::Handle<std::vector<PileupSummaryInfo>> hPU;
    ev.getByToken(puTok_, hPU);
    if (hPU.isValid() && h_evt_nTrueInt_) {
      for (const auto& pu : *hPU) {
        if (pu.getBunchCrossing() == 0) { h_evt_nTrueInt_->Fill(toF(pu.getTrueNumInteractions())); break; }
      }
    }
  }

  // ----------------------------
  // 1) GEN (MC only)
  // ----------------------------
  bool hasH=false, hasG=false, hasPiPi=false, hasHcand=false;
  const reco::GenParticle* H=nullptr;
  SimpleCand g; RhoCand rho; HCand hc;
  std::vector<SimpleCand> genPionsAll;

  bool f_hasH=false, f_hasG=false, f_hasPiPi=false, f_hasHcand=false;
  const reco::GenParticle* f_H=nullptr;
  SimpleCand f_g; RhoCand f_rho; HCand f_hc;
  std::vector<SimpleCand> genPionsFid;

  if (isMC_) {
    buildGenCandidates_(ev, false, hasH, hasG, hasPiPi, hasHcand, H, g, rho, hc, &genPionsAll);

    // gen/all: keep inclusive Higgs if found; others if present
    if (hasH && H) genAll_.fillH(*H);
    genAll_.fillPions(genPionsAll);
    if (hasG)     genAll_.fillGamma(g);
    if (hasPiPi)  genAll_.fillRho(rho);
    if (hasHcand) genAll_.fillHcand(hc);

    // gen/fid_strict: fill ONLY if hcand exists in fid selection
    buildGenCandidates_(ev, true, f_hasH, f_hasG, f_hasPiPi, f_hasHcand, f_H, f_g, f_rho, f_hc, &genPionsFid);
    if (f_hasHcand) {
      if (f_hasH && f_H) genFid_.fillH(*f_H);
      genFid_.fillPions(genPionsFid);
      genFid_.fillGamma(f_g);
      genFid_.fillRho(f_rho);
      genFid_.fillHcand(f_hc);
    }
  }

  // ----------------------------
  // 2) SCOUT handles (always try)
  // ----------------------------
  edm::Handle<std::vector<Run3ScoutingPhoton>>   hPho;
  edm::Handle<std::vector<Run3ScoutingParticle>> hPFC;
  edm::Handle<std::vector<Run3ScoutingPFJet>>    hJet;
  ev.getByToken(scoutPhoTok_, hPho);
  ev.getByToken(scoutPFCTok_, hPFC);
  ev.getByToken(scoutJetTok_, hJet);

  // Build inclusive rho from PFCands (best OS pair)
  bool hasSRho=false; RhoCand sRho; int sIdxP=-1, sIdxM=-1;
  std::vector<SimpleCand> scoutPionsAll;
  if (hPFC.isValid()) {
    scoutPionsAll.reserve(hPFC->size());
    if (bestPFCandPairToRho_(*hPFC, sIdxP, sIdxM, sRho, &scoutPionsAll)) {
      hasSRho = true;
      scoutIncl_.fillPions(scoutPionsAll);
      scoutIncl_.fillRho(sRho);
    } else {
      scoutIncl_.fillPions(scoutPionsAll);
    }
  }

  // Photon (prefer conditioned on rho)
  bool hasSG=false; SimpleCand sG;
  if (hPho.isValid()) {
    if (hasSRho) {
      if (selectBestScoutPhotonGivenRho_(*hPho, sRho.p4, sG)) { hasSG = true; scoutIncl_.fillGamma(sG); }
      else if (selectBestScoutPhoton_(*hPho, sG))            { hasSG = true; scoutIncl_.fillGamma(sG); }
    } else {
      if (selectBestScoutPhoton_(*hPho, sG)) { hasSG = true; scoutIncl_.fillGamma(sG); }
    }
  }

  // Inclusive hcand (from best rho + chosen gamma)
  bool hasSH=false; HCand sH;
  if (hasSG && hasSRho) {
    TLorentzVector vg; setP4_massless(vg, sG.pt, sG.eta, sG.phi);
    sH.gamma = sG;
    sH.rho = sRho;
    sH.p4 = vg + sRho.p4;
    sH.mass = toF(sH.p4.M());
    hasSH = true;
    scoutIncl_.fillHcand(sH);
  }

  // ----------------------------
  // 2A) Best-candidate stage (kinematic + vertex-ish cuts) using the SAME best rho pair
  // ----------------------------
  if (hasSH && hPFC.isValid() && sIdxP >= 0 && sIdxM >= 0) {
    const auto& piP = hPFC->at(sIdxP);
    const auto& piM = hPFC->at(sIdxM);

    const bool sameV = (piP.vertex() == piM.vertex());
    const float dzDiff  = std::abs(toF(piP.dz())  - toF(piM.dz()));
    const float dxyDiff = std::abs(toF(piP.dxy()) - toF(piM.dxy()));
    const float dr2pp   = deltaR2(sRho.piPlus.eta, sRho.piPlus.phi, sRho.piMinus.eta, sRho.piMinus.phi);
    const float dphiGR  = absDeltaPhi(sG.phi, sRho.phi());
    const float ptbal   = (sG.pt > 0.f) ? (sRho.pt() / sG.pt) : 0.f;

    bool passVtx = true;
    if (pfRequireSameVtx_ && !sameV) passVtx = false;
    if (dzDiff  > toF(pfMaxAbsDzDiff_))  passVtx = false;
    if (dxyDiff > toF(pfMaxAbsDxyDiff_)) passVtx = false;

    bool passKin = true;
    if (dr2pp > toF(pfMaxDR2PiPi_)) passKin = false;
    if (sRho.pt() < toF(pfRhoMinPt_)) passKin = false;
    if (dphiGR < toF(pfMinAbsDPhiGRho_)) passKin = false;
    if (ptbal < toF(pfMinPtBal_) || ptbal > toF(pfMaxPtBal_)) passKin = false;

    if (pfRhoMassMax_ > 0.0) {
      if (sRho.mass < toF(pfRhoMassMin_) || sRho.mass > toF(pfRhoMassMax_)) passKin = false;
    }

    if (passVtx && passKin) {
      scoutBest_.fill(sH);
    }
  }

  // ----------------------------
  // 3A) Truth matching (MC only): GEN gamma + GEN pions -> SCOUT gamma + SCOUT PFCands
  // ----------------------------
  if (isMC_ && doTruthMatching_ && hPho.isValid() && hPFC.isValid()) {
    // denominators
    if (hasG && match_.h_eff_g_den) match_.h_eff_g_den->Fill(g.pt);
    if (hasPiPi && match_.h_eff_rho_den) match_.h_eff_rho_den->Fill(rho.pt());
    if (f_hasHcand) {
      if (match_.h_eff_g_den_fid)   match_.h_eff_g_den_fid->Fill(f_g.pt);
      if (match_.h_eff_rho_den_fid) match_.h_eff_rho_den_fid->Fill(f_rho.pt());
    }

    // match GEN gamma -> SCOUT photon
    bool matchedG=false; SimpleCand mG;
    {
      int best=-1; float bestDR2=matchPhoDR2_;
      for (int i=0;i<(int)hPho->size();++i) {
        const auto& p=hPho->at(i);
        const float dr2 = deltaR2(g.eta, g.phi, toF(p.eta()), toF(p.phi()));
        if (dr2 < bestDR2) { bestDR2=dr2; best=i; }
      }
      if (match_.h_dr2_g && bestDR2 < 998.f) match_.h_dr2_g->Fill(bestDR2);
      if (best >= 0) {
        const auto& p=hPho->at(best);
        mG = SimpleCand{toF(p.pt()), toF(p.eta()), toF(p.phi()), 0.f, 0};
        matchedG = true;
        if (match_.h_eff_g_num) match_.h_eff_g_num->Fill(g.pt);
        if (match_.h_resp_g_pt && g.pt > 0.f) match_.h_resp_g_pt->Fill(mG.pt / g.pt);
        if (f_hasHcand && match_.h_eff_g_num_fid) match_.h_eff_g_num_fid->Fill(f_g.pt);
      }
    }

    // match GEN pi+ / pi- -> SCOUT pfCand (independent best matches)
    bool matchedPiP=false, matchedPiM=false;
    int mPiP=-1, mPiM=-1;
    float dr2_p=999.f, dr2_m=999.f;

    auto matchOnePion = [&](const SimpleCand& genPi, int wantCharge, bool& ok, int& outIdx, float& outDR2) {
      ok = false; outIdx = -1; outDR2 = 999.f;
      int best = -1; float bestDR2 = matchPionDR2_;
      for (int i=0;i<(int)hPFC->size();++i) {
        const auto& p = hPFC->at(i);
        if (pfPionOnly_ && std::abs(p.pdgId()) != 211) continue;
        if (!passScoutPFQuality_(p)) continue;
        const int q = chargeFromPdgId(p.pdgId());
        if (q != wantCharge) continue;
        const float eta = pf_eta(p, pfUseTrkVars_);
        const float phi = pf_phi(p, pfUseTrkVars_);
        const float dr2 = deltaR2(genPi.eta, genPi.phi, eta, phi);
        if (dr2 < bestDR2) { bestDR2 = dr2; best = i; }
      }
      outDR2 = bestDR2;
      if (best >= 0) { ok = true; outIdx = best; }
    };

    bool matchedRho=false; RhoCand mRho;
    if (hasPiPi) {
      matchOnePion(rho.piPlus,  +1, matchedPiP, mPiP, dr2_p);
      matchOnePion(rho.piMinus, -1, matchedPiM, mPiM, dr2_m);
      if (match_.h_dr2_piP && dr2_p < 998.f) match_.h_dr2_piP->Fill(dr2_p);
      if (match_.h_dr2_piM && dr2_m < 998.f) match_.h_dr2_piM->Fill(dr2_m);

      if (matchedPiP && matchedPiM) {
        const auto& piP = hPFC->at(mPiP);
        const auto& piM = hPFC->at(mPiM);

        const float ptP = pf_pt(piP, pfUseTrkVars_), etaP = pf_eta(piP, pfUseTrkVars_), phiP = pf_phi(piP, pfUseTrkVars_);
        const float ptM = pf_pt(piM, pfUseTrkVars_), etaM = pf_eta(piM, pfUseTrkVars_), phiM = pf_phi(piM, pfUseTrkVars_);

        TLorentzVector vP, vM;
        setP4_pion(vP, ptP, etaP, phiP);
        setP4_pion(vM, ptM, etaM, phiM);
        const TLorentzVector vR = vP + vM;

        mRho.piPlus  = SimpleCand{ptP, etaP, phiP, M_PI_CH, +1};
        mRho.piMinus = SimpleCand{ptM, etaM, phiM, M_PI_CH, -1};
        mRho.p4 = vR;
        mRho.mass = toF(vR.M());
        matchedRho = true;

        if (match_.h_eff_rho_num) match_.h_eff_rho_num->Fill(rho.pt());
        if (match_.h_resp_rho_m)  match_.h_resp_rho_m->Fill(mRho.mass - rho.mass);
        if (f_hasHcand && match_.h_eff_rho_num_fid) match_.h_eff_rho_num_fid->Fill(f_rho.pt());
      }
    }

    // build matched hcand and fill scoutTruth_
    if (matchedG && matchedRho) {
      TLorentzVector vG; setP4_massless(vG, mG.pt, mG.eta, mG.phi);
      HCand mH;
      mH.gamma = mG;
      mH.rho = mRho;
      mH.p4 = vG + mRho.p4;
      mH.mass = toF(mH.p4.M());

      scoutTruth_.fillGamma(mG);
      scoutTruth_.fillRho(mRho);
      scoutTruth_.fillHcand(mH);

      if (hasHcand && match_.h_resp_h_m) match_.h_resp_h_m->Fill(mH.mass - hc.mass);
    }
  }

  // ----------------------------
  // 3B) L1 Stage-2 (RECO): EGamma + Jet (no rho/pions, just dumb EG+Jet)
  // ----------------------------
  if (doL1_) {
    edm::Handle<l1t::EGammaBxCollection> hL1EG;
    edm::Handle<l1t::JetBxCollection>    hL1J;
    ev.getByToken(l1EGTok_, hL1EG);
    ev.getByToken(l1JetTok_, hL1J);

    // We only use BX=0 by default (most meaningful offline)
    std::vector<const l1t::EGamma*> egs;
    std::vector<const l1t::Jet*>    jets;
    egs.reserve(32);
    jets.reserve(64);

    if (hL1EG.isValid()) {
      for (int bx = hL1EG->getFirstBX(); bx <= hL1EG->getLastBX(); ++bx) {
        if (l1UseBX0Only_ && bx != 0) continue;
        for (auto it = hL1EG->begin(bx); it != hL1EG->end(bx); ++it) {
          const auto& eg = *it;
          const float pt  = toF(eg.pt());
          const float eta = toF(eg.eta());
          if (pt < toF(l1EGMinPt_)) continue;
          if (l1EGMaxAbsEta_ > 0.0 && std::abs(eta) > toF(l1EGMaxAbsEta_)) continue;
          egs.push_back(&eg);

          if (l1_.h_eg_pt)  l1_.h_eg_pt->Fill(pt);
          if (l1_.h_eg_eta) l1_.h_eg_eta->Fill(eta);
          if (l1_.h_eg_phi) l1_.h_eg_phi->Fill(toF(eg.phi()));
        }
      }
    }

    if (hL1J.isValid()) {
      for (int bx = hL1J->getFirstBX(); bx <= hL1J->getLastBX(); ++bx) {
        if (l1UseBX0Only_ && bx != 0) continue;
        for (auto it = hL1J->begin(bx); it != hL1J->end(bx); ++it) {
          const auto& j = *it;
          const float pt  = toF(j.pt());
          const float eta = toF(j.eta());
          if (pt < toF(l1JetMinPt_)) continue;
          if (l1JetMaxAbsEta_ > 0.0 && std::abs(eta) > toF(l1JetMaxAbsEta_)) continue;
          jets.push_back(&j);

          if (l1_.h_j_pt)  l1_.h_j_pt->Fill(pt);
          if (l1_.h_j_eta) l1_.h_j_eta->Fill(eta);
          if (l1_.h_j_phi) l1_.h_j_phi->Fill(toF(j.phi()));
        }
      }
    }

    if (l1_.h_neg)  l1_.h_neg->Fill(toF((int)egs.size()));
    if (l1_.h_njet) l1_.h_njet->Fill(toF((int)jets.size()));

    // Best EG and best Jet by pt (BX0)
    const l1t::EGamma* bestEG = nullptr;
    const l1t::Jet*    bestJ  = nullptr;

    for (auto* eg : egs) {
      if (!bestEG || eg->pt() > bestEG->pt()) bestEG = eg;
    }
    for (auto* j : jets) {
      if (!bestJ || j->pt() > bestJ->pt()) bestJ = j;
    }

    if (bestEG && bestJ) {
      const float eg_pt  = toF(bestEG->pt());
      const float eg_eta = toF(bestEG->eta());
      const float eg_phi = toF(bestEG->phi());

      const float j_pt  = toF(bestJ->pt());
      const float j_eta = toF(bestJ->eta());
      const float j_phi = toF(bestJ->phi());

      TLorentzVector vEG, vJ;
      setP4_massless(vEG, eg_pt, eg_eta, eg_phi);
      setP4_massless(vJ,  j_pt,  j_eta,  j_phi);
      const TLorentzVector v = vEG + vJ;

      const float dr   = std::sqrt(deltaR2(eg_eta, eg_phi, j_eta, j_phi));
      const float dphi = absDeltaPhi(eg_phi, j_phi);
      const float ptbal = (eg_pt > 0.f) ? (j_pt / eg_pt) : 0.f;

      if (l1_.h_m_egj)      l1_.h_m_egj->Fill(toF(v.M()));
      if (l1_.h_pt_egj)     l1_.h_pt_egj->Fill(toF(v.Pt()));
      if (l1_.h_dr_eg_j)    l1_.h_dr_eg_j->Fill(dr);
      if (l1_.h_dphi_eg_j)  l1_.h_dphi_eg_j->Fill(dphi);
      if (l1_.h_ptbal_j_eg) l1_.h_ptbal_j_eg->Fill(ptbal);

      // Optional: match GEN gamma -> L1 EG (only if MC and GEN gamma exists from your gen/all build)
      if (isMC_ && hasG) {
        const float dr2 = deltaR2(g.eta, g.phi, eg_eta, eg_phi);
        if (l1_.h_dr2_genG_l1EG) l1_.h_dr2_genG_l1EG->Fill(dr2);
        if (dr2 < matchL1EG_DR2_ && l1_.h_resp_eg_pt && g.pt > 0.f) {
          l1_.h_resp_eg_pt->Fill(eg_pt / g.pt);
        }
      }
    }
  }

  // ----------------------------
  // 4) Scouting jets (inclusive diagnostics + toy "jet-as-rho")
  // ----------------------------
  if (hJet.isValid()) {
    std::vector<int> sel;
    sel.reserve(hJet->size());
    double ht = 0.0;

    for (int i=0;i<(int)hJet->size();++i) {
      const auto& j = hJet->at(i);
      const float pt = toF(j.pt());
      const float eta = toF(j.eta());
      if (pt < toF(jetMinPt_)) continue;
      if (std::abs(eta) > toF(jetMaxAbsEta_)) continue;
      sel.push_back(i);
      if (pt >= toF(htJetMinPt_)) ht += pt;
    }

    if (jets_.h_njet) jets_.h_njet->Fill(toF((int)sel.size()));
    if (jets_.h_ht)   jets_.h_ht->Fill(toF(ht));

    if (!sel.empty()) {
      std::sort(sel.begin(), sel.end(), [&](int a, int b){ return toF(hJet->at(a).pt()) > toF(hJet->at(b).pt()); });
      const auto& j1 = hJet->at(sel[0]);
      if (jets_.h_j1_pt)   jets_.h_j1_pt->Fill(toF(j1.pt()));
      if (jets_.h_j1_eta)  jets_.h_j1_eta->Fill(toF(j1.eta()));
      if (jets_.h_j1_phi)  jets_.h_j1_phi->Fill(toF(j1.phi()));
      if (jets_.h_j1_m)    jets_.h_j1_m->Fill(toF(j1.m()));
      if (jets_.h_j1_area) jets_.h_j1_area->Fill(toF(j1.jetArea()));

      // min dr2 to gamma/rho (use the already-chosen sG/sRho)
      if (hasSG && jets_.h_min_dr2_g_j) {
        float best = 1e9f;
        for (int idx : sel) {
          const auto& j = hJet->at(idx);
          best = std::min(best, deltaR2(sG.eta, sG.phi, toF(j.eta()), toF(j.phi())));
        }
        jets_.h_min_dr2_g_j->Fill(best);
      }
      if (hasSRho && jets_.h_min_dr2_rho_j) {
        float best = 1e9f;
        for (int idx : sel) {
          const auto& j = hJet->at(idx);
          best = std::min(best, deltaR2(sRho.eta(), sRho.phi(), toF(j.eta()), toF(j.phi())));
        }
        jets_.h_min_dr2_rho_j->Fill(best);
      }

      // toy: treat leading jet as "rho proxy" and build jet+gamma mass
      if (hasSG) {
        TLorentzVector vJ, vG;
        vJ.SetPtEtaPhiM(toF(j1.pt()), toF(j1.eta()), toF(j1.phi()), toF(j1.m()));
        setP4_massless(vG, sG.pt, sG.eta, sG.phi);
        const TLorentzVector vH = vJ + vG;
        if (jets_.h_h_m_jetAsRho)  jets_.h_h_m_jetAsRho->Fill(toF(vH.M()));
        if (jets_.h_h_pt_jetAsRho) jets_.h_h_pt_jetAsRho->Fill(toF(vH.Pt()));
      }
    }
  }
}

DEFINE_FWK_MODULE(ScoutingHiggsRho);