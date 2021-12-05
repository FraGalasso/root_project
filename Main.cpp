#include <cstring>
#include <iostream>
#include <vector>

#include "Particle.hpp"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom.h"
// #include "TCanvas.h"

R__LOAD_LIBRARY(ParticleType_cpp.so);
R__LOAD_LIBRARY(ResonanceType_cpp.so);
R__LOAD_LIBRARY(Particle2_cpp.so);

int main() {
  std::vector<Particle> particle_vec{{}};
  std::vector<Particle> generated_particle_vec{{}};
  Particle::AddParticleType("pion+", 0.13957, 1);
  Particle::AddParticleType("pion-", 0.13957, -1);
  Particle::AddParticleType("kaon+", 0.49367, 1);
  Particle::AddParticleType("kaon-", 0.49367, -1);
  Particle::AddParticleType("proton", 0.93827, 1);
  Particle::AddParticleType("antiproton", 0.93827, -1);
  Particle::AddParticleType("k*", 0.89166, 0, 0.05);

  Particle::PrintParticleTypes();

  // 0 means pion+, 1 means pion-, 2 means kaon+, 3 means kaon-, 4 means proton,
  // 5 means antiproton, 6 means k*
  std::vector<int> particles_generated = {0, 0, 0, 0, 0, 0, 0};

  // rootfile where I can save my histos
  TFile *file = new TFile("chediocelamandibuona.root", "RECREATE");

  // histo for phi distribution
  TH1F *phi_distrib =
      new TH1F("phi_distrib", "Distribution of phi", 400, 0, 2 * TMath::Pi());

  // histo for theta distribution
  TH1F *theta_distrib =
      new TH1F("theta_distrib", "Distribution of theta", 400, 0, TMath::Pi());

  // histo for momentum distribution
  TH1F *momentum_distrib =
      new TH1F("momentum_distrib", "Distribution of momentum", 1000, 0, 6);

  // histo for transverse momentum distribution
  TH1F *trans_momentum_distrib =
      new TH1F("trans_momentum_distrib", "Distribution of transverse momentum",
               1000, 0, 6);

  // histo for energy distribution
  TH1F *energy_distrib =
      new TH1F("energy_distrib", "Distribution of energy", 10000, 0, 8);

  // histo for invariant mass distribution
  TH1F *inv_mass =
      new TH1F("inv_mass", "Distribution of invariant mass", 10000, 0, 8);
  inv_mass->Sumw2();

  // histo for invariant mass distribution with different sign
  TH1F *diff_inv_mass = new TH1F(
      "diff_inv_mass",
      "Distribution of invariant mass of particles with different signs", 10000,
      0, 8);
  diff_inv_mass->Sumw2();

  // histo for invariant mass distribution with same sign
  TH1F *same_inv_mass =
      new TH1F("same_inv_mass",
               "Distribution of invariant mass of particles with same signs",
               10000, 0, 8);
  same_inv_mass->Sumw2();

  // histo for invariant mass of kaons and pions with same sign
  TH1F *same_inv_mass_ka_pi = new TH1F(
      "same_sign_inv_mass_ka_pi",
      "Distribution of invariant mass of kaons and pions with same signs",
      10000, 0, 8);
  same_inv_mass_ka_pi->Sumw2();

  // histo for invariant mass of kaons and pions with different sign
  TH1F *diff_inv_mass_ka_pi = new TH1F(
      "diff_inv_mass_ka_pi",
      "Distribution of invariant mass of kaons and pions with different signs",
      10000, 0, 8);
  diff_inv_mass_ka_pi->Sumw2();

  // histo for invariant mass of particles generated from a decayment
  TH1F *part_gen = new TH1F(
      "part_gen",
      "Distribution of invariant mass of particles generated from a decayment",
      10000, 0, 8);
  part_gen->Sumw2();

  Particle p("name", 0, 0, 0);
  Particle dau1("name", 0, 0, 0);
  Particle dau2("name", 0, 0, 0);
  double phi{0};
  double theta{0};
  double momentum{0};
  double px{0};
  double py{0};
  double pz{0};
  double transverse_momentum{0};

  for (int i = 0; i < 100000; ++i) {
    for (int j = 0; j < 100; ++j) {
      phi = (gRandom->Uniform(2 * TMath::Pi()));
      phi_distrib->Fill(phi);
      theta = (gRandom->Uniform(TMath::Pi()));
      theta_distrib->Fill(theta);
      momentum = (gRandom->Exp(1));
      momentum_distrib->Fill(momentum);
      px = momentum * TMath::Sin(theta) * TMath::Cos(phi);
      py = momentum * TMath::Sin(theta) * TMath::Sin(phi);
      pz = momentum * TMath::Cos(theta);
      p.SetP(px, py, pz);
      transverse_momentum = TMath::Sqrt(px * px + py * py);
      trans_momentum_distrib->Fill(transverse_momentum);

      double x = gRandom->Uniform(1);
      if (x < 0.4) {
        p.SetIndex("pion+");

        ++particles_generated[0];
      } else if (x < 0.8) {
        p.SetIndex("pion-");
        ++particles_generated[1];
      } else if (x < 0.85) {
        p.SetIndex("kaon+");
        ++particles_generated[2];
      } else if (x < 0.9) {
        p.SetIndex("kaon-");
        ++particles_generated[3];
      } else if (x < 0.945) {
        p.SetIndex("proton");
        ++particles_generated[4];
      } else if (x < 0.99) {
        p.SetIndex("antiproton");
        ++particles_generated[5];
      } else {
        p.SetIndex("k*");
        ++particles_generated[6];
        double y = gRandom->Uniform(1);
        if (y < 0.5) {
          dau1.SetIndex("pion+");
          dau2.SetIndex("kaon-");
        } else {
          dau1.SetIndex("pion-");
          dau2.SetIndex("kaon+");
        };
        p.Decay2body(dau1, dau2);
        generated_particle_vec.push_back(dau1);
        generated_particle_vec.push_back(dau1);
        part_gen->Fill(dau1.InvMass(dau2));
      }
      particle_vec.push_back(p);
      energy_distrib->Fill(p.Energy());
    };
    particle_vec.clear();
    generated_particle_vec.clear();
  };
  for (int k = 0; k < int(particle_vec.size()); ++k) {
    if ((std::strcmp(particle_vec[k].GetParticleName(), "k*")) == 0) {
      continue;
    }
    for (int l = k + 1; l < int(particle_vec.size()); ++l) {
      if ((std::strcmp(particle_vec[l].GetParticleName(), "k*")) == 0) {
        continue;
      }
      double m = particle_vec[k].InvMass(particle_vec[l]);
      inv_mass->Fill(m);
      if (particle_vec[k].GetParticleCharge() !=
          particle_vec[l].GetParticleCharge()) {
        diff_sign_inv_mass->Fill(m);
      } else {
        same_sign_inv_mass->Fill(m);
      }
      if ((std::strcmp(particle_vec[k].GetParticleName(), "pion+") == 0 &&
           std::strcmp(particle_vec[k].GetParticleName(), "kaon+") == 0) ||
          (std::strcmp(particle_vec[k].GetParticleName(), "pion-") == 0 &&
           std::strcmp(particle_vec[k].GetParticleName(), "kaon-") == 0)) {
        same_sign_inv_mass_ka_pi->Fill(m);
      }
      if ((std::strcmp(particle_vec[k].GetParticleName(), "pion+") == 0 &&
           std::strcmp(particle_vec[k].GetParticleName(), "kaon-") == 0) ||
          (std::strcmp(particle_vec[k].GetParticleName(), "pion-") == 0 &&
           std::strcmp(particle_vec[k].GetParticleName(), "kaon+") == 0)) {
        diff_sign_inv_mass_ka_pi->Fill(m);
      }
    }
    for (int l = 0; l < int(generated_particle_vec.size()); ++l) {
      if ((std::strcmp(particle_vec[l].GetParticleName(), "k*")) == 0) {
        continue;
      };
      double m = particle_vec[k].InvMass(particle_vec[l]);
      inv_mass->Fill(m);
      if (particle_vec[k].GetParticleCharge() !=
          particle_vec[l].GetParticleCharge()) {
        diff_sign_inv_mass->Fill(m);
      } else {
        same_sign_inv_mass->Fill(m);
      };
      if ((std::strcmp(particle_vec[k].GetParticleName(), "pion+") == 0 &&
           std::strcmp(particle_vec[k].GetParticleName(), "kaon+") == 0) ||
          (std::strcmp(particle_vec[k].GetParticleName(), "pion-") == 0 &&
           std::strcmp(particle_vec[k].GetParticleName(), "kaon-") == 0)) {
        same_sign_inv_mass_ka_pi->Fill(m);
      };
      if ((std::strcmp(particle_vec[k].GetParticleName(), "pion+") == 0 &&
           std::strcmp(particle_vec[k].GetParticleName(), "kaon-") == 0) ||
          (std::strcmp(particle_vec[k].GetParticleName(), "pion-") == 0 &&
           std::strcmp(particle_vec[k].GetParticleName(), "kaon+") == 0)) {
        diff_sign_inv_mass_ka_pi->Fill(m);
      }
    }
    for (int a = 0; a < int(generated_particle_vec.size()); ++a) {
      for (int b = a + 1; b < int(generated_particle_vec.size()); ++b) {
        if ((std::strcmp(generated_particle_vec[b].GetParticleName(), "k*")) ==
            0) {
          continue;
        };
        double mm =
            generated_particle_vec[a].InvMass(generated_particle_vec[b]);
        inv_mass->Fill(mm);
        if (generated_particle_vec[a].GetParticleCharge() !=
            generated_particle_vec[b].GetParticleCharge()) {
          diff_sign_inv_mass->Fill(mm);
        } else {
          same_sign_inv_mass->Fill(mm);
        };
        if ((std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "pion+") == 0 &&
             std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "kaon+") == 0) ||
            (std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "pion-") == 0 &&
             std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "kaon-") == 0)) {
          same_sign_inv_mass_ka_pi->Fill(mm);
        };
        if ((std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "pion+") == 0 &&
             std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "kaon-") == 0) ||
            (std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "pion-") == 0 &&
             std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "kaon+") == 0)) {
          diff_sign_inv_mass_ka_pi->Fill(mm);
        }
      };
    };

    for (int a = 0; a < int(generated_particle_vec.size()); ++a) {
      for (int b = a + 1; b < int(generated_particle_vec.size()); ++b) {
        if ((std::strcmp(generated_particle_vec[b].GetParticleName(), "k*")) ==
            0) {
          continue;
        };
        double mm =
            generated_particle_vec[a].InvMass(generated_particle_vec[b]);
        inv_mass->Fill(mm);
        if (generated_particle_vec[a].GetParticleCharge() !=
            generated_particle_vec[b].GetParticleCharge()) {
          diff_sign_inv_mass->Fill(mm);
        } else {
          same_sign_inv_mass->Fill(mm);
        };
        if ((std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "pion+") == 0 &&
             std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "kaon+") == 0) ||
            (std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "pion-") == 0 &&
             std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "kaon-") == 0)) {
          same_sign_inv_mass_ka_pi->Fill(mm);
        };
        if ((std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "pion+") == 0 &&
             std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "kaon-") == 0) ||
            (std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "pion-") == 0 &&
             std::strcmp(generated_particle_vec[a].GetParticleName(),
                         "kaon+") == 0)) {
          diff_sign_inv_mass_ka_pi->Fill(mm);
        }
      }
    }

    TH1F *histoparticles =
        new TH1F("histoparticles", "Particles generated", 7, 1, 7);

    for (int n = 0; n < 7; ++n) {
      histoparticles->Fill(n + 1, particles_generated[n]);
    };

    TCanvas *can = new TCanvas("can", "Many Histos", 200, 10, 600, 400);

    can->cd(1);
    histoparticles->DrawCopy("H");
    histoparticles->DrawCopy("E, P, same");
    can->cd(2);
    phi_distrib->DrawCopy("H");
    phi_distrib->DrawCopy("E, P, same");
    can->cd(3);
    theta_distrib->DrawCopy("H");
    theta_distrib->DrawCopy("E, P, same");
    can->cd(4);
    momentum_distrib->DrawCopy("H");
    momentum_distrib->DrawCopy("E, P, same");
    can->cd(5);
    trans_momentum_distrib->DrawCopy("H");
    trans_momentum_distrib->DrawCopy("E, P, same");
    can->cd(6);
    energy_distrib->DrawCopy("H");
    energy_distrib->DrawCopy("E, P, same");
    can->cd(7);
    inv_mass->DrawCopy("H");
    inv_mass->DrawCopy("E, P, same");
    can->cd(8);
    diff_sign_inv_mass->DrawCopy("H");
    diff_sign_inv_mass->DrawCopy("E, P, same");
    can->cd(9);
    same_sign_inv_mass->DrawCopy("H");
    same_sign_inv_mass->DrawCopy("E, p, same");
    can->cd(10);
    diff_sign_inv_mass_ka_pi->DrawCopy("H");
    diff_sign_inv_mass_ka_pi->DrawCopy("E, P, same");
    can->cd(11);
    same_sign_inv_mass_ka_pi->DrawCopy("H");
    same_sign_inv_mass_ka_pi->DrawCopy("E, P, same");
    can->cd(12);
    part_gen->DrawCopy("H");
    part_gen->DrawCopy("E, P, same");

    file->Write();
    file->Close();
  }
}
