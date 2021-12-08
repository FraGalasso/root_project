#include <cstring>
#include <iostream>
#include <vector>

#include "Particle.hpp"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TRandom.h"

R__LOAD_LIBRARY(ParticleType_cpp.so);
R__LOAD_LIBRARY(ResonanceType_cpp.so);
R__LOAD_LIBRARY(Particle_cpp.so);

int main() {
  // vectors where I will store my particles
  std::vector<Particle> particle_vec{};
  std::vector<Particle> generated_particle_vec{};

  // initializing different particle types
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
               1000, 0, 4);

  // histo for energy distribution
  TH1F *energy_distrib =
      new TH1F("energy_distrib", "Distribution of energy", 10000, 0, 5);

  // histo for invariant mass distribution
  TH1F *inv_mass =
      new TH1F("inv_mass", "Distribution of invariant mass", 10000, 0, 5);
  inv_mass->Sumw2();

  // histo for invariant mass distribution with different sign
  TH1F *diff_inv_mass = new TH1F(
      "diff_inv_mass",
      "Distribution of invariant mass of particles with different signs", 10000,
      0, 5);
  diff_inv_mass->Sumw2();

  // histo for invariant mass distribution with same sign
  TH1F *same_inv_mass =
      new TH1F("same_inv_mass",
               "Distribution of invariant mass of particles with same signs",
               10000, 0, 5);
  same_inv_mass->Sumw2();

  // histo for invariant mass of kaons and pions with same sign
  TH1F *same_inv_mass_ka_pi = new TH1F(
      "same_inv_mass_ka_pi",
      "Distribution of invariant mass of kaons and pions with same signs",
      10000, 0, 5);
  same_inv_mass_ka_pi->Sumw2();

  // histo for invariant mass of kaons and pions with different sign
  TH1F *diff_inv_mass_ka_pi = new TH1F(
      "diff_inv_mass_ka_pi",
      "Distribution of invariant mass of kaons and pions with different signs",
      10000, 0, 5);
  diff_inv_mass_ka_pi->Sumw2();

  // histo for invariant mass of particles generated from a decayment
  TH1F *part_gen = new TH1F(
      "part_gen",
      "Distribution of invariant mass of particles generated from a decayment",
      10000, 0, 5);
  part_gen->Sumw2();

  /* inizializing particles as protons, just to avoid errors from using
     particles with uninitialized names, I will change them later on*/
  Particle p("proton", 0, 0, 0);
  Particle dau1("proton", 0, 0, 0);
  Particle dau2("proton", 0, 0, 0);
  double phi{0};
  double theta{0};
  double momentum{0};
  double px{0};
  double py{0};
  double pz{0};
  double transverse_momentum{0};
  gRandom->SetSeed();

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

      // generating a random particle, proportions are:
      // 80% pions, 10% kaons, 9% protons, 1% k*(resonance)
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
        // if a particle is a k* it will decay, proportions are:
        // 50% pion+/kaon-, 50% pion-/kaon+
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
        generated_particle_vec.push_back(dau2);
        part_gen->Fill(dau1.InvMass(dau2));
      }
      particle_vec.push_back(p);
      energy_distrib->Fill(p.Energy());
    }
    // filling up all invariant mass histos
    for (int k = 0; k < int(particle_vec.size()); ++k) {
      if ((std::strcmp(particle_vec[k].GetParticleName(), "k*") == 0)) {
        continue;
      }
      for (int l = k + 1; l < int(particle_vec.size()); ++l) {
        if ((std::strcmp(particle_vec[l].GetParticleName(), "k*")) == 0) {
          continue;
        }
        double const m = particle_vec[k].InvMass(particle_vec[l]);
        inv_mass->Fill(m);
        if (particle_vec[k].GetParticleCharge() !=
            particle_vec[l].GetParticleCharge()) {
          diff_inv_mass->Fill(m);
        } else {
          same_inv_mass->Fill(m);
        }
        if (std::strcmp(particle_vec[k].GetParticleName(), "pion+") == 0) {
          if (std::strcmp(particle_vec[l].GetParticleName(), "kaon+") == 0) {
            same_inv_mass_ka_pi->Fill(m);
          }
          if (std::strcmp(particle_vec[l].GetParticleName(), "kaon-")) {
            diff_inv_mass_ka_pi->Fill(m);
          }
        }
        if (std::strcmp(particle_vec[k].GetParticleName(), "pion-") == 0) {
          if (std::strcmp(particle_vec[l].GetParticleName(), "kaon-") == 0) {
            same_inv_mass_ka_pi->Fill(m);
          }
          if (std::strcmp(particle_vec[l].GetParticleName(), "kaon+") == 0) {
            diff_inv_mass_ka_pi->Fill(m);
          }
        }
        if (std::strcmp(particle_vec[k].GetParticleName(), "kaon+") == 0) {
          if (std::strcmp(particle_vec[l].GetParticleName(), "pion+") == 0) {
            same_inv_mass_ka_pi->Fill(m);
          }
          if (std::strcmp(particle_vec[l].GetParticleName(), "pion-") == 0) {
            diff_inv_mass_ka_pi->Fill(m);
          }
        }
        if (std::strcmp(particle_vec[k].GetParticleName(), "kaon-") == 0) {
          if (std::strcmp(particle_vec[l].GetParticleName(), "pion-") == 0) {
            same_inv_mass_ka_pi->Fill(m);
          }
          if (std::strcmp(particle_vec[l].GetParticleName(), "pion+") == 0) {
            diff_inv_mass_ka_pi->Fill(m);
          }
        }
      }
      for (int f = 0; f < int(generated_particle_vec.size()); ++f) {
        double const m = particle_vec[k].InvMass(generated_particle_vec[f]);
        inv_mass->Fill(m);
        if (particle_vec[k].GetParticleCharge() !=
            particle_vec[f].GetParticleCharge()) {
          diff_inv_mass->Fill(m);
        } else {
          same_inv_mass->Fill(m);
        };

        if (std::strcmp(particle_vec[k].GetParticleName(), "pion+") == 0) {
          if (std::strcmp(generated_particle_vec[f].GetParticleName(),
                          "kaon+") == 0) {
            same_inv_mass_ka_pi->Fill(m);
          }
          if (std::strcmp(generated_particle_vec[f].GetParticleName(),
                          "kaon-") == 0) {
            diff_inv_mass_ka_pi->Fill(m);
          }
        }
        if (std::strcmp(particle_vec[k].GetParticleName(), "pion-") == 0) {
          if (std::strcmp(generated_particle_vec[f].GetParticleName(),
                          "kaon-") == 0) {
            same_inv_mass_ka_pi->Fill(m);
          }
          if (std::strcmp(generated_particle_vec[f].GetParticleName(),
                          "kaon+") == 0) {
            diff_inv_mass_ka_pi->Fill(m);
          }
        }
        if (std::strcmp(particle_vec[k].GetParticleName(), "kaon+") == 0) {
          if (std::strcmp(generated_particle_vec[f].GetParticleName(),
                          "pion+") == 0) {
            same_inv_mass_ka_pi->Fill(m);
          }
          if (std::strcmp(generated_particle_vec[f].GetParticleName(),
                          "pion-") == 0) {
            diff_inv_mass_ka_pi->Fill(m);
          }
        }
        if (std::strcmp(particle_vec[k].GetParticleName(), "kaon-") == 0) {
          if (std::strcmp(generated_particle_vec[f].GetParticleName(),
                          "pion-") == 0) {
            same_inv_mass_ka_pi->Fill(m);
          }
          if (std::strcmp(generated_particle_vec[f].GetParticleName(),
                          "pion+") == 0) {
            diff_inv_mass_ka_pi->Fill(m);
          }
        }
      }
    }
    for (int a = 0; a < int(generated_particle_vec.size()); ++a) {
      for (int b = a + 1; b < int(generated_particle_vec.size()); ++b) {
        double const m =
            generated_particle_vec[a].InvMass(generated_particle_vec[b]);
        inv_mass->Fill(m);
        if (generated_particle_vec[a].GetParticleCharge() !=
            generated_particle_vec[b].GetParticleCharge()) {
          diff_inv_mass->Fill(m);
        } else {
          same_inv_mass->Fill(m);
        }

        if (std::strcmp(generated_particle_vec[a].GetParticleName(), "pion+") ==
            0) {
          if (std::strcmp(generated_particle_vec[b].GetParticleName(),
                          "kaon+") == 0) {
            same_inv_mass_ka_pi->Fill(m);
          }
          if (std::strcmp(generated_particle_vec[b].GetParticleName(),
                          "kaon-") == 0) {
            diff_inv_mass_ka_pi->Fill(m);
          }
        }
        if (std::strcmp(generated_particle_vec[a].GetParticleName(), "pion-") ==
            0) {
          if (std::strcmp(generated_particle_vec[b].GetParticleName(),
                          "kaon-") == 0) {
            same_inv_mass_ka_pi->Fill(m);
          }
          if (std::strcmp(generated_particle_vec[b].GetParticleName(),
                          "kaon+") == 0) {
            diff_inv_mass_ka_pi->Fill(m);
          }
        }
        if (std::strcmp(generated_particle_vec[a].GetParticleName(), "kaon+") ==
            0) {
          if (std::strcmp(generated_particle_vec[b].GetParticleName(),
                          "pion+") == 0) {
            same_inv_mass_ka_pi->Fill(m);
          }
          if (std::strcmp(generated_particle_vec[b].GetParticleName(),
                          "pion-") == 0) {
            diff_inv_mass_ka_pi->Fill(m);
          }
        }
        if (std::strcmp(generated_particle_vec[a].GetParticleName(), "kaon-") ==
            0) {
          if (std::strcmp(generated_particle_vec[b].GetParticleName(),
                          "pion-") == 0) {
            same_inv_mass_ka_pi->Fill(m);
          }
          if (std::strcmp(generated_particle_vec[b].GetParticleName(),
                          "pion+") == 0) {
            diff_inv_mass_ka_pi->Fill(m);
          }
        }
      }
    }
    /* when I'm done with all the invariant mass stuff, I have to clear the 2
       vectors and start again*/
    particle_vec.clear();
    generated_particle_vec.clear();
  }

  // histo for particle type distribution
  TH1F *histoparticles =
      new TH1F("histoparticles", "Particles generated", 7, 1, 7);
  for (int n = 0; n < int(particles_generated.size()); ++n) {
    histoparticles->SetBinContent(n + 1, (particles_generated[n]) / total);
    std::cout << "Particle Number: " << n + 1
              << " | Entries: " << particles_generated[n] << '\n';
  };

  // rootfile where I can save my histos
  TFile *file = new TFile("analysis.root", "RECREATE");

  histoparticles->Write("histoparticles");
  phi_distrib->Write("phi_distrib");
  theta_distrib->Write("theta_distrib");
  momentum_distrib->Write("momentum_distrib");
  trans_momentum_distrib->Write("trans_momentum_distrib");
  energy_distrib->Write("energy_distrib");
  inv_mass->Write("inv_mass");
  diff_inv_mass->Write("diff_inv_mass");
  same_inv_mass->Write("same_inv_mass");
  diff_inv_mass_ka_pi->Write("diff_inv_mass_ka_pi");
  same_inv_mass_ka_pi->Write("same_inv_mass_ka_pi");
  part_gen->Write("part_gen");

  file->Close();

  delete histoparticles;
  delete phi_distrib;
  delete theta_distrib;
  delete momentum_distrib;
  delete trans_momentum_distrib;
  delete energy_distrib;
  delete inv_mass;
  delete diff_inv_mass;
  delete same_inv_mass;
  delete diff_inv_mass_ka_pi;
  delete same_inv_mass_ka_pi;
  delete part_gen;
  delete file;
}
