#include "Particle.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "ParticleType.hpp"
#include "ResonanceType.hpp"

std::vector<ParticleType *> Particle::fParticleType;

// ParticleType *Particle::fParticleType[fMaxNumParticleType];

int Particle::fNParticleType{0};

Particle::Particle() : fIndex(-1), p_x_(0), p_y_(0), p_z_(0){};
Particle::Particle(const char *name, double p_x, double p_y, double p_z) {
  p_x_ = p_x;
  p_y_ = p_y;
  p_z_ = p_z;
  const int index = FindParticle(name);
  if (index == -1) {
    std::cout << "There is no particle with name " << name << ".\n";
  };
  fIndex = index;
};

int Particle::FindParticle(const char *name) {
  for (int i = 0; i < fNParticleType; ++i) {
    if (std::strcmp(fParticleType[i]->GetName(), name) == 0) {
      return i;  // controllo se ci sono particelle con questo nome
    };
  };
  return -1;  // se alla fine del ciclo non ne ho trovate ritorno un errore
};

int Particle::GetIndex() { return fIndex; };
double Particle::GetPx() const { return p_x_; };
double Particle::GetPy() const { return p_y_; };
double Particle::GetPz() const { return p_z_; };
double Particle::GetParticleMass() const {
  return fParticleType[fIndex]->GetMass();
};
int Particle::GetParticleCharge() const {
  return fParticleType[fIndex]->GetCharge();
};
const char *Particle::GetParticleName() {
  return fParticleType[fIndex]->GetName();
};

void Particle::AddParticleType(const char *name, double mass, int charge,
                               double width) {
  if (fNParticleType < fMaxNumParticleType) {
    if (FindParticle(name) == -1) {
      if (width == 0) {
        // fParticleType[fNParticleType] = new ParticleType{name, mass, charge};
        ParticleType *pParticle = new ParticleType{name, mass, charge};
        fParticleType.push_back(pParticle);
        ++fNParticleType;
      } else {
        // fParticleType[fNParticleType] = new ResonanceType{name, mass, charge,
        // width};
        ResonanceType *rParticle = new ResonanceType{name, mass, charge, width};
        fParticleType.push_back(rParticle);
        ++fNParticleType;
      }
    } else {
      std::cout << "This particle already is in the ParticeType array.\n";
    }
  } else {
    std::cerr << "I don't know how you did this, but it is extremely wrong";
  }
};

int Particle::Decay2body(Particle &dau1, Particle &dau2) const {
  if (GetParticleMass() == 0.0) {
    std::cout << "Decayment cannot be preformed if mass is zero\n";
    return 1;
  }

  double massMot = GetParticleMass();
  double massDau1 = dau1.GetParticleMass();
  double massDau2 = dau2.GetParticleMass();

  if (fIndex > -1) {
    float x1, x2, w, y1, y2;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);
    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;
    massMot += fParticleType[fIndex]->GetWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    std::cout << "Decayment cannot be preformed because mass is too low in "
                 "this channel\n";
    return 2;
  }

  double pout =
      sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  dau1.SetP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi),
            pout * cos(theta));
  dau2.SetP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi),
            -pout * cos(theta));

  double energy =
      sqrt(p_x_ * p_x_ + p_y_ * p_y_ + p_z_ * p_z_ + massMot * massMot);

  double bx = p_x_ / energy;
  double by = p_y_ / energy;
  double bz = p_z_ / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
};

void Particle::Boost(double bx, double by, double bz) {
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx * p_x_ + by * p_y_ + bz * p_z_;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  p_x_ += gamma2 * bp * bx + gamma * bx * Energy();
  p_y_ += gamma2 * bp * by + gamma * by * Energy();
  p_z_ += gamma2 * bp * bz + gamma * bz * Energy();
};

double Particle::Momentum2() const {
  return p_x_ * p_x_ + p_y_ * p_y_ + p_z_ * p_z_;
};

double Particle::Momentum() const { return sqrt(Momentum2()); };

double Particle::Energy() const {
  return sqrt(GetParticleMass() * GetParticleMass() + Momentum2());
};

double Particle::InvMass(Particle &p) const {
  double const Energy2 = pow(Energy() + p.Energy(), 2);
  double const Mom2 = pow(Momentum() + p.Momentum(), 2);
  /*ATTENZIONE, controllare se non devo fare la somma vettoriale*/

  return sqrt(Energy2 - Mom2);
};

void Particle::SetIndex(int index) {
  if (index < 0 || index >= fNParticleType) {
    std::cout << "Invalid particle index\n";
    return;
  }

  fIndex = index;
};

void Particle::SetIndex(const char *name) { SetIndex(FindParticle(name)); };

void Particle::SetP(double px, double py, double pz) {
  p_x_ = px;
  p_y_ = py;
  p_z_ = pz;
};

void Particle::PrintParticleTypes() {
  for (int i = 0; i < fNParticleType; ++i) {
    std::cout << "Particle Type: " << i + 1;
    fParticleType[i]->Print();
    std::cout << '\n';
  }
};

void Particle::PrintParticle() {
  std::cout << "Type index: " << fIndex
            << "\nParticle name: " << fParticleType[fIndex]->GetName()
            << "\nParticle momentum: (" << p_x_ << ", " << p_y_ << ", " << p_z_
            << ")\n";
};
