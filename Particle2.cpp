#include "Particle2.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>

std::vector<ParticleType *> Particle::fParticleType;

Particle::Particle() : fIndex(-1), p_x_(0), p_y_(0), p_z_(0){};
Particle::Particle(const char *name, double p_x, double p_y, double p_z) {
  p_x_ = p_x;
  p_y_ = p_y;
  p_z_ = p_z;
  fIndex = FindParticle(name);
  if (fIndex == -1) {
    std::cout << "There is no particle with name " << name << ".\n";
  }
};

int Particle::FindParticle(const char *name) {
  for (int i = 0; i < int(fParticleType.size()); ++i) {
    if (std::strcmp(fParticleType.at(i)->GetName(), name) == 0) {
      return i;  // checking whether there are particles with this name
    };
  };
  return -1;  // return an error if particle cannot be found
};

int Particle::GetIndex() { return fIndex; };
double Particle::GetPx() const { return p_x_; };
double Particle::GetPy() const { return p_y_; };
double Particle::GetPz() const { return p_z_; };
double Particle::GetParticleMass() const {
  if (fIndex == -1) {
    std::cout << "You didn't initialize a particle. a\n";
    return 0;
  }
  return fParticleType.at(fIndex)->GetMass();
};
int Particle::GetParticleCharge() const {
  if (fIndex == -1) {
    std::cout << "You didn't initialize a particle. b\n";
    return 0;
  }
  return fParticleType.at(fIndex)->GetCharge();
};
const char *Particle::GetParticleName() {
  // if (fIndex == -1) {
  //   std::cout << "You didn't initialize a particle. c\n";
  //   return " ";
  // }
  return fParticleType.at(fIndex)->GetName();
};

void Particle::AddParticleType(const char *name, double mass, int charge,
                               double width) {
  if (int(fParticleType.size()) >= fMaxNumParticleType) {
    std::cout << "I didn't manage to add " << name
              << ".\nPlease don't add too many particles: only "
              << fMaxNumParticleType << " are allowed.\n";
    return;
  }
  if (FindParticle(name) == -1) {
    if (width == 0) {
      ParticleType *pParticle = new ParticleType{name, mass, charge};
      fParticleType.push_back(pParticle);
    } else {
      ResonanceType *rParticle = new ResonanceType{name, mass, charge, width};
      fParticleType.push_back(rParticle);
    }
  } else {
    std::cout << "This particle already is in the ParticeType array.\n"
              << "Why are you trying to add " << name << " again?\n";
  }
}

void Particle::AddParticleType(const char *name, double mass, int charge) {
  AddParticleType(name, mass, charge, 0.0);
};

int Particle::Decay2body(Particle &dau1, Particle &dau2) const {
  if (GetParticleMass() == 0.0) {
    std::cout << "Decayment cannot be performed if mass is zero\n";
    return 1;
  }

  double massMot = GetParticleMass();
  double massDau1 = dau1.GetParticleMass();
  double massDau2 = dau2.GetParticleMass();

  if (fIndex > -1) {
    float x1, x2, w, y1;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);
    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    massMot += fParticleType.at(fIndex)->GetWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    std::cout << "Decayment cannot be performed because mass is too low in "
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
  if (index < 0 || index >= int(fParticleType.size())) {
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
  for (int i = 0; i < int(fParticleType.size()); ++i) {
    std::cout << "Particle Type: " << i + 1;
    fParticleType.at(i)->Print();
    std::cout << '\n';
  }
};

void Particle::PrintParticle() {
  std::cout << "Type index: " << fIndex
            << "\nParticle name: " << fParticleType.at(fIndex)->GetName()
            << "\nParticle momentum: (" << p_x_ << ", " << p_y_ << ", " << p_z_
            << ")\n";
};
