#ifndef PARTICLE2_HPP
#define PARTICLE2_HPP

#include <vector>

#include "ResonanceType.hpp"

class Particle {
 public:
  Particle();
  Particle(const char *name, double p_x, double p_y, double p_z);

  int GetIndex();
  double GetPx() const;
  double GetPy() const;
  double GetPz() const;
  double GetParticleMass() const;
  int GetParticleCharge() const;
  const char *GetParticleName();

  /*adds new types of particles, prints errors if particles are already existing
    or if there are too many particles*/
  void static AddParticleType(const char *name, double mass, int charge,
                              double width);
  void static AddParticleType(const char *name, double mass, int charge);

  // particle decay in dau1 and dau2
  int Decay2body(Particle &dau1, Particle &dau2) const;

  double Momentum2() const;
  double Momentum() const;

  // particle energy
  double Energy() const;

  // invariant mass computed with another particle p
  double InvMass(Particle &p) const;

  /*sets manually the particle index, you can use this member function with name
    or with particle code*/
  void SetIndex(int index);
  void SetIndex(const char *name);
  void SetP(double px, double py, double pz);

  // prints the whole particle types vector
  static void PrintParticleTypes();

  void PrintParticle();

 private:
  double p_x_;
  double p_y_;
  double p_z_;

  // index to evaluate the type of particle
  int fIndex;

  // max possible number of particle
  static const int fMaxNumParticleType = 10;

  // vector with existing types of particles
  static std::vector<ParticleType *> fParticleType;

  // find types of particles by name, returns position of particles
  // returns -1 if there isn't any particle with this name
  int static FindParticle(const char *name);

  // useful for AddParticleType
  void Boost(double bx, double by, double bz);
};

#endif