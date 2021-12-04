#ifndef PARTICLE_TYPE
#define PARTICLE_TYPE

#include <iostream>

class ParticleType {
 public:
  // parametric constructor
  ParticleType(const char* name, double mass, int charge);

  const char* GetName() const;
  double GetMass() const;
  int GetCharge() const;

  /*I need it to be virtual to use it both on a particleType and on a
    resonanceType*/
  virtual double GetWidth() const;

  // prints attributes
  virtual void Print() const;

 private:
  const char* fName;
  double const fMass;
  int const fCharge;
};

#endif