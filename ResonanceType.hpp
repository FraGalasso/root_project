#ifndef RESONANCE_TYPE
#define RESONANCE_TYPE

#include "ParticleType.hpp"

class ResonanceType : public ParticleType {
 public:
  // parametric constructor
  ResonanceType(const char* name, double mass, int charge, double width);

  double GetWidth() const;

  // similar to the one in the base class
  void Print() const;

 private:
  double const fWidth;
};

#endif