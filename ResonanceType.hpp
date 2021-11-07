#ifndef RESONANCE_TYPE
#define RESONANCE_TYPE

#include "ParticleType.hpp"

class ResonanceType : public ParticleType {
 public:
  // costruttore parametrico
  ResonanceType(const char* name, double mass, int charge, double width);

  double GetWidth() const;

  // metodo print come nella classe madre
  void Print() const;

 private:
  double const fWidth;
};

#endif