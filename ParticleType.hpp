#ifndef PARTICLE_TYPE
#define PARTICLE_TYPE

#include <iostream>

class ParticleType {
 public:
  // costruttore parametrico
  ParticleType(const char* name, double mass, int charge);

  const char* GetName() const;
  double GetMass() const;
  int GetCharge() const;

  /*mi serve impementarlo virtual per poterlo chiamare su un generico
    particletype*/
  virtual double GetWidth() const;

  // metodo che stampa a schermo gli attributi
  virtual void Print() const;

 private:
  const char* fName;
  double const fMass;
  int const fCharge;
};

#endif