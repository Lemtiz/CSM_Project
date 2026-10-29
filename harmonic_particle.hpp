#ifndef HARMONIC_PARTICLE_HPP
#define HARMONIC_PARTICLE_HPP

#include "./Labwork/pvector.hpp"

struct harmonic_particle
{
  pvector<double, 3> r;
  pvector<double, 3> v;
  pvector<double, 3> f;
  pvector<double, 3> f_lj;
  pvector<double, 3> f_harmonic;
};

#endif
