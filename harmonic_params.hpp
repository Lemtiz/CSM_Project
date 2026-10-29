#ifndef HARMONIC_PARAMS_HPP
#define HARMONIC_PARAMS_HPP

#include <cmath>
#include <stdexcept>

struct harmonic_params
{
  int nx = 5;
  int ny = 5;
  int nz = 5;
  double temperature = 1.0;
  double density = 0.5;
  double sigma = 1.0;
  double epsilon = 1.0;
  double mass = 1.0;
  double rc = 2.5;
  double k = 100.0;
  double dt_slow = 0.005;
  int seed = 0;

  int particle_count() const
  {
    return nx * ny * nz;
  }

  double box_side() const
  {
    return std::cbrt(static_cast<double>(particle_count()) / density);
  }

  void validate() const
  {
    if (nx <= 0 || ny <= 0 || nz <= 0 || density <= 0.0 ||
        temperature < 0.0 || sigma <= 0.0 || epsilon <= 0.0 ||
        mass <= 0.0 || rc <= 0.0 || k <= 0.0 || dt_slow <= 0.0)
      throw std::invalid_argument("Invalid simulation parameters");
    if (rc >= 0.5 * box_side())
      throw std::invalid_argument("rc must be smaller than half the box");
  }
};

#endif
