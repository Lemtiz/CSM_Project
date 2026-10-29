#ifndef HARMONIC_SIM_HPP
#define HARMONIC_SIM_HPP

#include <cmath>
#include <vector>

#include "./harmonic_params.hpp"
#include "./harmonic_particle.hpp"

class harmonic_sim
{
protected:
  harmonic_params pars;
  pvector<double, 3> L;
  std::vector<harmonic_particle> parts;

  pvector<double, 3> minimum_image(
      const pvector<double, 3>& displacement) const
  {
    return displacement - mulcw(L, rint(divcw(displacement, L)));
  }

  void pbc(int i)
  {
    parts[i].r -= mulcw(L, rint(divcw(parts[i].r, L)));
  }

  void prepare_initial_conf()
  {
    const double side = pars.box_side();
    L = {side, side, side};
    parts.assign(static_cast<std::size_t>(pars.particle_count()),
                 harmonic_particle{});

    int index = 0;
    for (int ix = 0; ix < pars.nx; ++ix)
      for (int iy = 0; iy < pars.ny; ++iy)
        for (int iz = 0; iz < pars.nz; ++iz)
        {
          parts[index].r = {
              (ix + 0.5) * side / pars.nx - 0.5 * side,
              (iy + 0.5) * side / pars.ny - 0.5 * side,
              (iz + 0.5) * side / pars.nz - 0.5 * side};
          ++index;
        }
    parts[0].r = {0.25, 0.0, 0.0};
  }

  void calc_lj_forces()
  {
    for (harmonic_particle& part : parts)
      part.f_lj = {0.0, 0.0, 0.0};

    for (std::size_t i = 0; i < parts.size(); ++i)
      for (std::size_t j = i + 1; j < parts.size(); ++j)
      {
        const pvector<double, 3> dr = minimum_image(parts[i].r - parts[j].r);
        const double r2 = dr * dr;
        if (r2 == 0.0)
          continue;
        const double r = std::sqrt(r2);
        if (r >= pars.rc)
          continue;

        const double sr2 = pars.sigma * pars.sigma / r2;
        const double sr6 = sr2 * sr2 * sr2;
        const double sr12 = sr6 * sr6;
        const double factor =
            24.0 * pars.epsilon * (2.0 * sr12 - sr6) / r2;
        const pvector<double, 3> fij = factor * dr;
        parts[i].f_lj += fij;
        parts[j].f_lj -= fij;
      }
  }

  double lj_energy(double r) const
  {
    if (r >= pars.rc)
      return 0.0;
    const double sr = pars.sigma / r;
    const double sr6 = std::pow(sr, 6.0);
    const double src = pars.sigma / pars.rc;
    const double src6 = std::pow(src, 6.0);
    return 4.0 * pars.epsilon *
           (sr6 * sr6 - sr6 - src6 * src6 + src6);
  }

  double potential_energy() const
  {
    double potential = 0.0;
    for (std::size_t i = 0; i < parts.size(); ++i)
    {
      if (i == 0)
        potential += 0.5 * pars.k * (parts[i].r * parts[i].r);
      for (std::size_t j = i + 1; j < parts.size(); ++j)
        potential += lj_energy(
            minimum_image(parts[i].r - parts[j].r).norm());
    }
    return potential;
  }

  double kinetic_energy() const
  {
    double kinetic = 0.0;
    for (const harmonic_particle& part : parts)
      kinetic += 0.5 * pars.mass * (part.v * part.v);
    return kinetic;
  }

  void calculate_forces()
  {
    calc_lj_forces();
    parts[0].f_harmonic = (-pars.k) * parts[0].r;
    for (std::size_t i = 1; i < parts.size(); ++i)
      parts[i].f_harmonic = {0.0, 0.0, 0.0};
    for (harmonic_particle& part : parts)
      part.f = part.f_lj + part.f_harmonic;
  }

public:
  explicit harmonic_sim(const harmonic_params& parameters)
      : pars(parameters)
  {
    pars.validate();
    prepare_initial_conf();
  }

  const harmonic_params& parameters() const { return pars; }
  const pvector<double, 3>& box() const { return L; }
  double total_energy() const { return kinetic_energy() + potential_energy(); }
};

#endif
