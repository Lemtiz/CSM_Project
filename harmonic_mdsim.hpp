#ifndef HARMONIC_MDSIM_HPP
#define HARMONIC_MDSIM_HPP

#include <algorithm>
#include <chrono>
#include <cmath>

#include "./harmonic_sim.hpp"
#include "./Labwork/randnumgen.hpp"

class harmonic_mdsim : public harmonic_sim
{
  static constexpr double stability_factor = 0.1;
  static constexpr double pi = 3.14159265358979323846;

  void initialize_velocities()
  {
    rng.seed(pars.seed);
    const double standard_deviation =
        std::sqrt(pars.temperature / pars.mass);
    pvector<double, 3> center_of_mass = {0.0, 0.0, 0.0};

    for (harmonic_particle& part : parts)
    {
      double u1 = 0.0;
      double u2 = 0.0;
      do
      {
        u1 = rng.ranf();
        u2 = rng.ranf();
      } while (u1 == 0.0);
      const double radius = std::sqrt(-2.0 * std::log(u1));
      part.v = {
          standard_deviation * radius * std::cos(2.0 * pi * u2),
          standard_deviation * radius * std::sin(2.0 * pi * u2),
          standard_deviation * (2.0 * rng.ranf() - 1.0)};
      center_of_mass += part.v;
    }

    center_of_mass /= static_cast<double>(parts.size());
    for (harmonic_particle& part : parts)
      part.v -= center_of_mass;
  }

  void reset()
  {
    prepare_initial_conf();
    initialize_velocities();
    calculate_forces();
  }

  void vv_step(double dt)
  {
    calculate_forces();
    for (std::size_t i = 0; i < parts.size(); ++i)
    {
      parts[i].r += parts[i].v * dt +
                    (0.5 * dt * dt / pars.mass) * parts[i].f;
      pbc(static_cast<int>(i));
    }
    calculate_forces();
    for (harmonic_particle& part : parts)
      part.v += (0.5 * dt / pars.mass) * part.f;
  }

  void slow_vv_step(double dt)
  {
    calc_lj_forces();
    for (std::size_t i = 0; i < parts.size(); ++i)
    {
      parts[i].r += parts[i].v * dt +
                    (0.5 * dt * dt / pars.mass) * parts[i].f_lj;
      pbc(static_cast<int>(i));
    }
    calc_lj_forces();
    for (harmonic_particle& part : parts)
      part.v += (0.5 * dt / pars.mass) * part.f_lj;
  }

  void fast_step(double dt)
  {
    parts[0].f_harmonic = (-pars.k) * parts[0].r;
    const pvector<double, 3> initial =
        parts[0].f_harmonic / pars.mass;
    parts[0].r += parts[0].v * dt + 0.5 * dt * dt * initial;
    pbc(0);
    parts[0].f_harmonic = (-pars.k) * parts[0].r;
    const pvector<double, 3> final =
        parts[0].f_harmonic / pars.mass;
    parts[0].v += 0.5 * dt * (initial + final);
  }

public:
  explicit harmonic_mdsim(const harmonic_params& parameters)
      : harmonic_sim(parameters)
  {
    reset();
  }

  double run_velocity_verlet(double total_time)
  {
    reset();
    const double dt = stability_factor /
                      std::sqrt(pars.k / pars.mass);
    const long long steps = std::max(
        1LL, static_cast<long long>(std::ceil(total_time / dt)));
    const double actual_dt = total_time / static_cast<double>(steps);
    const auto start = std::chrono::steady_clock::now();
    for (long long step = 0; step < steps; ++step)
      vv_step(actual_dt);
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(end - start).count();
  }

  int napa_substeps() const
  {
    const double fast_dt = stability_factor /
                           std::sqrt(pars.k / pars.mass);
    return std::max(
        1, static_cast<int>(std::ceil(pars.dt_slow / fast_dt)));
  }

  double run_r_napa(double total_time)
  {
    reset();
    const long long steps = std::max(
        1LL, static_cast<long long>(
                 std::ceil(total_time / pars.dt_slow)));
    const double slow_dt = total_time / static_cast<double>(steps);
    const int substeps = napa_substeps();
    const double fast_dt = slow_dt / static_cast<double>(substeps);
    const auto start = std::chrono::steady_clock::now();
    for (long long step = 0; step < steps; ++step)
    {
      slow_vv_step(0.5 * slow_dt);
      for (int substep = 0; substep < substeps; ++substep)
        fast_step(fast_dt);
      slow_vv_step(0.5 * slow_dt);
    }
    const auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(end - start).count();
  }
};

#endif
