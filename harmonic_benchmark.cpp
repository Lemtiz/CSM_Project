#include "./harmonic_benchmark.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

benchmark_result benchmark_k(
    double k,
    double total_time,
    double dt_slow,
    int seed)
{
  harmonic_params parameters;
  parameters.k = k;
  parameters.dt_slow = dt_slow;
  parameters.seed = seed;
  harmonic_mdsim simulation(parameters);

  benchmark_result result;
  result.k = k;
  result.dt_slow = dt_slow;
  result.substeps = simulation.napa_substeps();
  result.dt_vv = 0.1 / std::sqrt(k / parameters.mass);
  result.vv_seconds = simulation.run_velocity_verlet(total_time);
  result.napa_seconds = simulation.run_r_napa(total_time);
  result.speedup = result.vv_seconds / result.napa_seconds;
  return result;
}

void write_benchmark_csv(
    const std::string& path,
    const std::vector<benchmark_result>& results)
{
  std::ofstream output(path);
  if (!output)
    throw std::runtime_error("Cannot open benchmark output: " + path);
  output << "k,dt_vv,dt_slow,substeps,vv_seconds,napa_seconds,speedup\n";
  output << std::setprecision(12);
  for (const benchmark_result& result : results)
    output << result.k << ',' << result.dt_vv << ','
           << result.dt_slow << ',' << result.substeps << ','
           << result.vv_seconds << ',' << result.napa_seconds << ','
           << result.speedup << '\n';
}
