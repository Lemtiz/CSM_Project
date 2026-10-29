#ifndef HARMONIC_BENCHMARK_HPP
#define HARMONIC_BENCHMARK_HPP

#include <string>
#include <vector>

#include "./harmonic_mdsim.hpp"

struct benchmark_result
{
  double k = 0.0;
  double dt_vv = 0.0;
  double dt_slow = 0.0;
  int substeps = 0;
  double vv_seconds = 0.0;
  double napa_seconds = 0.0;
  double speedup = 0.0;
};

benchmark_result benchmark_k(
    double k,
    double total_time,
    double dt_slow,
    int seed);

void write_benchmark_csv(
    const std::string& path,
    const std::vector<benchmark_result>& results);

#endif
