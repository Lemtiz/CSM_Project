#include "./harmonic_benchmark.hpp"

#include <cstdlib>

int main(int argc, char** argv)
{
  const double total_time = argc > 1 ? std::atof(argv[1]) : 10.0;
  const double dt_slow = argc > 2 ? std::atof(argv[2]) : 0.005;
  const std::string output =
      argc > 3 ? argv[3] : "benchmark.csv";
  const std::vector<double> k_values = {
      25.0, 50.0, 100.0, 200.0, 400.0, 800.0,
      1600.0, 3200.0, 6400.0, 12800.0, 25600.0, 51200.0,
      102400.0};

  std::vector<benchmark_result> results;
  results.reserve(k_values.size());
  for (double k : k_values)
    results.push_back(benchmark_k(k, total_time, dt_slow, 0));

  write_benchmark_csv(output, results);
  std::cout << "k dt_vv dt_slow substeps VV(s) r-NAPA(s) speedup\n";
  for (const benchmark_result& result : results)
  {
    std::cout << result.k << ' '
              << result.dt_vv << ' '
              << result.dt_slow << ' '
              << result.substeps << ' '
              << result.vv_seconds << ' '
              << result.napa_seconds << ' '
              << result.speedup << '\n';
  }
  return 0;
}
