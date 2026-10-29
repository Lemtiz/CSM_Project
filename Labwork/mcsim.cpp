#include "./sim.hpp"
int main(int argc, char **argv)
{
  mcsim<double, particleLJ<double>> mc;
  mc.init_rng(0); //  initialize the random number generator
  mc.prepare_initial_conf(); // create initial configuration
  mc.run(); // launch MC simulation
  return 0;
}
