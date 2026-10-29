#include "./sim.hpp"
int main(int argc, char **argv)
{
  mdsim<double,particleLJ<double>> md;
  md.init_rng(0); // inizializzo il generatore di numeri casuali
  md.prepare_initial_conf(); // creo la configurazione iniziale
  md.run(); // lancio la simulazione MC
  return 0;
}
