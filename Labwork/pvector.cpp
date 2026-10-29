#include "./pvector.hpp"
int main(void)
{
  pvector<double,3> P, P1 = {1.3,2.7,3.4};
  pvector<float, 4> V;
  double L=4.5;
  P(2)=3.0;
  std::cout << "P(2)=" << P(2) << "\n"; 
 
  const pvector<double,3> PC = {0.2,0.3,0.4};
  std::cout << "PC=(" << PC(0) << "," << PC(1) << "," << PC(2) << ")\n";

#if 1
  P = L*rint(P1/L);
  P.show("P=");
  // Cases for which it is necessary to have ntype operator()(int idx) const { ... } overloaded method
 #endif
  return 0;
}
