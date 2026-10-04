// Test harness for transform_core.h (compiled with g++ outside R)
#include <cstdio>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include "../../../src/transform_core.h"
int main(int argc, char** argv) {
  // input: params.txt, grid.txt, inputs.txt (z1 w u_type e_pps u_ttr per line)
  std::ifstream fp(argv[1]);
  coronco::ArmParams p;
  fp >> p.lam_p >> p.theta >> p.c_resp >> p.z_tau >> p.tau >> p.os_model >> p.pi_death >> p.gam0 >> p.gam1
     >> p.lam_o >> p.c_dec >> p.grid_h >> p.ttr_mode >> p.ttr_rate;
  std::vector<double> g;
  std::ifstream fg(argv[2]);
  double v;
  while (fg >> v) g.push_back(v);
  p.grid_g = g.data(); p.grid_n = (int) g.size();
  std::ifstream fi(argv[3]);
  std::FILE* fo = std::fopen(argv[4], "w");
  double z1, w, ut, ep, uttr;
  while (fi >> z1 >> w >> ut >> ep >> uttr) {
    double pfs, os, ttr; int prog, resp;
    coronco::transform_one(z1, w, ut, ep, uttr, p, pfs, os, prog, resp, ttr);
    std::fprintf(fo, "%.17g %.17g %d %d %.17g\n", pfs, os, prog, resp, ttr);
  }
  std::fclose(fo);
  // accrual check
  double at[3] = {0.0, 6.0, 24.0};
  double ac[3] = {0.0, 0.25, 1.0};
  std::FILE* fa = std::fopen(argv[5], "w");
  for (int i = 0; i <= 20; ++i) {
    double u = i / 20.0;
    std::fprintf(fa, "%.17g %.17g\n", u, coronco::accrual_inverse(u, at, ac, 2));
  }
  std::fclose(fa);
  return 0;
}
