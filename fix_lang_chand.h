#ifdef FIX_CLASS
// Registrar el estilo de fix y vincular el nombre del comando
FixStyle(langevin/chand, FixLangChand)

#else

#ifndef LMP_FIX_LANG_CHAND_H
#define LMP_FIX_LANG_CHAND_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLangChand : public Fix {
 public:
  FixLangChand(class LAMMPS *, int, char **);
  ~FixLangChand() override;
  int setmask() override;
  void init() override;
  void setup(int vflag) override;
  void initial_integrate(int vflag) override;
  void final_integrate() override; // <- Sin el argumento vflag
  void reset_dt() override;
  void grow_arrays(int nmax) override;
  double memory_usage() override;

 private:
  double t_start, t_stop, damp, t_target;
  int seed;
  double c0, c1, c2;

  double **f_old;
  double **Vn_store;

  class RanMars *random;

  void calculate_coefficients();
  void generate_correlated_noise(double m, double &Xn, double &Vn);
};

}

#endif
#endif
