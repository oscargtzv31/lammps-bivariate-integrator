#include "fix_lang_chand.h"
#include <cmath>
#include <algorithm>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "random_mars.h"
#include "memory.h"
#include "utils.h"

using namespace LAMMPS_NS;

FixLangChand::FixLangChand(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), f_old(nullptr), Vn_store(nullptr), random(nullptr)
{
  if (narg < 7) 
    error->all(FLERR, "Sintaxis: fix ID group langevin/chand Tstart Tstop damp seed");

  time_integrate = 1;

  t_start = utils::numeric(FLERR, arg[3], false, lmp);
  t_stop  = utils::numeric(FLERR, arg[4], false, lmp);
  damp    = utils::numeric(FLERR, arg[5], false, lmp);
  seed    = utils::inumeric(FLERR, arg[6], false, lmp);

  if (damp <= 0.0) error->all(FLERR, "Fix langevin/chand: damp should be > 0.0");
  if (seed <= 0)   error->all(FLERR, "Fix langevin/chand: seed should be > 0");

  // Seed processor independent 
  random = new RanMars(lmp, seed + comm->me);

  grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);
}

FixLangChand::~FixLangChand()
{
  delete random;
  memory->destroy(f_old);
  memory->destroy(Vn_store);
  if (atom) atom->delete_callback(id, Atom::GROW);
}

int FixLangChand::setmask()
{
  int mask = 0;
  mask |= FixConst::INITIAL_INTEGRATE;
  mask |= FixConst::FINAL_INTEGRATE;
  return mask;
}

void FixLangChand::init()
{
  calculate_coefficients();
}

void FixLangChand::setup(int vflag)
{
  calculate_coefficients();
}

void FixLangChand::grow_arrays(int nmax)
{
  memory->grow(f_old, nmax, 3, "fix_lang_chand:f_old");
  memory->grow(Vn_store, nmax, 3, "fix_lang_chand:Vn_store");
}

double FixLangChand::memory_usage()
{
  double bytes = 0.0;
  bytes += (double) atom->nmax * 3 * sizeof(double);
  bytes += (double) atom->nmax * 3 * sizeof(double);
  return bytes;
}

void FixLangChand::calculate_coefficients()
{
  double dt = update->dt;
  double gamma = 1.0 / damp; // xi = 1/damp
  double x = gamma * dt;     // x = xi * dt (adimensional)

  // Rampa de temperatura entre t_start y t_stop
  double delta = update->ntimestep - update->beginstep;
  if (update->endstep == update->beginstep) {
    t_target = t_start;
  } else {
    double delta_total = update->endstep - update->beginstep;
    t_target = t_start + (t_stop - t_start) * (delta / delta_total);
  }

  // Preventing numerical unstabillity with Taylor series 
  if (x < 1.0e-4) {
    c0 = 1.0 - x + 0.5 * x * x - (1.0 / 6.0) * x * x * x + (1.0 / 24.0) * x * x * x * x;
    c1 = 1.0 - 0.5 * x + (1.0 / 6.0) * x * x - (1.0 / 24.0) * x * x * x;
    c2 = 0.5 - (1.0 / 6.0) * x + (1.0 / 24.0) * x * x - (1.0 / 120.0) * x * x * x;
  } else {
    c0 = std::exp(-x);
    c1 = (1.0 - c0) / x;
    c2 = (1.0 - c1) / x;
  }
}

void FixLangChand::generate_correlated_noise(double m, double &Xn, double &Vn)
{
  double dt = update->dt;
  double gamma = 1.0 / damp;
  double x = gamma * dt;
  
  double boltz = force->boltz;
  double ftm2v = force->ftm2v;
  
  // Factor de escala térmico-cinético exacto por unidad de masa k_B * T / m
  double kBT_m = (boltz * t_target / m) * ftm2v;

  double var_X, var_V, cov_XV;

  // Evaluación numéricamente estable de los segundos momentos (Eq. 25 y 26)
  if (x < 1.0e-4) {
    var_X  = (2.0 / 3.0) * kBT_m * x * dt * dt * (1.0 - 0.75 * x + 0.35 * x * x);
    var_V  = 2.0 * kBT_m * x * (1.0 - x + (2.0 / 3.0) * x * x);
    cov_XV = kBT_m * x * dt * (1.0 - x + (7.0 / 12.0) * x * x);
  } else {
    var_X  = (kBT_m / (gamma * gamma)) * (2.0 * x + 4.0 * std::exp(-x) - 3.0 - std::exp(-2.0 * x));
    var_V  = kBT_m * (1.0 - std::exp(-2.0 * x));
    cov_XV = (kBT_m / gamma) * std::pow(1.0 - std::exp(-x), 2);
  }

  double sigma_x = std::sqrt(std::max(0.0, var_X));
  double z1 = random->gaussian();
  double z2 = random->gaussian();

  Xn = sigma_x * z1;
  
  double var_v_cond = std::max(0.0, var_V - (cov_XV * cov_XV) / std::max(1.0e-18, var_X));
  Vn = (cov_XV / std::max(1.0e-18, sigma_x)) * z1 + std::sqrt(var_v_cond) * z2;
}

void FixLangChand::initial_integrate(int vflag)
{
  if (t_start != t_stop) calculate_coefficients();

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double ftm2v = force->ftm2v;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double m = rmass ? rmass[i] : mass[type[i]];
      double dtm = ftm2v / m;

      f_old[i][0] = f[i][0];
      f_old[i][1] = f[i][1];
      f_old[i][2] = f[i][2];

      double Xn_x, Vn_x, Xn_y, Vn_y, Xn_z, Vn_z;

      generate_correlated_noise(m, Xn_x, Vn_x);
      generate_correlated_noise(m, Xn_y, Vn_y);
      generate_correlated_noise(m, Xn_z, Vn_z);

      Vn_store[i][0] = Vn_x;
      Vn_store[i][1] = Vn_y;
      Vn_store[i][2] = Vn_z;

      // Actualización de posiciones
      x[i][0] += c1 * v[i][0] * dt + c2 * f[i][0] * dtm * dt * dt + Xn_x;
      x[i][1] += c1 * v[i][1] * dt + c2 * f[i][1] * dtm * dt * dt + Xn_y;
      x[i][2] += c1 * v[i][2] * dt + c2 * f[i][2] * dtm * dt * dt + Xn_z;
    }
  }
}

void FixLangChand::final_integrate()
{
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = update->dt;
  double ftm2v = force->ftm2v;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double m = rmass ? rmass[i] : mass[type[i]];
      double dtm = ftm2v / m;

      // Actualización de velocidades
      v[i][0] = c0 * v[i][0] + (c1 - c2) * dtm * dt * f_old[i][0] + c2 * dtm * dt * f[i][0] + Vn_store[i][0];
      v[i][1] = c0 * v[i][1] + (c1 - c2) * dtm * dt * f_old[i][1] + c2 * dtm * dt * f[i][1] + Vn_store[i][1];
      v[i][2] = c0 * v[i][2] + (c1 - c2) * dtm * dt * f_old[i][2] + c2 * dtm * dt * f[i][2] + Vn_store[i][2];
    }
  }
}

void FixLangChand::reset_dt()
{
  calculate_coefficients();
}
