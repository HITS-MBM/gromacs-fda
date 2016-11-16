/*
 * Functions for per-atom interactions.
 *
 * Copyright Bogdan Costescu 2011-2013
 */

#include <cmath>
#include "fda.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

using namespace fda;

void FDA::per_atom_real_write_frame(FILE *f, std::vector<real> const& force_per_node, gmx_bool no_end_zeros)
{
  int i, j;
  gmx_bool first_on_line = TRUE;

  j = len;
  if (no_end_zeros) {
    /* detect the last non-zero item */
    for (j = len; j > 0; j--)
      if (force[j - 1] != 0.0)
        break;
  }
  /* j holds the index of first zero item or the length of force */
  for (i = 0; i < j; i++) {
    if (first_on_line) {
      fprintf(f, "%e", force[i]);
      first_on_line = FALSE;
    } else {
      fprintf(f, " %e", force[i]);
    }
  }
  fprintf(f, "\n");
}

static inline real pf_vector2unsignedscalar(const rvec v, const int i, const int j, const rvec *x) {
  rvec r;
  real p;

  rvec_sub(x[j], x[i], r);
  p = cos_angle(r, v) * norm(v);
  if (p < 0)
    return -p;
  return p;
}

void FDA::per_atom_sum(std::vector<real>& force_per_node, t_pf_atom_summed *atoms, int atoms_len, const rvec *x, int Vector2Scalar)
{
  int i, ii, j, jj;
  t_pf_interaction_array_summed ias;
  t_pf_interaction_summed is;
  real scalar_force = 0.0;

  pf_per_atom_real_set(per_atom_sum, 0.0);
  for (i = 0; i < atoms_len; i++) {
    ii = atoms[i].nr;
    ias = atoms[i].interactions;
    for (j = 0; j < ias.len; j++) {
      is = ias.array[j];
      jj = is.jjnr;
      switch (Vector2Scalar) {
        case PF_VECTOR2SCALAR_NORM:
          scalar_force = norm(is.force);
          break;
        case PF_VECTOR2SCALAR_PROJECTION:
          scalar_force = pf_vector2unsignedscalar(is.force, ii, jj, x);
          break;
        default:
      	  gmx_fatal(FARGS, "Unknown option for Vector2Scalar.\n");
          break;
      }
      per_atom_sum->force[ii] += scalar_force;
      per_atom_sum->force[jj] += scalar_force;
    }
  }
}

void FDA::write_atom_virial_sum(FILE *f, tensor *atom_vir, int natoms)
{
  int i;
  gmx_bool first_on_line = TRUE;

  for (i = 0; i < natoms; i++) {
    if (first_on_line) {
      fprintf(f, "%e %e %e %e %e %e", -atom_vir[i][XX][XX], -atom_vir[i][YY][YY], -atom_vir[i][ZZ][ZZ], -atom_vir[i][XX][YY], -atom_vir[i][XX][ZZ], -atom_vir[i][YY][ZZ]);
      first_on_line = FALSE;
    } else {
      fprintf(f, " %e %e %e %e %e %e", -atom_vir[i][XX][XX], -atom_vir[i][YY][YY], -atom_vir[i][ZZ][ZZ], -atom_vir[i][XX][YY], -atom_vir[i][XX][ZZ], -atom_vir[i][YY][ZZ]);
    }
  }
  fprintf(f, "\n");
}

static gmx_inline real tensor_to_vonmises(tensor t)
{
  real txy, tyz, tzx;

  txy = t[XX][XX]-t[YY][YY];
  tyz = t[YY][YY]-t[ZZ][ZZ];
  tzx = t[ZZ][ZZ]-t[XX][XX];
  return sqrt(0.5 * (txy*txy + tyz*tyz + tzx*tzx +
    6 * (t[XX][YY]*t[XX][YY] + t[XX][ZZ]*t[XX][ZZ] + t[YY][ZZ]*t[YY][ZZ])));
}

void FDA::write_atom_virial_sum_von_mises(FILE *f, tensor *atom_vir, int natoms)
{
  int i;
  gmx_bool first_on_line = TRUE;

  for (i = 0; i < natoms; i++) {
    if (first_on_line) {
      fprintf(f, "%e", tensor_to_vonmises(atom_vir[i]));
      first_on_line = FALSE;
    } else {
      fprintf(f, " %e", tensor_to_vonmises(atom_vir[i]));
    }
  }
  fprintf(f, "\n");
}
