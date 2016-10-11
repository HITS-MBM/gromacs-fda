/*
 * Functions for per-atom interactions.
 *
 * Copyright Bogdan Costescu 2011-2013
 */

#include <math.h>

#include "gromacs/fda/pf_array.h"
#include "gromacs/fda/pf_per_atom.h"
#include "gromacs/fda/types/pf_array_summed.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

void pf_per_atom_real_set(t_pf_per_atom_real *per_atom_real, real val) {
  int i;

  for (i = 0; i < per_atom_real->len; i++)
    per_atom_real->force[i] = val;
}

void pf_per_atom_real_init(t_pf_per_atom_real **per_atom_real, int len, real val) {
  t_pf_per_atom_real *p;
  int i;

  snew(p, 1);
  p->len = len;
  snew(p->force, len);
  pf_per_atom_real_set(p, val);
  *per_atom_real = p;
}

void pf_per_atom_real_int_set(t_pf_per_atom_real_int *per_atom_real_int, real val_real, int val_int) {
  int i;

  for (i = 0; i < per_atom_real_int->len; i++) {
    per_atom_real_int->force[i] = val_int;
    per_atom_real_int->interactions[i] = val_real;
  }
}

void pf_per_atom_real_int_init(t_pf_per_atom_real_int **per_atom_real_int, int len, real val_real, int val_int) {
  t_pf_per_atom_real_int *p;
  int i;

  snew(p, 1);
  p->len = len;
  snew(p->force, len);
  snew(p->interactions, len);
  pf_per_atom_real_int_set(p, val_real, val_int);
  *per_atom_real_int = p;
}

void pf_per_atom_real_write_frame(FILE *f, real *force, int len, gmx_bool no_end_zeros) {
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

void pf_per_atom_sum(t_pf_per_atom_real *per_atom_sum, t_pf_atom_summed *atoms, int atoms_len, const rvec *x, int Vector2Scalar) {
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

/* findmax parameter: TRUE = find max; FALSE = find min */
static gmx_inline void pf_per_atom_minmax_op(real *force, real val, gmx_bool findmax) {
  /* val is a norm of a vector so it's always non-negative */
  //fprintf(stderr, "force=%f, val=%f\n", *force, val);
  if (findmax) {
    /* looking for maximum */
    if (val > (*force))
      *force = val;
  } else {
    /*looking for minimum */
    if (((*force) == 0.0) || (val < (*force)))
      *force = val;
  }
}

/* findmax parameter: TRUE = find max; FALSE = find min */
void pf_per_atom_minmax(t_pf_per_atom_real *per_atom_real, t_pf_atom_summed *atoms, int atoms_len, gmx_bool findmax, const rvec *x, int Vector2Scalar) {
  int i, ii, j, jj;
  t_pf_interaction_array_summed ias;
  t_pf_interaction_summed is;
  real scalar_force = 0.0;

  //fprintf(stderr, "pf_per_atom_minmax, findmax=%d, atoms_len=%d\n", findmax ? 1 : 0, atoms_len);
  pf_per_atom_real_set(per_atom_real, 0.0);
  for (i = 0; i < atoms_len; i++) {
    ii = atoms[i].nr;
    ias = atoms[i].interactions;
    //fprintf(stderr, "ii=%d,  nr. interactions=%d\n", ii, ias.len);
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
          break;
      }
      //fprintf(stderr, "force=%32f\n", scalar_force);
      pf_per_atom_minmax_op(&(per_atom_real->force[ii]), scalar_force, findmax);
      pf_per_atom_minmax_op(&(per_atom_real->force[jj]), scalar_force, findmax);
    }
  }
}

static gmx_inline real pf_per_atom_average_compute(real r, int i) {
  switch (i) {
    case 0:
      /* there were no interactions, so return 0 */
      return 0.0;
      break;
    case 1:
      /* avoid division by 1 */
      return r;
      break;
    default:
      return (r / i);
      break;
  }
}

/* computes averages in place and returns a t_pf_per_atom_real pointing to it;
 * this way, only one large memory allocation is made at beginning for the pf->global->per_atom_real_int
 */
void pf_per_atom_average(t_pf_per_atom_real_int *per_atom_average, t_pf_atom_summed *atoms, int atoms_len, const rvec *x, int Vector2Scalar) {
  int i, ii, j, jj;
  t_pf_interaction_array_summed ias;
  t_pf_interaction_summed is;
  real scalar_force = 0.0;

  /* a counter for interactions is needed because not all interactions are accounted for a certain atom;
   * when atoms i and j interact, the interaction is only stored once (for whichever i and j is smaller), but should be accounted for both i and j
   */
  pf_per_atom_real_int_set(per_atom_average, 0.0, 0);
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
          break;
      }
      per_atom_average->force[ii] += scalar_force;
      per_atom_average->force[jj] += scalar_force;
      (per_atom_average->interactions[ii])++;
      (per_atom_average->interactions[jj])++;
    }
  }
  /* summing is done, now compute the average */
  for (i = 0; i < per_atom_average->len; i++) {
    per_atom_average->force[i] = pf_per_atom_average_compute(per_atom_average->force[i], per_atom_average->interactions[i]);
  }
}

/**
 * The stress is the negative atom_vir value.
 */
void pf_write_atom_virial_sum(FILE *f, tensor *atom_vir, int natoms)
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

/**
 * For von Mises no negative values are needed, since all items are squared.
 */
void pf_write_atom_virial_sum_von_mises(FILE *f, tensor *atom_vir, int natoms)
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
