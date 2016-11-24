/*
 * I/O
 * 
 * Copyright Bogdan Costescu 2011-2013
 */

#include "fda.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

using namespace fda;

/* below a few definitions from the original PF implementation needed for compatibility mode */

/* original PF implementation [Stacklies] had the following values:
 * enum {iNone, iBond, iAngle, iDihedral, iPolar, iLJ, iCoul, iAll};
 * each pair interacts only through one type due to the approximations for angles and dihedrals
 * and the j<i interaction holds Coulomb (bottom left triangle) while i>j holds LJ (top right triangle);
 * to make the mapping possible, the checking is done in a particular order - this means that for
 * a pair which interacts through both bond and angle, the result will be marked as only bond; also
 * for a pair which interacts through both LJ and Coulomb, the result will be makes as only LJ
 */
enum {
  pf_compat_iNone,
  pf_compat_iBond,
  pf_compat_iAngle,
  pf_compat_iDihedral,
  pf_compat_iPolar,
  pf_compat_iLJ,
  pf_compat_iCoul,
  pf_compat_iAll
};

/* Mark the end of an entry in binary output files */
#define PF_COMPAT_NEW_ENTRY -280480

/* original PF implementation stored interaction as a char */
static inline char pf_compatibility_map_interaction(int type) {
  if (type & PF_INTER_BOND)
    return pf_compat_iBond;
  if (type & PF_INTER_ANGLE)
    return pf_compat_iAngle;
  if (type & PF_INTER_DIHEDRAL)
    return pf_compat_iDihedral;
  if (type & PF_INTER_POLAR)
    return pf_compat_iPolar;
  if (type & PF_INTER_LJ)
    return pf_compat_iLJ;
  if (type & PF_INTER_COULOMB)
    return pf_compat_iCoul;
  if (type & PF_INTER_NB14)
    return pf_compat_iDihedral;
  gmx_fatal(FARGS, "Could not map interaction %d.", type);
  return pf_compat_iNone;	/* just to make the compiler happy */
}

gmx_int64_t pf_atoms_summed_to_fm(t_pf_atoms *atoms, int **fmatoms) {
  t_pf_atom_summed i_atom_summed;
  gmx_int64_t interactions_count;
  int *fm;
  int i;

  snew(fm, atoms->len);
  //fprintf(stderr, "pf_atoms_summed_to_fm: atoms->len=%d\n", atoms->len);
  /* find the total nr. of interactions, needed to size the interaction arrays */
  interactions_count = 0;
  for (i = 0; i < atoms->len; i++) {
    i_atom_summed = atoms->summed[i];
    fm[i] = i_atom_summed.nr;
    //fprintf(stderr, "fm[%d]=%d\n", i, fm[i]);
    interactions_count += i_atom_summed.interactions.len;
    //fprintf(stderr, "len=%d, i=%d, anr=%d, ilen=%d\n", atoms->len, i, i_atom_summed.nr, i_atom_summed.interactions.len);
  }
  if (interactions_count == 0) {
    sfree(fm);
    *fmatoms = NULL;
  } else {
    *fmatoms = fm;
  }
  return interactions_count;
}

gmx_int64_t pf_atoms_scalar_to_fm(t_pf_atoms *atoms, int **fmatoms) {
  t_pf_atom_scalar i_atom_scalar;
  gmx_int64_t interactions_count;
  int *fm;
  int i;

  snew(fm, atoms->len);
  //fprintf(stderr, "pf_atoms_scalar_to_fm: atoms->len=%d\n", atoms->len);
  /* find the total nr. of interactions, needed to size the interaction arrays */
  interactions_count = 0;
  for (i = 0; i < atoms->len; i++) {
    i_atom_scalar = atoms->scalar[i];
    fm[i] = i_atom_scalar.nr;
    //fprintf(stderr, "fm[%d]=%d\n", i, fm[i]);
    interactions_count += i_atom_scalar.interactions.len;
    //fprintf(stderr, "len=%d, i=%d, anr=%d, ilen=%d\n", atoms->len, i, i_atom_scalar.nr, i_atom_scalar.interactions.len);
  }
  if (interactions_count == 0) {
    sfree(fm);
    *fmatoms = NULL;
  } else {
    *fmatoms = fm;
  }
  return interactions_count;
}

/* for each atom in src, adds it's interactions to the same atom in dst */
void pf_atoms_merge(t_pf_atoms *dst, t_pf_atoms *src, int OnePair) {
  int i;

  if (src->len != dst->len)
    gmx_fatal(FARGS, "Mismatch in atoms list: %d vs. %d\n", src->len, dst->len);
  switch (OnePair) {
    case PF_ONEPAIR_DETAILED:
      for (i = 0; i < src->len; i++) {
        pf_atom_detailed_merge(&(dst->detailed[i]), &(src->detailed[i]));
      }
      break;
    case PF_ONEPAIR_SUMMED:
      for (i = 0; i < src->len; i++) {
        pf_atom_summed_merge(&(dst->summed[i]), &(src->summed[i]));
      }
      break;
  }
}

/* for each interaction, divide the force by divisor */
void pf_atoms_divide(t_pf_atoms *atoms, int OnePair, real divisor) {
  int i;

  switch (OnePair) {
    case PF_ONEPAIR_DETAILED:
      for (i = 0; i < atoms->len; i++) {
        pf_atom_detailed_real_divide(&(atoms->detailed[i]), divisor);
      }
      break;
    case PF_ONEPAIR_SUMMED:
      for (i = 0; i < atoms->len; i++) {
        pf_atom_summed_real_divide(&(atoms->summed[i]), divisor);
      }
      break;
  }
}
