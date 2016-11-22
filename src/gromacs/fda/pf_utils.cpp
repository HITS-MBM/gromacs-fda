/*
 * Utils
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#include <stdio.h>
#include "fda.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

/* strings below should not contain capital letters, as they are compared against a tolower() of user provided input */
/* this can't be turned into a simple const char * array as the ones below because the PF_INTER_* values are not defined in enum */
t_pf_interactions_type pf_interactions_type[] = {
  {PF_INTER_BOND, "bond"},
  {PF_INTER_ANGLE, "angle"},
  {PF_INTER_DIHEDRAL, "dihedral"},
  {PF_INTER_POLAR, "polar"},
  {PF_INTER_COULOMB, "coulomb"},
  {PF_INTER_LJ, "lj"},
  {PF_INTER_NB14, "nb14"},
  {PF_INTER_BONDED, "bonded"},
  {PF_INTER_NONBONDED, "nonbonded"},
  {PF_INTER_ALL, "all"}
};
#define PF_INTERACTIONS_NR asize(pf_interactions_type)

void pf_atom_summed_merge_to_scalar(t_pf_atom_summed *src, t_pf_atom_scalar *dst, const rvec *x, int Vector2Scalar) {
  t_pf_interaction_array_summed *ia;
  t_pf_interaction_summed *i;
  int j;

  if (src->nr != dst->nr)
    gmx_fatal(FARGS, "Mismatch of PF atom id: summed %d vs scalar %d\n", src->nr, dst->nr);
  ia = &(src->interactions);
  for (j = 0; j < ia->len; j++) {
    i = &(ia->array[j]);
    pf_atom_scalar_add(dst, i->jjnr, i->type, pf_vector2signedscalar(i->force, x[src->nr], x[i->jjnr], Vector2Scalar));
  }
}

void pf_atoms_summed_merge_to_scalar(t_pf_atoms *atoms, const rvec *x, int Vector2Scalar) {
  int i;

  for (i = 0; i < atoms->len; i++) {
    pf_atom_summed_merge_to_scalar(&(atoms->summed[i]), &(atoms->scalar[i]), x, Vector2Scalar);
  }
}

/* takes a force vector v and returns its norm (magnitude) along with a sign
 * or the projection on the position vector formed by atoms i and j;
 * the sign will be negative (=attractive force) when the vector is in the same direction
 * as the position vector formed by atoms i and j and positive otherwise
 */
real pf_vector2signedscalar(const rvec v, const rvec xi, const rvec xj, int Vector2Scalar) {
  rvec r;
  real c;

  /* r is the position vector */
  rvec_sub(xj, xi, r);
  c = cos_angle(r, v);
  switch (Vector2Scalar) {
    case PF_VECTOR2SCALAR_NORM:
      /* based on the cos of the angle between vectors:
       * if it's positive, the vector v is oriented within -90 to 90 degress with respect to r - considered to be same direction = attractive = negative;
       * if it's negative, the vector v is oriented within 90 to 270 degrees with respect to r - considered to be in opposite direction = repulsvive = positive;
       * if it's zero, the vector v is oriented at exactly 90 or 270 degrees with respect to r - the force can be considered neither attractive nor repulsvive, so it's set to zero
       */
      //fprintf(stderr, "v2ss: cos=%f norm=%f\n", c, norm(v));
      //fprintf(stderr, "v2ss: cos=%f, norm=%f, v[0]=%f, v[1]=%f, v[2]=%f, r[0]=%f, r[1]=%f, r[2]=%f\n", c, norm(v), v[0], v[1], v[2], r[0], r[1], r[2]);
      //if (c == 1.0) fprintf(stderr, "v2ss: cos=%f, norm=%f, v[0]=%f, v[1]=%f, v[2]=%f\n", c, norm(v), v[0], v[1], v[2]);
      if (c > 0) {
        return -norm(v);
      } else if (c < 0) {
        return norm(v);
      } else {
        return 0.0;
      }
      break;
    case PF_VECTOR2SCALAR_PROJECTION:
      /* it's minus c to go along the norm based calculation above: positive cos = in the same direction = attractive = negative */
      return -c * norm(v);
      break;
    default:
      /* this should not happen - make it visible! */
      return GMX_FLOAT_MAX;
      break;
  }
}

/* residue nr. doesn't hold much importance in GROMACS, it's only stored to be used when writing back
 * a structure file (f.e. PDB). Because of this, there is no residue nr. equivalent of the atom nr.;
 * the residue nr. which is held in atoms->resinfo[].nr is just what was in the PDB that was read, f.e.
 * if the PDB had residues 2 and 3 defined (no 1), atoms->resinfo[].nr will also contain 2 and 3;
 * as a result, it's not possible to loop over residues with something like:
 * for (r = 0; r < residues; r++)
 * NOTE: due to use of a global atom index, the mapping of atoms to residues won't work with DD !!!
 *
 * Things are evern worse because, with the above scheme, there is no guarantee that the residue
 * number is unique! For doing operations over residues, such guarantee is essential; otherwise
 * data for 2 different residues but which happen to have the same nr. will be added together - this
 * can happen f.e. if the structure is too big and residue nr. has to wraparound or if there are
 * several chains for which residue numbering starts from 1. GROMACS includes a residue renumbering
 * facility, but this destroys the correspondence atom nr. to residue nr. the user expects.
 * The logic below implements the following scheme:
 * - get both the residue nr. and the renumbered nr. and store them as a pair
 * - as long as the pair appears again, everything is fine, the residue nr. can be used and user
 * doesn't see any difference
 * - if another pair appears, let the user know and switch to renumbered nr.
 *
 * in output, residue nr. will be kept, but renumbered nr. will be -1, so that it starts from zero;
 * this way it fits with VMD's resid vs. residue numbering
 * 
 * fills in fda->atom2residue and fda->syslen_residues
 */

/* inspired by src/gmxlib/readir.c::do_egp_flag() which reads multiple energy exclusion groups */
int pf_interactions_type_str2val(char *typestr) {
  /* The maximum number of energy group pairs would be MAXPTR*(MAXPTR+1)/2.
   * But since this is much larger than STRLEN, such a line can not be parsed.
   * The real maximum is the number of names that fit in a string: STRLEN/2.
   */
#define EGP_MAX (STRLEN/2)
  int nelem, i, j;
  int type = PF_INTER_NONE;
  char *names[EGP_MAX];

  nelem = str_nelem(typestr, EGP_MAX, names);
  for (i = 0; i < nelem; i++) {
    lowstring(names[i]);
    j = 0;
    while ((j < PF_INTERACTIONS_NR) && gmx_strcasecmp(names[i], pf_interactions_type[j].str))
      j++;
    if (j == PF_INTERACTIONS_NR)
      gmx_fatal(FARGS,"%s is not a recognized interaction type\n", names[i]);
    type |= pf_interactions_type[j].val;
  }
  return type;
}

/* this function will allocate memory in strdup(), so a memory leak will happen if used in a loop with no free() around */
char *pf_interactions_type_val2str(int type) {
  /* as there are only a few interaction types possible and they have short names, it's safe to assume that the whole will fit in STRLEN */
  char tmpstr[STRLEN], tmpstr2[STRLEN];
  int i;
  gmx_bool first = TRUE;

  tmpstr[0] = '\0';
  for (i = 0; i < PF_INTERACTIONS_NR; i++)
    if ((type & pf_interactions_type[i].val) == pf_interactions_type[i].val) {
      if (first) {
        strcat(tmpstr, pf_interactions_type[i].str);
        first = FALSE;
      } else {
        sprintf(tmpstr2, " %s", pf_interactions_type[i].str);
        strcat(tmpstr, tmpstr2);
      }
    }
  return strdup(tmpstr);
}
