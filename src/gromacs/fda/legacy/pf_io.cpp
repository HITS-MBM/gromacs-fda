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

static inline void pf_write_space_separated_int_array(FILE *f, int *x, int len) {
  int i;

  for (i = 0; i < len; i++)
    if (i == 0)
      fprintf(f, "%d", x[i]);
    else
      fprintf(f, " %d", x[i]);
  fprintf(f, "\n");
}

static inline void pf_write_space_separated_char_array(FILE *f, char *x, int len) {
  int i;

  for (i = 0; i < len; i++)
    if (i == 0)
      fprintf(f, "%d", x[i]);
    else
      fprintf(f, " %d", x[i]);
  fprintf(f, "\n");
}

static inline void pf_write_space_separated_real_array(FILE *f, real *x, int len) {
  int i;

  for (i = 0; i < len; i++)
    if (i == 0)
      fprintf(f, "%e", x[i]);
    else
      fprintf(f, " %e", x[i]);
  fprintf(f, "\n");
}

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
/* version of force matrix implementation */
#define PF_COMPAT_FM_VERSION "1.5"

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

void pf_write_frame_atoms_compat(DistributedForces const& forces, FILE *f, int *framenr, int interactions_count, int *fmatoms, int *index, real *force, char *interaction, gmx_bool ascii) {
  int i;

  if (ascii) {
    fprintf(f, "<begin_block>\n");
    fprintf(f, "%d\n", *framenr);
    fprintf(f, "%d\n", interactions_count);
    pf_write_space_separated_int_array(f, fmatoms, atoms->len);
    pf_write_space_separated_int_array(f, index, interactions_count);
    pf_write_space_separated_real_array(f, force, interactions_count);
    pf_write_space_separated_char_array(f, interaction, interactions_count);
    fprintf(f, "<end_block>\n");
  }
  else {
    fwrite(framenr, sizeof(*framenr), 1, f);
    fwrite(&interactions_count, sizeof(interactions_count), 1, f);
    fwrite(fmatoms, sizeof(fmatoms[0]), atoms->len, f);
    fwrite(index, sizeof(index[0]), interactions_count, f);
    fwrite(force, sizeof(force[0]), interactions_count, f);
    fwrite(interaction, sizeof(interaction[0]), interactions_count, f);
    i = PF_COMPAT_NEW_ENTRY;
    fwrite(&i, sizeof(i), 1, f);
  }
}

void pf_write_frame_atoms_summed_compat(DistributedForces const& forces, FILE *f, int *framenr, const rvec *x, gmx_bool ascii, int Vector2Scalar)
{
  t_pf_atom_summed i_atom_summed;
  t_pf_interaction_summed j_atom_summed;
  t_pf_interaction_array_summed *ia_summed;
  int i, j, ii, jj, c_ix;
  gmx_int64_t interactions_count;
  int *index;
  real *force;
  char *interaction;
  int *fmatoms;

  interactions_count = pf_atoms_summed_to_fm(atoms, &fmatoms);
  if (interactions_count == 0)
    return;

  snew(index, interactions_count);
  snew(force, interactions_count);
  snew(interaction, interactions_count);

  /* map between the arrays and a square matrix of size atoms->len;
   * this assumes that g1 and g2 are identical because in the original PF implementation
   * there was only one group and the file format was designed for only one group;
   * this assumption is checked during pf_init
   */
  c_ix = 0;
  for (i = 0; i < atoms->len; i++) {
    i_atom_summed = atoms->summed[i];
    ia_summed = &i_atom_summed.interactions;
    /* i should be here the same as the atoms->sys2pf[i_atom_summed.nr], but this might not hold when running in parallel */
    ii = atoms->sys2pf[i_atom_summed.nr];
    for (j = 0; j < ia_summed->len; j++) {
      j_atom_summed = ia_summed->array[j];
      jj = atoms->sys2pf[j_atom_summed.jjnr];
      //fprintf(stderr, "fm i=%d, j=%d, ii=%d, ii.len=%d, jj=%d, type=%d, fx=%e\n", i_atom_summed.nr, j_atom_summed.jjnr, ii, ia_summed->len, jj, j_atom_summed.type, j_atom_summed.force[0]);
      /* try to simulate the half-matrix storage of the original PF implementation */
      index[c_ix] = (ii > jj) ? (jj * atoms->len) + ii : (ii * atoms->len) + jj;
      force[c_ix] = pf_vector2signedscalar(j_atom_summed.force, x[i_atom_summed.nr], x[j_atom_summed.jjnr], Vector2Scalar);
      interaction[c_ix] = pf_compatibility_map_interaction(j_atom_summed.type);
      c_ix++;
    }
  }

  pf_write_frame_atoms_compat(atoms, f, framenr, interactions_count, fmatoms, index, force, interaction, ascii);
  
  sfree(fmatoms);
  sfree(index);
  sfree(force);
  sfree(interaction);

  (*framenr)++;
}

/* computes the COM for residues in system;
 * only the atoms for which sys_in_g is non-zero are considered, such that the COM might
 * not express the COM of the whole residue but the COM of the atoms of the residue which
 * are interesting for PF
 */
rvec *pf_residues_com(FDA *fda, gmx_mtop_t *top_global, const rvec *x)


/* writes a header as in original PF implementation;
 * as the original PF implementation calculated everything then wrote out everything,
 * nsteps was known at the time when the file was written; however, to make it work
 * when data is written as it comes (frame by frame), the header is written once when
 * the file is open, but with nsteps=1 and space for up to 8 digits; when the file is
 * closed, the header is written again, this time with the correct nsteps
 */
void pf_write_compat_header(t_pf_atoms *atoms, FILE* f, int nsteps, int syslen, int fmdim, char *groupname) {
    fprintf(f, "<begin_block>\n");
    fprintf(f, "; Forcemat version %s\n", PF_COMPAT_FM_VERSION);
    fprintf(f, "; Matrix containing pairwise forces.\n");
    fprintf(f, "; Matrix dimension %i x %i\n", atoms->len, atoms->len);
    fprintf(f, "version=%s\n", PF_COMPAT_FM_VERSION);
    fprintf(f, "groupname=%s\n", groupname);
    fprintf(f, "writefreq=%i\n", 1);
    /* reserve place for up to 8 digits for nsteps; this can be done with %08i or %8i
     * depending on whether the space should be reserved by filling with 0 or blanks
     */
    fprintf(f, "nsteps=%08i\n", nsteps);
    fprintf(f, "sysanr=%i\n", syslen);
    fprintf(f, "fmdim=%i\n", fmdim);
    fprintf(f, "intsize=%i\n", (int)sizeof(int));
    fprintf(f, "realsize=%i\n", (int)sizeof(real));
    fprintf(f, "<end_block>\n");
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

/* this will be used only together with scalar interactions, so it's safe to use atoms->scalar */
void pf_x_inc(t_pf_atoms *atoms, rvec *pf_x, const rvec *x) {
  int i;

  for (i = 0; i < atoms->len; i++)
    rvec_inc(pf_x[i], x[atoms->scalar[i].nr]);
}

void pf_x_real_div(rvec *pf_x, int pf_x_len, real divisor) {
  int i;

  for (i = 0; i < pf_x_len; i++)
    svdiv(divisor, pf_x[i]);
}
