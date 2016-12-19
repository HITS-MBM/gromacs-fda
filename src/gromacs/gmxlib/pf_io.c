/*
 * I/O
 * 
 * Copyright Bogdan Costescu 2011-2013
 */

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/pf_array_detailed.h"
#include "gromacs/legacyheaders/pf_array_summed.h"
#include "gromacs/legacyheaders/pf_interactions.h"
#include "gromacs/legacyheaders/pf_per_atom.h"
#include "gromacs/legacyheaders/pf_utils.h"
#include "gromacs/legacyheaders/types/pf_array.h"
#include "gromacs/legacyheaders/types/pf_array_detailed.h"
#include "gromacs/legacyheaders/types/pf_array_scalar.h"
#include "gromacs/legacyheaders/types/pf_array_summed.h"
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
enum {pf_compat_iNone, pf_compat_iBond, pf_compat_iAngle, pf_compat_iDihedral, pf_compat_iPolar, pf_compat_iLJ, pf_compat_iCoul, pf_compat_iAll};
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
  atom_id i;

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
  atom_id i;

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

void pf_write_frame_atoms_compat(t_pf_atoms *atoms, FILE *f, int *framenr, int interactions_count, int *fmatoms, int *index, real *force, char *interaction, gmx_bool ascii) {
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

void pf_write_frame_atoms_summed_compat(t_pf_atoms *atoms, FILE *f, int *framenr, const rvec *x, gmx_bool ascii, int Vector2Scalar) {
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
rvec *pf_residues_com(t_pf_global *pf_global, gmx_mtop_t *top_global, const rvec *x) {
  int moltype_index, mol_index, d;
  atom_id i, atom_index, atom_global_index, residue_global_index;
  t_atoms *atoms;
  t_atom *atom_info;
  gmx_molblock_t *mb;
  rvec r;
  real *mass;
  rvec *com;

  snew(com, pf_global->syslen_residues);
  snew(mass, pf_global->syslen_residues);

  for (i = 0; i < pf_global->syslen_residues; i++) {
    mass[i] = 0.0;
    for (d = 0; d < DIM; d++)
      com[i][d] = 0.0;
  }

  atom_global_index = 0;
  for (moltype_index = 0; moltype_index < top_global->nmolblock; moltype_index++) {	//enum all molecule types
    mb = &top_global->molblock[moltype_index];
    for (mol_index = 0; mol_index < mb->nmol; mol_index++) {				//enum all instances of each molecule type
    atoms = &top_global->moltype[mb->type].atoms;
      for(atom_index = 0; atom_index < atoms->nr; atom_index++) {			//enum atoms
        if ((pf_global->sys_in_g1[atom_global_index]) || (pf_global->sys_in_g2[atom_global_index])) {
          residue_global_index = pf_global->atom2residue[atom_global_index];
          atom_info=&atoms->atom[atom_index];
          mass[residue_global_index] += atom_info->m;
          svmul(atom_info->m, x[atom_global_index], r);
          rvec_inc(com[residue_global_index], r);
        }
        atom_global_index++;
      }
    }
  } /* for (moltype_index = 0... */

  /* there might be residues with no interesting atoms, so their mass would be set to the initial 0.0;
   * float/double comparison can be tricky in general, but here it's used only to prevent division by 0
   */
  for (i = 0; i < pf_global->syslen_residues; i++) {
    if (mass[i] != 0.0)
      for (d = 0; d < DIM; d++)
        com[i][d] /= mass[i];
  }

  sfree(mass);
  return com;
}

void pf_write_frame_detailed(t_pf_atoms *atoms, FILE* f, int *framenr, const rvec *x, gmx_bool bVector, int Vector2Scalar) {
  atom_id i, ii, j, jj;
  t_pf_interaction_array_detailed *iad;
  t_pf_interaction_detailed id;
  int type;

  fprintf(f, "frame %ld\n", (long int)(*framenr));
  for (i = 0; i < atoms->len; i++) {
    ii = atoms->detailed[i].nr;
    for (type = 1; type < PF_INTER_ALL; type <<= 1) {
      iad = pf_interaction_array_by_type(&(atoms->detailed[i].interactions), type);
      for (j = 0 ; j < iad->len; j++) {
        id = iad->array[j];
        jj = id.jjnr;
        //fprintf(stderr, "i=%ld j=%ld type=%s\n", ii, jj, pf_interactions_type_val2str(type));
        if (bVector)
          fprintf(f, "%ld %ld %e %e %e %d\n", (long int)ii, (long int)jj, id.force[XX], id.force[YY], id.force[ZZ], type);
        else
          fprintf(f, "%ld %ld %e %d\n", (long int)ii, (long int)jj, pf_vector2signedscalar(id.force, x[ii], x[jj], Vector2Scalar), type);
      }
    }
  }

  (*framenr)++;
}

void pf_write_frame_summed(t_pf_atoms *atoms, FILE* f, int *framenr, const rvec *x, gmx_bool bVector, int Vector2Scalar) {
  atom_id i, ii, j, jj;
  t_pf_interaction_array_summed ias;
  t_pf_interaction_summed is;

  fprintf(f, "frame %ld\n", (long int)(*framenr));
  for (i = 0; i < atoms->len; i++) {
    ii = atoms->summed[i].nr;
    ias = atoms->summed[i].interactions;
    for (j = 0 ; j < ias.len; j++) {
      is = ias.array[j];
      jj = is.jjnr;
      //fprintf(stderr, "i=%ld j=%ld type=%s\n", ii, jj, pf_interactions_type_val2str(is.type));
      if (bVector)
        fprintf(f, "%ld %ld %e %e %e %d\n", (long int)ii, (long int)jj, is.force[XX], is.force[YY], is.force[ZZ], is.type);
      else
        fprintf(f, "%ld %ld %e %d\n", (long int)ii, (long int)jj, pf_vector2signedscalar(is.force, x[ii], x[jj], Vector2Scalar), is.type);
    }
  }

  (*framenr)++;
}

void pf_write_frame_scalar(t_pf_atoms *atoms, FILE* f, int *framenr) {
  atom_id i, ii, j, jj;
  t_pf_interaction_array_scalar ias;
  t_pf_interaction_scalar is;

  fprintf(f, "frame %ld\n", (long int)(*framenr));
  for (i = 0; i < atoms->len; i++) {
    ii = atoms->scalar[i].nr;
    ias = atoms->scalar[i].interactions;
    for (j = 0 ; j < ias.len; j++) {
      is = ias.array[j];
      jj = is.jjnr;
      //fprintf(stderr, "i=%ld j=%ld type=%s\n", ii, jj, pf_interactions_type_val2str(is.type));
      fprintf(f, "%ld %ld %e %d\n", (long int)ii, (long int)jj, is.force, is.type);
    }
  }

  (*framenr)++;
}

void pf_write_frame(t_pf_global *pf_global, const rvec *x,  gmx_mtop_t *top_global) {
  rvec *com = NULL;	/* for residue based summation */

  /* this function is called from md.c, needs to check whether PF was initialized... */
  if (!pf_global->bInitialized)
    return;

  if (pf_file_out_PF_or_PS(pf_global->ResidueBased))
    com = pf_residues_com(pf_global, top_global, x);

  switch (pf_global->OnePair) {
    case PF_ONEPAIR_DETAILED:
      switch (pf_global->AtomBased) {
        case FILE_OUT_NONE:
          break;
        case FILE_OUT_PAIRWISE_FORCES_VECTOR:
          pf_write_frame_detailed(pf_global->atoms, pf_global->of_atoms, &pf_global->nsteps_atoms, x, TRUE, pf_global->Vector2Scalar);
          break;
        case FILE_OUT_PAIRWISE_FORCES_SCALAR:
          pf_write_frame_detailed(pf_global->atoms, pf_global->of_atoms, &pf_global->nsteps_atoms, x, FALSE, pf_global->Vector2Scalar);
          break;
        case FILE_OUT_VIRIAL_STRESS:
          gmx_fatal(FARGS, "Per atom detailed output not supported for virial stress.\n");
          break;
        case FILE_OUT_VIRIAL_STRESS_VON_MISES:
          gmx_fatal(FARGS, "Per atom detailed output not supported for virial stress von mises.\n");
          break;
        case FILE_OUT_COMPAT_BIN:
          break;
        case FILE_OUT_COMPAT_ASCII:
          gmx_fatal(FARGS, "Can't map detailed pairwise interactions in compatibility mode.\n");
          break;
        default:
          gmx_fatal(FARGS, "Per atom detailed output not implemented.\n");
          break;
      }
      switch (pf_global->ResidueBased) {
        case FILE_OUT_NONE:
          break;
        case FILE_OUT_PAIRWISE_FORCES_VECTOR:
          pf_write_frame_detailed(pf_global->residues, pf_global->of_residues, &pf_global->nsteps_residues, (const rvec*)com, TRUE, pf_global->Vector2Scalar);
          break;
        case FILE_OUT_PAIRWISE_FORCES_SCALAR:
          pf_write_frame_detailed(pf_global->residues, pf_global->of_residues, &pf_global->nsteps_residues, (const rvec*)com, FALSE, pf_global->Vector2Scalar);
          break;
        case FILE_OUT_VIRIAL_STRESS:
          gmx_fatal(FARGS, "Per residue detailed output not supported for virial stress.\n");
          break;
        case FILE_OUT_VIRIAL_STRESS_VON_MISES:
          gmx_fatal(FARGS, "Per residue detailed output not supported for virial stress von mises.\n");
          break;
        case FILE_OUT_COMPAT_BIN:
          break;
        case FILE_OUT_COMPAT_ASCII:
          gmx_fatal(FARGS, "Can't map detailed pairwise interactions in compatibility mode.\n");
          break;
        default:
          gmx_fatal(FARGS, "Per residue detailed output not implemented.\n");
          break;
      }
      break; /* PF_ONEPAIR_DETAILED */
    case PF_ONEPAIR_SUMMED:
      switch (pf_global->AtomBased) {
        case FILE_OUT_NONE:
          break;
        case FILE_OUT_PAIRWISE_FORCES_VECTOR:
          pf_write_frame_summed(pf_global->atoms, pf_global->of_atoms, &pf_global->nsteps_atoms, x, TRUE, pf_global->Vector2Scalar);
          break;
        case FILE_OUT_PAIRWISE_FORCES_SCALAR:
          pf_write_frame_summed(pf_global->atoms, pf_global->of_atoms, &pf_global->nsteps_atoms, x, FALSE, pf_global->Vector2Scalar);
          break;
        case FILE_OUT_PUNCTUAL_STRESS:
          pf_per_atom_sum(pf_global->per_atom_real, pf_global->atoms->summed, pf_global->atoms->len, x, pf_global->Vector2Scalar);
          pf_per_atom_real_write_frame(pf_global->of_atoms, pf_global->per_atom_real->force, pf_global->per_atom_real->len, pf_global->no_end_zeros);
          break;
        case FILE_OUT_VIRIAL_STRESS:
          pf_write_atom_virial_sum(pf_global->of_atoms, pf_global->atom_vir, pf_global->syslen_atoms);
          break;
        case FILE_OUT_VIRIAL_STRESS_VON_MISES:
          pf_write_atom_virial_sum_von_mises(pf_global->of_atoms, pf_global->atom_vir, pf_global->syslen_atoms);
          break;
        case FILE_OUT_COMPAT_BIN:
          break;
        case FILE_OUT_COMPAT_ASCII:
          pf_write_frame_atoms_summed_compat(pf_global->atoms, pf_global->of_atoms, &pf_global->nsteps_atoms, x, (pf_global->AtomBased == FILE_OUT_COMPAT_ASCII), pf_global->Vector2Scalar);
          break;
        default:
          gmx_fatal(FARGS, "Per atom summed output not implemented.\n");
          break;
      }
      switch (pf_global->ResidueBased) {
        case FILE_OUT_NONE:
          break;
        case FILE_OUT_PAIRWISE_FORCES_VECTOR:
          pf_write_frame_summed(pf_global->residues, pf_global->of_residues, &pf_global->nsteps_residues, (const rvec*)com, TRUE, pf_global->Vector2Scalar);
          break;
        case FILE_OUT_PAIRWISE_FORCES_SCALAR:
          pf_write_frame_summed(pf_global->residues, pf_global->of_residues, &pf_global->nsteps_residues, (const rvec *)com, FALSE, pf_global->Vector2Scalar);
          break;
        case FILE_OUT_PUNCTUAL_STRESS:
          pf_per_atom_sum(pf_global->per_residue_real, pf_global->residues->summed, pf_global->residues->len, (const rvec*)com, pf_global->Vector2Scalar);
          pf_per_atom_real_write_frame(pf_global->of_residues, pf_global->per_residue_real->force, pf_global->per_residue_real->len, pf_global->no_end_zeros);
          break;
        case FILE_OUT_VIRIAL_STRESS:
          gmx_fatal(FARGS, "Per residue summed output not supported for virial stress.\n");
          break;
        case FILE_OUT_VIRIAL_STRESS_VON_MISES:
          gmx_fatal(FARGS, "Per residue summed output not supported for virial stress von mises.\n");
          break;
        case FILE_OUT_COMPAT_BIN:
          break;
        case FILE_OUT_COMPAT_ASCII:
          pf_write_frame_atoms_summed_compat(pf_global->residues, pf_global->of_residues, &pf_global->nsteps_residues, (const rvec*)com, (pf_global->ResidueBased == FILE_OUT_COMPAT_ASCII), pf_global->Vector2Scalar);
          break;
        default:
          gmx_fatal(FARGS, "Per residue summed output not implemented.\n");
          break;
      }
      break; /* PF_ONEPAIR_SUMMED */
    default: /* switch (pf_global->OnePair) */
      break;
  }

  if (pf_file_out_PF_or_PS(pf_global->ResidueBased)) sfree(com);
}

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

void pf_open(t_pf_global *pf_global)
{
  /* this function is called from md.c, needs to check whether PF was initialized... */
  if (!pf_global->bInitialized) return;

  if (pf_file_out_PF_or_PS(pf_global->AtomBased)) {
    make_backup(pf_global->ofn_atoms);
    pf_global->of_atoms = gmx_fio_fopen(pf_global->ofn_atoms, "w+");
    if (pf_in_compatibility_mode(pf_global->AtomBased))
      pf_write_compat_header(pf_global->atoms, pf_global->of_atoms, 1, pf_global->syslen_atoms, pf_global->atoms->len, pf_global->groupname);
  }

  if (pf_file_out_PF_or_PS(pf_global->ResidueBased)) {
    make_backup(pf_global->ofn_residues);
    pf_global->of_residues = gmx_fio_fopen(pf_global->ofn_residues, "w+");
    if (pf_in_compatibility_mode(pf_global->ResidueBased))
      pf_write_compat_header(pf_global->residues, pf_global->of_residues, 1, pf_global->syslen_residues, pf_global->residues->len, pf_global->groupname);
  }

  if (pf_file_out_VS(pf_global->AtomBased)) {
    make_backup(pf_global->ofn_atoms);
    pf_global->of_atoms = gmx_fio_fopen(pf_global->ofn_atoms, "w+");
  }

  if (pf_file_out_VS(pf_global->ResidueBased)) {
    make_backup(pf_global->ofn_residues);
    pf_global->of_atoms = gmx_fio_fopen(pf_global->ofn_atoms, "w+");
  }
}

void pf_close(t_pf_global *pf_global)\
{
  /* this function is called from md.c, needs to check whether PF was initialized... */
  if (!pf_global->bInitialized) return;

  if (pf_file_out_PF_or_PS(pf_global->AtomBased)) {
    if (pf_in_compatibility_mode(pf_global->AtomBased)) {
      fprintf(pf_global->of_atoms, "EOF");
      fflush(pf_global->of_atoms);
      /* need to write header again to contain the correct nr. of steps */
      rewind(pf_global->of_atoms);
      pf_write_compat_header(pf_global->atoms, pf_global->of_atoms, pf_global->nsteps_atoms, pf_global->syslen_atoms, pf_global->atoms->len, pf_global->groupname);
    }
    gmx_fio_fclose(pf_global->of_atoms);
  }

  if (pf_file_out_PF_or_PS(pf_global->ResidueBased)) {
    if (pf_in_compatibility_mode(pf_global->ResidueBased)) {
      fprintf(pf_global->of_residues, "EOF");
      fflush(pf_global->of_residues);
      /* need to write header again to contain the correct nr. of steps */
      rewind(pf_global->of_residues);
      pf_write_compat_header(pf_global->residues, pf_global->of_residues, pf_global->nsteps_residues, pf_global->syslen_residues, pf_global->residues->len, pf_global->groupname);
    }
    gmx_fio_fclose(pf_global->of_residues);
  }

  if (pf_file_out_VS(pf_global->AtomBased)) {
    gmx_fio_fclose(pf_global->of_atoms);
  }

  if (pf_file_out_VS(pf_global->ResidueBased)) {
    gmx_fio_fclose(pf_global->of_residues);
  }
}

/* for each atom in src, adds it's interactions to the same atom in dst */
void pf_atoms_merge(t_pf_atoms *dst, t_pf_atoms *src, int OnePair) {
  atom_id i;

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
  atom_id i;

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
  atom_id i;

  for (i = 0; i < atoms->len; i++)
    rvec_inc(pf_x[i], x[atoms->scalar[i].nr]);
}

void pf_x_real_div(rvec *pf_x, atom_id pf_x_len, real divisor) {
  atom_id i;

  for (i = 0; i < pf_x_len; i++)
    svdiv(divisor, pf_x[i]);
}

void pf_write_frame_atoms_scalar_compat(t_pf_atoms *atoms, FILE *f, int *framenr, gmx_bool ascii) {
  t_pf_atom_scalar i_atom_scalar;
  t_pf_interaction_scalar j_atom_scalar;
  t_pf_interaction_array_scalar *ia_scalar;
  int i, j, ii, jj, c_ix;
  gmx_int64_t interactions_count;
  int *index;
  real *force;
  char *interaction;
  int *fmatoms;

  interactions_count = pf_atoms_scalar_to_fm(atoms, &fmatoms);
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
    i_atom_scalar = atoms->scalar[i];
    ia_scalar = &i_atom_scalar.interactions;
    /* i should be here the same as the atoms->sys2pf[i_atom_scalar.nr], but this might not hold when running in parallel */
    ii = atoms->sys2pf[i_atom_scalar.nr];
    for (j = 0; j < ia_scalar->len; j++) {
      j_atom_scalar = ia_scalar->array[j];
      jj = atoms->sys2pf[j_atom_scalar.jjnr];
      //fprintf(stderr, "fm i=%d, j=%d, ii=%d, ii.len=%d, jj=%d, type=%d, fx=%e\n", i_atom_scalar.nr, j_atom_scalar.jjnr, ii, ia_scalar->len, jj, j_atom_scalar.type, j_atom_scalar.force);
      /* try to simulate the half-matrix storage of the original PF implementation */
      index[c_ix] = (ii > jj) ? (jj * atoms->len) + ii : (ii * atoms->len) + jj;
      force[c_ix] = j_atom_scalar.force;
      interaction[c_ix] = pf_compatibility_map_interaction(j_atom_scalar.type);
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

/* write scalar time averages; this is similar to pf_write_frame, except that time averages are used */
void pf_write_scalar_time_averages(t_pf_global *pf_global) {
  /* there are several cases:
   * steps == 0 - there was no frame saved at all or there are no frames saved since the last writing, so no writing needed
   * steps > 0, period == 0 - no frames written so far, but frames saved, so writing one frame
   * steps > 0, period != 0 - frames saved, so writing the last frame
   *
   * if this is called from pf_write_and_save... then steps is certainly > 0
   * period is certainly != 1, otherwise the averages functions are not called
   */

  /* this function is called from md.c, needs to check whether PF was initialized... */
  if (!pf_global->bInitialized)
    return;

  if (pf_global->time_averages->steps == 0)
    return;

  //fprintf(stderr, "steps=%d\n", pf_global->time_averages->steps);
  /* atoms */
  if (pf_file_out_PF_or_PS(pf_global->AtomBased)) {
    pf_atoms_scalar_real_divide(pf_global->atoms, (real)pf_global->time_averages->steps);
    if (pf_in_compatibility_mode(pf_global->AtomBased))
      pf_write_frame_atoms_scalar_compat(pf_global->atoms, pf_global->of_atoms, &pf_global->nsteps_atoms, (pf_global->AtomBased == FILE_OUT_COMPAT_ASCII));
    else
      pf_write_frame_scalar(pf_global->atoms, pf_global->of_atoms, &pf_global->nsteps_atoms);
    pf_atoms_scalar_init(pf_global->atoms);
  }
  /* residues */
  if (pf_file_out_PF_or_PS(pf_global->ResidueBased)) {
    pf_atoms_scalar_real_divide(pf_global->residues, (real)pf_global->time_averages->steps);
    pf_x_real_div(pf_global->time_averages->com, pf_global->syslen_residues, (real)pf_global->time_averages->steps);
    if (pf_in_compatibility_mode(pf_global->ResidueBased))
      pf_write_frame_atoms_scalar_compat(pf_global->residues, pf_global->of_residues, &pf_global->nsteps_residues, (pf_global->ResidueBased == FILE_OUT_COMPAT_ASCII));
    else
      pf_write_frame_scalar(pf_global->residues, pf_global->of_residues, &pf_global->nsteps_residues);
    pf_atoms_scalar_init(pf_global->residues);
    clear_rvecs(pf_global->syslen_residues, pf_global->time_averages->com);
  }

  pf_global->time_averages->steps = 0;
}

/* main function for scalar time averages; saves data and decides when to write it out
 * 
 * dealing with residues is more complicated because COMs have to be averaged over time;
 * it's not possible to average the atom positions and calculate COMs only once before
 * writing because for this pf_global->atoms would need to be initialized (to get the atom
 * number or to get the sys2ps mapping) which only happens when AtomBased is non-zero
 */
void pf_save_and_write_scalar_time_averages(t_pf_global *pf_global, const rvec *x,  gmx_mtop_t *top_global) {
  rvec *com;    /* for residue based summation */

  /* this function is called from md.c, needs to check whether PF was initialized... */
  if (!pf_global->bInitialized)
    return;

  //fprintf(stderr, "save and write scalar data, steps=%d\n", pf_global->time_averages->steps);
  /* first save the data */
  if (pf_file_out_PF_or_PS(pf_global->AtomBased))
    pf_atoms_summed_merge_to_scalar(pf_global->atoms, x, pf_global->Vector2Scalar);
  if (pf_file_out_PF_or_PS(pf_global->ResidueBased)) {
    com = pf_residues_com(pf_global, top_global, x);
    pf_atoms_summed_merge_to_scalar(pf_global->residues, (const rvec *)com, pf_global->Vector2Scalar);
    pf_x_inc(pf_global->residues, pf_global->time_averages->com, (const rvec *)com);
    sfree(com);
  }

  pf_global->time_averages->steps++;
  if ((pf_global->time_averages->period != 0) && (pf_global->time_averages->steps >= pf_global->time_averages->period))
    pf_write_scalar_time_averages(pf_global);
}
