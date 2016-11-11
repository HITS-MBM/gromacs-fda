/*
 * Utils
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#include <stdio.h>

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
#include "pf_array.h"
#include "pf_array_detailed.h"
#include "pf_array_scalar.h"
#include "pf_array_summed.h"
#include "pf_interactions.h"
#include "pf_per_atom.h"
#include "pf_utils.h"
#include "types/pf_array_scalar.h"
#include "types/pf_array_summed.h"

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

t_pf_int_list *pf_int_list_alloc(int len) {
  t_pf_int_list *p;

  snew(p, 1);
  p->len = len;
  snew(p->list, p->len);
  return p;
}

void pf_int_list_free(t_pf_int_list *p) {
  sfree(p->list);
  sfree(p);
}

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

/* fills an indexing table: real atom nr. to pf number based on data from group atom lists*/
void pf_fill_sys2pf(int *sys2pf, int *len, t_pf_int_list *p) {
  int i;

  for (i = 0; i < p->len; i++) {
    if (sys2pf[p->list[i]] == -1) {
      sys2pf[p->list[i]] = *len;
      (*len)++;
    }
  }
}

char *pf_make_sys_in_group(int syslen, t_pf_int_list *p) {
  char *sys_in_group;
  int i;

  if (p->len == 0)
    return NULL;
  snew(sys_in_group, syslen);
  for (i = 0; i < p->len; i++) {
    sys_in_group[p->list[i]] = 1;
  }
  return sys_in_group;
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

/* returns the global residue number - based on mtop_util.c::gmx_mtop_atominfo_global(),
 * equivalent to a call to it with mtop->maxres_renum = INT_MAX */
int pf_get_global_residue_number(const gmx_mtop_t *mtop, const int atnr_global) {
  int mb;
  int a_start, a_end, maxresnr, at_loc;
  t_atoms *atoms=NULL;

  mb = -1;
  a_end = 0;
  maxresnr = 0;

  /* find the molecule block to which this atom belongs; maxresnr counts the nr. of residues encountered along the way */
  do
  {
    if (mb >= 0)
          maxresnr += mtop->molblock[mb].nmol*atoms->nres;
    mb++;
    atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
    a_start = a_end;
    a_end = a_start + mtop->molblock[mb].nmol*atoms->nr;
  }
  while (atnr_global >= a_end);

  /* get the location of the atom inside the current block;
   * as the block could appear repeated multiple times (repetition count stored in mtop->molblock[mb].nmol),
   * a simple substraction is not enough, modulo the total nr. of atoms in this block needs to be performed afterwards */
  at_loc = (atnr_global - a_start) % atoms->nr;
  /* (atnr_global - a_start)/atoms = in which repetition of the molblock this atom is found;
   * above *atoms->nres = the residues corresponding to these repetitions;
  * atoms->atom[at_loc].resind = residue index from the start of the molblock;
   * original function had maxresnr + 1 + (...), to make residue numbering start from 1 */
  return maxresnr + (atnr_global - a_start)/atoms->nr*atoms->nres + atoms->atom[at_loc].resind;
}

void pf_fill_atom2residue(FDA *fda, gmx_mtop_t *top_global) {
  int moltype_index, mol_index, atom_index, atom_global_index, resnrmax, i;
  t_atoms *atoms;
  t_atom *atom_info;
  gmx_molblock_t *mb;
  int *a2r_resnr, *a2r_renum;       /* atom 2 residue correspondence tables, both are filled, one will be used in the end */
  int *resnr2renum;                 /* residue nr. to renumbered nr. corespondence */
  int resnr;                            /* residue nr.; set to atoms->resinfo[].nr => type int */
  int renum;                        /* renumbered residue nr.; increased monotonically, so could theoretically be as large as the nr. of atoms => type int */
  gmx_bool bResnrCollision = FALSE;     /* TRUE if residue number collision */

  /* fill in the map between atom nr. and residue nr.; adapted from the code fragment in Data_Structures page on GROMACS website */
  snew(a2r_resnr, top_global->natoms);
  snew(a2r_renum, top_global->natoms);
  /* get the maximum number of residues in any molecule block;
   * any residue nr. obtained via atoms->resinfo[].nr is guaranteed to be smaller than this;
   * this will be used for residue nr. collision detection - initialized to -1,
   * will be set to the renumbered value on the first encounter, then a new renumbered value means a collision */
  resnrmax = gmx_mtop_maxresnr(top_global, 0);
  snew(resnr2renum, resnrmax + 1);
  //fprintf(stderr, "resnrmax=%d\n",resnrmax);
  for (i = 0; i <= resnrmax; i++)
    resnr2renum[i] = -1;
  atom_global_index = 0;
  renum = 0;
  for (moltype_index = 0; moltype_index < top_global->nmolblock; moltype_index++) {	//enum all molecule types
    mb = &top_global->molblock[moltype_index];
    for (mol_index = 0; mol_index < mb->nmol; mol_index++) {				//enum all instances of each molecule type
      atoms = &top_global->moltype[mb->type].atoms;
      for(atom_index = 0; atom_index < atoms->nr; atom_index++) {			//enum atoms
        atom_info=&atoms->atom[atom_index];
        resnr = atoms->resinfo[atom_info->resind].nr;
        renum = pf_get_global_residue_number(top_global, atom_global_index);
        //fprintf(stderr, "atom=%d, resnr=%d, renum=%d\n", atom_global_index, resnr, renum);
        //fprintf(stderr, "pair %d:%d\n", resnr, resnr2renum[resnr]);
        if ((resnr2renum[resnr] != renum) && (resnr2renum[resnr] != -1)) {
          bResnrCollision = TRUE;
          //fprintf(stderr, "Renumbering because %d:%d should be %d:%d\n", resnr, renum, resnr, resnr2renum[resnr]);
        }
        resnr2renum[resnr] = renum;
        a2r_resnr[atom_global_index] = resnr;
        a2r_renum[atom_global_index] = renum;
        atom_global_index++;
      }
    }
  } /* for (moltype_index = 0... */
  sfree(resnr2renum);

  /* renum is set to the residue number of the last atom, so should be largest value so far */
  switch (fda->ResiduesRenumber) {
    case PF_RESIDUE_RENUMBER_AUTO:
      if (bResnrCollision) {
        fprintf(stderr, "Residue number collision detected, residues will be renumbered.\n");
        fda->atom2residue = a2r_renum;
        fda->syslen_residues = renum + 1;
        sfree(a2r_resnr);
      }
      else {
        fprintf(stderr, "Residues will not be renumbered.\n");
        fda->atom2residue = a2r_resnr;
        fda->syslen_residues = resnrmax + 1;
        sfree(a2r_renum);
      }
      break;
    case PF_RESIDUE_RENUMBER_DO:
      fprintf(stderr, "Residues will be renumbered.\n");
      fda->atom2residue = a2r_renum;
      fda->syslen_residues = renum + 1;
      sfree(a2r_resnr);
      break;
    case PF_RESIDUE_RENUMBER_DONT:
      if (bResnrCollision)
        fprintf(stderr, "Residue number collision detected, residues will NOT be renumbered.\n");
      else
        fprintf(stderr, "Residues will not be renumbered.\n");
      fda->atom2residue = a2r_resnr;
      fda->syslen_residues = resnrmax + 1;
      sfree(a2r_renum);
      break;
  }
  //fprintf(stderr, "fda->syslen_residues=%d\n", fda->syslen_residues);
}

/* makes a list from the atoms numbers in the group
 * 
 * this can also be used later to map atoms to different nodes in a parallel run
 */
t_pf_int_list *pf_group2atoms(int len, int *list){
  t_pf_int_list *p;
  int i;

  if (len == 0)
    return NULL;
  p = pf_int_list_alloc(len);
  for (i = 0; i < p->len; i++)
    p->list[i] = list[i];
  return p;
}

/* makes a list of residue numbers based on atom numbers of this group;
 * this is slightly more complex than needed to allow the residue numbers to retain the ordering given to atoms
 */
t_pf_int_list *pf_groupatoms2residues(t_pf_int_list *atoms, FDA *fda) {
  int i, r;
  int nr_residues_in_g = 0;
  int *group_residues;
  gmx_bool *in_g;
  t_pf_int_list *p;

  /* as the list length is unknown at this point, make it as large as possible */
  snew(group_residues, fda->syslen_residues);
  /* store information about being already put in the list, to avoid lookup in a loop */
  snew(in_g, fda->syslen_residues);
  for (i = 0; i < fda->syslen_residues; i++)
    in_g[i] = FALSE;
  
  for (i = 0; i < atoms->len; i++) {
    r = fda->atom2residue[atoms->list[i]];
    if (!in_g[r]) {
      group_residues[nr_residues_in_g] = r;
      nr_residues_in_g++;
      in_g[r] = TRUE;
    }
  }

  /* as the length is now known, copy the residue ids to a properly sized list */
  p = pf_int_list_alloc(nr_residues_in_g);
  for (i = 0; i < p->len; i++)
    p->list[i] = group_residues[i];
  sfree(in_g);
  sfree(group_residues);
  return p;
}

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
