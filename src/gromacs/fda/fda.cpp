/*
 * Handling of arrays of pairwise forces.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#include "fda.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "pf_array.h"
#include "pf_array_detailed.h"
#include "pf_array_summed.h"
#include "pf_interactions.h"
#include "pf_utils.h"

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

static const real THIRD   = 1.0 / 3.0;
static const real QUARTER = 0.25;

/* these arrays should correspond to the defines in include/fda.h */
/* they are used both for reading input from user and in output for confirmation */
const char *pf_file_out_atom_option[FILE_OUT_NR+1] = {
  "no",
  "pairwise_forces_vector",
  "pairwise_forces_scalar",
  "punctual_stress",
  "virial_stress",
  "virial_stress_von_mises",
  "compat_bin",
  "compat_ascii",
  NULL
};

const char *pf_file_out_residue_option[FILE_OUT_NR+1] = {
  "no",
  "pairwise_forces_vector",
  "pairwise_forces_scalar",
  "punctual_stress",
  "virial_stress",
  "virial_stress_von_mises",
  "compat_bin",
  "compat_ascii",
  NULL
};

const char *pf_residue_renumber_option[PF_RESIDUE_RENUMBER_NR+1] = {
  "auto",
  "yes",
  "no",
  NULL
};

const char *pf_vector2scalar_option[PF_VECTOR2SCALAR_NR+1] = {
  "norm",
  "projection",
  NULL
};

FDA::FDA(int nfile, const t_filenm fnm[], gmx_mtop_t *top_global)
 : bInitialized(false),
   AtomBased(0),
   ResidueBased(0),
   PFPS(0),
   VS(0),
   OnePair(0),
   Vector2Scalar(0),
   syslen_atoms(0),
   syslen_residues(0),
   sys_in_g1(nullptr),
   sys_in_g2(nullptr),
   atoms(nullptr),
   residues(nullptr),
   atom2residue(nullptr),
   ResiduesRenumber(0),
   type(0),
   ofn_atoms(nullptr),
   of_atoms(nullptr),
   ofn_residues(nullptr),
   of_residues(nullptr),
   groupname(nullptr),
   nsteps_atoms(0),
   nsteps_residues(0),
   time_averages(nullptr),
   per_atom_real(nullptr),
   per_atom_real_int(nullptr),
   per_residue_real(nullptr),
   per_residue_real_int(nullptr),
   no_end_zeros(false),
   atom_vir(nullptr)
{
  t_inpfile *inp;
  int ninp, i;
  const char *pf_file_in;
  const char *tmp;			/* needs to be const char * to avoid "warning assignment discards qualifiers from pointer target type" compiler messages */
  char g1name[STRLEN], g2name[STRLEN];	/* readinp.c::get_estr() function uses a char[32] buffer */
  char tmpstr[STRLEN];
  warninp_t wi;

  snew(atoms, 1);
  snew(residues, 1);

  /* check for the pf configuration file (specified with -pfi option);
   * if it doesn't exist, return NULL to specify that no pf handling is done;
   * otherwise, check also for specification of the index file (-pfn)
   */
  if (opt2bSet("-pfi", nfile, fnm)) {
	if (!opt2bSet("-pfn", nfile, fnm)) {
	  gmx_fatal(FARGS, "-pfi option (pairwise forces configuration) specified, an index file (-pfn) is also needed.\n");
	  }
	bInitialized = TRUE;
  } else {
	bInitialized = FALSE;
	gmx_fatal(FARGS, "No pairwise forces input file, no sense to compute pairwise forces.\n");
	return;
  }

  pf_file_in = opt2fn("-pfi", nfile, fnm);
  wi = init_warning(FALSE,0);
  /* use GROMACS function to read lines of form key = value from input file */
  inp = read_inpfile(pf_file_in, &ninp, wi);
/*  {
  int i;
  for (i=0; i < ninp; i++) {
	fprintf(stderr, "%d, %s, %s\n", i, inp[i].name, inp[i].value);
	}
  fprintf(stderr, "%s\n", get_estr(&ninp, &inp, "group1", "Protein"));
  fprintf(stderr, "%s\n", get_estr(&ninp, &inp, "group2", "Protein"));
  }
*/

  /* OnePair has to be initialized before the atoms/residues are initialized
   * because the data structures used for storing atoms/residues depend on it
   */
  OnePair = PF_ONEPAIR_SUMMED;
  STYPE("onepair", tmpstr, "summed");
  if (strcasecmp(tmpstr, "detailed") == 0)
	OnePair = PF_ONEPAIR_DETAILED;

  /* obtain in g1name and g2name the names of the 2 groups */
  STYPE("group1", g1name, "Protein");
  STYPE("group2", g2name, "Protein");
  /* check that there is an index file */
  if (!opt2bSet("-pfn", nfile, fnm))
	gmx_fatal(FARGS, "No index file (-pfn) for pairwise forces.\n");
  /* and finally read the index file and initialize the atoms lists */
  fprintf(stderr, "Pairwise forces for groups: %s and %s\n", g1name, g2name);

  STYPE("atombased", tmpstr, "no");
  AtomBased = FILE_OUT_NONE;
  for (i = 0; i < FILE_OUT_NR; i++)
	if (strcasecmp(tmpstr, pf_file_out_atom_option[i]) == 0)
	  AtomBased = i;

  STYPE("residuebased", tmpstr, "no");
  ResidueBased = FILE_OUT_NONE;
  for (i = 0; i < FILE_OUT_NR; i++)
	if (strcasecmp(tmpstr, pf_file_out_residue_option[i]) == 0)
	  ResidueBased = i;

  PFPS = pf_file_out_PF_or_PS(AtomBased) || pf_file_out_PF_or_PS(ResidueBased) ? 1 : 0;
  VS = pf_file_out_VS(AtomBased) || pf_file_out_VS(ResidueBased) ? 1 : 0;

  fprintf(stderr, "AtomBased: %s\n", pf_file_out_atom_option[AtomBased]);
  fprintf(stderr, "ResidueBased: %s\n", pf_file_out_residue_option[ResidueBased]);

  /* if using compatibility mode, there should be only one group */
  if ((pf_in_compatibility_mode(AtomBased)) || (pf_in_compatibility_mode(ResidueBased))) {
	if (gmx_strcasecmp(g1name, g2name))
	  gmx_fatal(FARGS, "When using compat mode, the two group names should the identical.\n");
	else
	  groupname = gmx_strdup(g1name);
  }

  if ((pf_file_out_stress(AtomBased) || pf_file_out_stress(ResidueBased)) && (OnePair != PF_ONEPAIR_SUMMED))
	gmx_fatal(FARGS, "Per atom data can only be computed from summed interactions.\n");

  STYPE("residuesrenumber", tmpstr, "auto");
  ResiduesRenumber = PF_RESIDUE_RENUMBER_AUTO;
  for (i = 0; i < PF_RESIDUE_RENUMBER_NR; i++)
	if (strcasecmp(tmpstr, pf_residue_renumber_option[i]) == 0)
	  ResiduesRenumber = i;
  fprintf(stderr, "ResidueRenumber: %s\n", pf_residue_renumber_option[ResiduesRenumber]);

  syslen_atoms = top_global->natoms;
  pf_fill_atom2residue(this, top_global);	/* also fills syslen_residues */

  if (AtomBased)
	atoms->sys2pf = pf_init_sys2pf(syslen_atoms);
  if (ResidueBased)
	residues->sys2pf = pf_init_sys2pf(syslen_residues);

  pf_read_group(this, opt2fn("-pfn", nfile, fnm), g1name, &sys_in_g1);
  pf_read_group(this, opt2fn("-pfn", nfile, fnm), g2name, &sys_in_g2);

  if (pf_file_out_PF_or_PS(AtomBased)) {
	pf_check_sys_in_g(this);
	pf_atoms_alloc(OnePair, atoms, syslen_atoms, "atoms");
	if (AtomBased == FILE_OUT_PUNCTUAL_STRESS) {
		pf_per_atom_real_init(&(per_atom_real), syslen_atoms, 0.0);
		ofn_atoms = gmx_strdup(opt2fn("-psa", nfile, fnm));
	} else {
		ofn_atoms = gmx_strdup(opt2fn("-pfa", nfile, fnm));
	}
	if (!ofn_atoms)
	  gmx_fatal(FARGS, "No file for writing out atom-based values.\n");
  }

  if (pf_file_out_PF_or_PS(ResidueBased)) {
	pf_check_sys_in_g(this);
	pf_atoms_alloc(OnePair, residues, syslen_residues, "residues");
	if (ResidueBased == FILE_OUT_PUNCTUAL_STRESS) {
		pf_per_atom_real_init(&(per_residue_real), syslen_residues, 0.0);
		ofn_residues = gmx_strdup(opt2fn("-psr", nfile, fnm));
	} else {
		ofn_residues = gmx_strdup(opt2fn("-pfr", nfile, fnm));
	}
	if (!ofn_residues)
	  gmx_fatal(FARGS, "No file for writing out residue-based values.\n");
  }

  if (pf_file_out_VS(AtomBased)) {
	snew(atom_vir, top_global->natoms);
	if (AtomBased == FILE_OUT_VIRIAL_STRESS)
		ofn_atoms = gmx_strdup(opt2fn("-vsa", nfile, fnm));
	if (AtomBased == FILE_OUT_VIRIAL_STRESS_VON_MISES)
		ofn_atoms = gmx_strdup(opt2fn("-vma", nfile, fnm));
	if (!ofn_atoms) gmx_fatal(FARGS, "No file for writing out virial stress.\n");
  }

  STYPE("type", tmpstr, "all");
  type = pf_interactions_type_str2val(tmpstr);
  if (type == PF_INTER_NONE)
	gmx_fatal(FARGS, "No interactions selected, no sense to compute pairwise forces.\n");
  else {
	fprintf(stderr, "Pairwise interactions selected: %s\n", pf_interactions_type_val2str(type));
  }

  STYPE("vector2scalar", tmpstr, "norm");
  Vector2Scalar = PF_VECTOR2SCALAR_NORM;
  for (i = 0; i < PF_VECTOR2SCALAR_NR; i++)
	if (strcasecmp(tmpstr, pf_vector2scalar_option[i]) == 0)
	  Vector2Scalar = i;
  fprintf(stderr, "Vector2Scalar: %s\n", pf_vector2scalar_option[Vector2Scalar]);
  if (((pf_in_compatibility_mode(AtomBased)) || (pf_in_compatibility_mode(ResidueBased))) &&
	(Vector2Scalar != PF_VECTOR2SCALAR_NORM))
	  gmx_fatal(FARGS, "When using compat mode, pf_vector2scalar should be set to norm.\n");

  /* initialization of time averages */
  snew(time_averages, 1);
  ITYPE("time_averages_period", time_averages->period, 1);
  if (time_averages->period < 0)
	gmx_fatal(FARGS, "Invalid value for time_averages_period: %d\n", time_averages->period);
  time_averages->steps = 0;
  if (time_averages->period != 1) {
	if (OnePair != PF_ONEPAIR_SUMMED)
	  gmx_fatal(FARGS, "Can only save scalar time averages from summed interactions.\n");
	/* for atoms/residues, pf_atoms_init is called for each step from md.c, but for time averages the initialization needs to be done only at the begining and after every data writing */
	if (pf_file_out_PF_or_PS(AtomBased)) {
	  if (!((pf_in_compatibility_mode(AtomBased)) || (AtomBased == FILE_OUT_PAIRWISE_FORCES_SCALAR)))
		gmx_fatal(FARGS, "Can only use time averages with scalar or compatibility output.\n");
	  pf_atoms_scalar_alloc(atoms, syslen_atoms, "atoms - time averages");
	  pf_atoms_scalar_init(atoms);
	}
	if (pf_file_out_PF_or_PS(ResidueBased)) {
	  if (!((pf_in_compatibility_mode(ResidueBased)) || (ResidueBased == FILE_OUT_PAIRWISE_FORCES_SCALAR)))
		gmx_fatal(FARGS, "Can only use time averages with scalar or compatibility output.\n");
	  pf_atoms_scalar_alloc(residues, syslen_residues, "residues - time averages");
	  pf_atoms_scalar_init(residues);
	  snew(time_averages->com, syslen_residues);
	  clear_rvecs(syslen_residues, time_averages->com);
	}
  }

  STYPE("no_end_zeros", tmpstr, "no");
  no_end_zeros = FALSE;
  if (strcasecmp(tmpstr, "no") != 0)
	no_end_zeros = TRUE;
}

void pf_atom_add_bonded_nocheck(FDA *fda, int i, int j, int type, rvec force) {
  t_pf_atoms *atoms;
  t_pf_atoms *residues;
  int ri = 0, rj = 0;       /* initialized to get rid of compiler warning, as they are only initialized later if ResidueBased is non-zero */
  rvec force_residue;

  atoms = fda->atoms;
  residues = fda->residues;
  //fprintf(stderr, "pf_atom_add_bonded_nocheck: adding i=%d, j=%d, type=%d\n", i, j, type);

  /* checking is symmetrical for atoms i and j; one of them has to be from g1, the other one from g2;
   * the check below makes the atoms equivalent, make them always have the same order (i,j) and not (j,i) where i < j;
   * force is the force atom j exerts on atom i; if i and j are switched, the force goes in opposite direction
   * it's possible that i > j, but ri < rj, so the force has to be handled separately for each of them
   * 
   * the logic is a bit complicated by the fact that AtomBased and ResidueBased are independent;
   * if ResidueBased part is done first, the AtomBased part can use i/j/force directly, without saving them
   * first in intermediate variables, as the initial values of i/j/force are no longer needed; if AtomBased
   * is done first (the original code), i/j/force are needed for the later atom->residue mapping
   * and saving in intermediate variables is needed
   */
  if (pf_file_out_PF_or_PS(fda->ResidueBased)) {
    /* the calling functions will not have i == j, but there is not such guarantee for ri and rj;
     * and it makes no sense to look at the interaction of a residue to itself
     */
    ri = fda->atom2residue[i];
    rj = fda->atom2residue[j];
    //fprintf(stderr, "pf_atom_add_bonded_nocheck: i=%d, j=%d, ri=%d, rj=%d, type=%d\n", i, j, ri, rj, type);
    if (ri != rj) {
      switch(fda->OnePair) {
        case PF_ONEPAIR_DETAILED:
          if (ri > rj) {
            int_swap(&ri, &rj);
            clear_rvec(force_residue);
            rvec_dec(force_residue, force);
            pf_atom_detailed_add(&residues->detailed[residues->sys2pf[ri]], rj, type, force_residue);
          } else {
            pf_atom_detailed_add(&residues->detailed[residues->sys2pf[ri]], rj, type, force);
          }
          break;
        case PF_ONEPAIR_SUMMED:
          if (ri > rj) {
            int_swap(&ri, &rj);
            clear_rvec(force_residue);
            rvec_dec(force_residue, force);
            pf_atom_summed_add(&residues->summed[residues->sys2pf[ri]], rj, type, force_residue);
          } else {
            pf_atom_summed_add(&residues->summed[residues->sys2pf[ri]], rj, type, force);
          }
          break;
        default:
          break;
      }
    }
  }

  if (pf_file_out_PF_or_PS(fda->AtomBased)) {
    //fprintf(stderr, "pf_atom_add_bonded_nocheck: i=%d, j=%d, type=%d\n", i, j, type);
    if (i > j) {
      int_swap(&i, &j);
      rvec_opp(force);
    }
    switch(fda->OnePair) {
      case PF_ONEPAIR_DETAILED:
        pf_atom_detailed_add(&atoms->detailed[atoms->sys2pf[i]], j, type, force);
        break;
      case PF_ONEPAIR_SUMMED:
        pf_atom_summed_add(&atoms->summed[atoms->sys2pf[i]], j, type, force);
        break;
      default:
        break;
    }
  }
}

void FDA::add_bonded(int i, int j, int type, rvec force)
{
  // leave early if the interaction is not interesting
  if (!(this->type & type)) return;
  if (!this->atoms_in_groups(i, j)) {
	//fprintf(stderr, "Warning: Unneeded pair in pf_atom_add_bonded i=%i, j=%i\n", i, j); fflush(stderr);
    return;
  }
  pf_atom_add_bonded_nocheck(this, i, j, type, force);
}

void FDA::add_virial_angle(int ai, int aj, int ak,
  rvec r_ij, rvec r_kj, rvec f_i, rvec f_k)
{
  tensor v;
  v[XX][XX] = r_ij[XX] * f_i[XX] + r_kj[XX] * f_k[XX];
  v[YY][YY] = r_ij[YY] * f_i[YY] + r_kj[YY] * f_k[YY];
  v[ZZ][ZZ] = r_ij[ZZ] * f_i[ZZ] + r_kj[ZZ] * f_k[ZZ];
  v[XX][YY] = r_ij[XX] * f_i[YY] + r_kj[XX] * f_k[YY];
  v[XX][ZZ] = r_ij[XX] * f_i[ZZ] + r_kj[XX] * f_k[ZZ];
  v[YY][ZZ] = r_ij[YY] * f_i[ZZ] + r_kj[YY] * f_k[ZZ];
  pf_atom_virial_add(this, ai, v, THIRD);
  pf_atom_virial_add(this, aj, v, THIRD);
  pf_atom_virial_add(this, ak, v, THIRD);
}

void FDA::add_virial_dihedral(int i, int j, int k, int l,
  rvec f_i, rvec f_k, rvec f_l, rvec r_ij, rvec r_kj, rvec r_kl)
{
  rvec r_lj;
  tensor v;
  rvec_sub(r_kj, r_kl, r_lj);
  v[XX][XX] = r_ij[XX] * f_i[XX] + r_kj[XX] * f_k[XX] + r_lj[XX] * f_l[XX];
  v[YY][YY] = r_ij[YY] * f_i[YY] + r_kj[YY] * f_k[YY] + r_lj[YY] * f_l[YY];
  v[ZZ][ZZ] = r_ij[ZZ] * f_i[ZZ] + r_kj[ZZ] * f_k[ZZ] + r_lj[ZZ] * f_l[ZZ];
  v[XX][YY] = r_ij[XX] * f_i[YY] + r_kj[XX] * f_k[YY] + r_lj[XX] * f_l[YY];
  v[XX][ZZ] = r_ij[XX] * f_i[ZZ] + r_kj[XX] * f_k[ZZ] + r_lj[XX] * f_l[ZZ];
  v[YY][ZZ] = r_ij[YY] * f_i[ZZ] + r_kj[YY] * f_k[ZZ] + r_lj[YY] * f_l[ZZ];
  pf_atom_virial_add(this, i, v, QUARTER);
  pf_atom_virial_add(this, j, v, QUARTER);
  pf_atom_virial_add(this, k, v, QUARTER);
  pf_atom_virial_add(this, l, v, QUARTER);
}
