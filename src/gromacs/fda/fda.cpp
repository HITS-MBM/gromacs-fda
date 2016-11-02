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

static const real HALF    = 1.0 / 2.0;
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
 : atom_based_forces(),
   residue_based_forces(),
   bInitialized(false),
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

  // allocates and fills with -1 the indexing table real atom nr. to pf number
  if (AtomBased)
	atom_based_forces.forces.resize(syslen_atoms);
  if (ResidueBased)
	residue_based_forces.forces.resize(syslen_residues);

  this->read_group(opt2fn("-pfn", nfile, fnm), g1name, &sys_in_g1);
  this->read_group(opt2fn("-pfn", nfile, fnm), g2name, &sys_in_g2);

  if (pf_file_out_PF_or_PS(AtomBased)) {
	pf_check_sys_in_g(this);
	pf_atoms_alloc(OnePair, atom_based_forces, syslen_atoms, "atoms");
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
	pf_atoms_alloc(OnePair, residue_based_forces, syslen_residues, "residues");
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
	  pf_atoms_scalar_alloc(atom_based_forces, syslen_atoms, "atoms - time averages");
	  pf_atoms_scalar_init(atom_based_forces);
	}
	if (pf_file_out_PF_or_PS(ResidueBased)) {
	  if (!((pf_in_compatibility_mode(ResidueBased)) || (ResidueBased == FILE_OUT_PAIRWISE_FORCES_SCALAR)))
		gmx_fatal(FARGS, "Can only use time averages with scalar or compatibility output.\n");
	  pf_atoms_scalar_alloc(residue_based_forces, syslen_residues, "residues - time averages");
	  pf_atoms_scalar_init(residue_based_forces);
	  snew(time_averages->com, syslen_residues);
	  clear_rvecs(syslen_residues, time_averages->com);
	}
  }

  STYPE("no_end_zeros", tmpstr, "no");
  no_end_zeros = FALSE;
  if (strcasecmp(tmpstr, "no") != 0)
	no_end_zeros = TRUE;
}

void FDA::add_bonded_nocheck(int i, int j, int type, rvec force)
{
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
  if (pf_file_out_PF_or_PS(ResidueBased)) {
    /* the calling functions will not have i == j, but there is not such guarantee for ri and rj;
     * and it makes no sense to look at the interaction of a residue to itself
     */
    int ri = atom2residue[i];
    int rj = atom2residue[j];
    //fprintf(stderr, "pf_atom_add_bonded_nocheck: i=%d, j=%d, ri=%d, rj=%d, type=%d\n", i, j, ri, rj, type);
    rvec force_residue;
    if (ri != rj) {
      switch(OnePair) {
        case PF_ONEPAIR_DETAILED:
          if (ri > rj) {
            int_swap(&ri, &rj);
            clear_rvec(force_residue);
            rvec_dec(force_residue, force);
            pf_atom_detailed_add(&residue_based_forces->detailed[residue_based_forces->sys2pf[ri]], rj, type, force_residue);
          } else {
            pf_atom_detailed_add(&residue_based_forces->detailed[residue_based_forces->sys2pf[ri]], rj, type, force);
          }
          break;
        case PF_ONEPAIR_SUMMED:
          if (ri > rj) {
            int_swap(&ri, &rj);
            clear_rvec(force_residue);
            rvec_dec(force_residue, force);
            pf_atom_summed_add(&residue_based_forces->summed[residue_based_forces->sys2pf[ri]], rj, type, force_residue);
          } else {
            pf_atom_summed_add(&residue_based_forces->summed[residue_based_forces->sys2pf[ri]], rj, type, force);
          }
          break;
        default:
          break;
      }
    }
  }

  if (pf_file_out_PF_or_PS(AtomBased)) {
    //fprintf(stderr, "pf_atom_add_bonded_nocheck: i=%d, j=%d, type=%d\n", i, j, type);
    if (i > j) {
      int_swap(&i, &j);
      rvec_opp(force);
    }
    switch(OnePair) {
      case PF_ONEPAIR_DETAILED:
        pf_atom_detailed_add(&atom_based_forces->detailed[atom_based_forces->sys2pf[i]], j, type, force);
        break;
      case PF_ONEPAIR_SUMMED:
        pf_atom_summed_add(&atom_based_forces->summed[atom_based_forces->sys2pf[i]], j, type, force);
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
  add_bonded_nocheck(i, j, type, force);
}
void FDA::add_nonbonded_single(int i, int j, int type, real force, real dx, real dy, real dz)
{
  rvec force_v;			/* vector force for interaction */

  /* first check that the interaction is interesting before doing computation and lookups */
  if (!(type & type))
    return;
  if (!atoms_in_groups(i, j)) {
	//fprintf(stderr, "Warning: Unneeded pair in pf_atom_add_nonbonded_single i=%i, j=%i\n", i, j); fflush(stderr);
    return;
  }

  force_v[0] = force * dx;
  force_v[1] = force * dy;
  force_v[2] = force * dz;
  add_bonded_nocheck(i, j, type, force_v);
}

void FDA::add_nonbonded(int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz)
{
  int ri = 0, rj = 0;
  real pf_lj_residue, pf_coul_residue, pf_lj_coul;
  rvec pf_lj_atom_v, pf_lj_residue_v, pf_coul_atom_v, pf_coul_residue_v;

  /* first check that the interaction is interesting before doing expensive calculations and atom lookup*/
  if (!(type & PF_INTER_COULOMB))
    if (!(type & PF_INTER_LJ))
      return;
    else {
      add_nonbonded_single(i, j, PF_INTER_LJ, pf_lj, dx, dy, dz);
      return;
    }
  else
    if (!(type & PF_INTER_LJ)) {
      add_nonbonded_single(i, j, PF_INTER_COULOMB, pf_coul, dx, dy, dz);
      return;
    }
  if (!atoms_in_groups(i, j)) {
	//fprintf(stderr, "Warning: Unneeded pair in pf_atom_add_nonbonded i=%i, j=%i\n", i, j); fflush(stderr);
    return;
  }

  /*fprintf(stderr, "Nonbonded interaction: ii=%d, jnr=%d\n", ii, jnr);*/

  /* checking is symmetrical for atoms i and j; one of them has to be from g1, the other one from g2
   * however, if only ResidueBased is non-zero, atoms won't be initialized... so the conversion to residue numebers needs to be done here already;
   * the check below makes the atoms equivalent, make them always have the same order (i,j) and not (j,i) where i < j;
   * force is the force atom j exerts on atom i; if i and j are switched, the force goes in opposite direction
   * it's possible that i > j, but ri < rj, so the force has to be handled separately for each of them
   */
  if (pf_file_out_PF_or_PS(ResidueBased)) {
    /* the calling functions will not have i == j, but there is not such guarantee for ri and rj;
     * and it makes no sense to look at the interaction of a residue to itself
     */
    ri = atom2residue[i];
    rj = atom2residue[j];
    //fprintf(stderr, "pf_atom_add_nonbonded: i=%d, j=%d, ri=%d, rj=%d\n", i, j, ri, rj);
    if (ri != rj) {
      if (ri > rj) {
        int tmp = rj; rj = ri; ri = tmp; // swap
        pf_lj_residue = -pf_lj;
        pf_coul_residue = -pf_coul;
      } else {
        pf_lj_residue = pf_lj;
        pf_coul_residue = pf_coul;
      }
      /* for detailed interactions, it's necessary to calculate and add separately, but for summed this can be simplified */
      switch(OnePair) {
        case PF_ONEPAIR_DETAILED:
          pf_coul_residue_v[0] = pf_coul_residue * dx;
          pf_coul_residue_v[1] = pf_coul_residue * dy;
          pf_coul_residue_v[2] = pf_coul_residue * dz;
          pf_lj_residue_v[0] = pf_lj_residue * dx;
          pf_lj_residue_v[1] = pf_lj_residue * dy;
          pf_lj_residue_v[2] = pf_lj_residue * dz;
          pf_atom_detailed_add(&residue_based_forces->detailed[residue_based_forces->sys2pf[ri]], rj, PF_INTER_LJ, pf_lj_residue_v);
          pf_atom_detailed_add(&residue_based_forces->detailed[residue_based_forces->sys2pf[ri]], rj, PF_INTER_COULOMB, pf_coul_residue_v);
          break;
        case PF_ONEPAIR_SUMMED:
          pf_lj_coul = pf_lj_residue + pf_coul_residue;
          pf_coul_residue_v[0] = pf_lj_coul * dx;
          pf_coul_residue_v[1] = pf_lj_coul * dy;
          pf_coul_residue_v[2] = pf_lj_coul * dz;
          pf_atom_summed_add(&residue_based_forces->summed[residue_based_forces->sys2pf[ri]], rj, PF_INTER_LJ | PF_INTER_COULOMB, pf_coul_residue_v);
          break;
        default:
          break;
      }
    }
  }

  /* i & j as well as pf_lj & pf_coul are not used after this point, so it's safe to operate on their values directly */
  if (pf_file_out_PF_or_PS(AtomBased)) {
    //fprintf(stderr, "pf_atom_add_nonbonded: i=%d, j=%d\n", i, j);
    if (i > j) {
      int tmp = j; j = i; i = tmp; // swap
      pf_lj = -pf_lj;
      pf_coul = -pf_coul;
    }
    /* for detailed interactions, it's necessary to calculate and add separately, but for summed this can be simplified */
    switch(OnePair) {
      case PF_ONEPAIR_DETAILED:
        pf_coul_atom_v[0] = pf_coul * dx;
        pf_coul_atom_v[1] = pf_coul * dy;
        pf_coul_atom_v[2] = pf_coul * dz;
        pf_lj_atom_v[0] = pf_lj * dx;
        pf_lj_atom_v[1] = pf_lj * dy;
        pf_lj_atom_v[2] = pf_lj * dz;
        pf_atom_detailed_add(&atom_based_forces->detailed[atom_based_forces->sys2pf[i]], j, PF_INTER_LJ, pf_lj_atom_v);
        pf_atom_detailed_add(&atom_based_forces->detailed[atom_based_forces->sys2pf[i]], j, PF_INTER_COULOMB, pf_coul_atom_v);
        break;
      case PF_ONEPAIR_SUMMED:
        pf_lj_coul = pf_lj + pf_coul;
        pf_coul_atom_v[0] = pf_lj_coul * dx;
        pf_coul_atom_v[1] = pf_lj_coul * dy;
        pf_coul_atom_v[2] = pf_lj_coul * dz;
        pf_atom_summed_add(&atom_based_forces->summed[atom_based_forces->sys2pf[i]], j, PF_INTER_LJ | PF_INTER_COULOMB, pf_coul_atom_v);
        break;
      default:
        break;
    }
  }
}

void FDA::add_virial(int ai, tensor v, real s)
{
    atom_vir[ai][XX][XX] += s * v[XX][XX];
    atom_vir[ai][YY][YY] += s * v[YY][YY];
    atom_vir[ai][ZZ][ZZ] += s * v[ZZ][ZZ];
    atom_vir[ai][XX][YY] += s * v[XX][YY];
    atom_vir[ai][XX][ZZ] += s * v[XX][ZZ];
    atom_vir[ai][YY][ZZ] += s * v[YY][ZZ];
}

void FDA::add_virial_bond(int ai, int aj, real f, real dx, real dy, real dz)
{
	tensor v;
	v[XX][XX] = dx * dx * f;
	v[YY][YY] = dy * dy * f;
	v[ZZ][ZZ] = dz * dz * f;
	v[XX][YY] = dx * dy * f;
	v[XX][ZZ] = dx * dz * f;
	v[YY][ZZ] = dy * dz * f;
	add_virial(ai, v, HALF);
	add_virial(aj, v, HALF);
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
  add_virial(ai, v, THIRD);
  add_virial(aj, v, THIRD);
  add_virial(ak, v, THIRD);
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
  add_virial(i, v, QUARTER);
  add_virial(j, v, QUARTER);
  add_virial(k, v, QUARTER);
  add_virial(l, v, QUARTER);
}
