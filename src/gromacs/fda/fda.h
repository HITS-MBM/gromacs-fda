/*
 * Arrays of pairwise forces. This is a compromise between memory required to
 * store the pairwise forces and the CPU time needed to access them (insert,
 * lookup).
 * For each atom i there is a struct containing several arrays; there are
 * separate arrays for each interaction type, so there is no need to keep an
 * extra 'type' field. Each interaction takes space for one atom index and a
 * real vector containing the value of the force.
 * The arrays are grown as needed, through memory reallocations, similar to the
 * nonbonded lists; to avoid too many reallocations (memory management can be
 * expensive especially with certain MPI implementations which register memory
 * to achieve zero copy), they are always done by increasing the current size
 * with a factor > 1. However, this factor should be kept small enough to
 * avoid allocating too much memory which will then remain unused.
 * For each array there is a counter for the number of items currently in use
 * and one for the number of items for which space was allocated. If they are
 * equal and there is a request for storing another interaction, a reallocation
 * will happen.
 * Lookup is done through a linear search in the array. This uses CPU time, but
 * saves memory - organizing the items in a binary tree would require 2
 * pointers per interaction and allocations cannot be done in bulk anymore.
 * Insertion can be done directly at the end if it's known that the jjnr does
 * not already exist in the list; otherwise a linear search is needed.
 *
 * For the case of summed up forces per pair, a different set of data structures
 * are needed; the main reason is the extra data present in the summed up
 * interaction structure (the type) - it makes no sense to use a single structure
 * and ignore this field for the detailed storage of interactions, as this will
 * increase significantly the storage requirements. Furthermore, when storing
 * the atoms in group1, it makes no sense to keep all those NULL pointers
 * towards lists of interactions.
 * 
 * NOTE: all structures are defined for atoms, but are then used for both atoms and residues;
 * it made no sense to define new ones as they are functionally the same
 * 
 * Copyright Bogdan Costescu 2010-2013
 */

#ifndef GMX_FDA_FDA_H
#define GMX_FDA_FDA_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "pf_array_detailed.h"
#include "pf_array_scalar.h"
#include "pf_array_summed.h"
#include "pf_per_atom.h"

#ifdef __cplusplus
#include <cstdio>
#include <vector>
#include "DistributedForces.h"
#include "FDASettings.h"
#endif

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif

namespace fda {

/* type for a list of int's, containing both the list itself and the
 * length of the list
 */
typedef struct {
  int *list;
  int len;
} t_pf_int_list;

typedef struct {
  int period;       /* nr. of steps to average before writing; if 1 (default), no averaging is done; if 0, averaging is done over all steps so only one frame is written at the end */
  int steps;        /* counter for current step, incremented for every call of pf_save_and_write_scalar_averages(); when it reaches time_averages_steps, data is written */
  rvec *com;        /* averaged residue COM coordinates over steps, needed for COM calculations; only initialized when fda->ResidueBased is non-zero */
} t_pf_time_averages;

#ifdef __cplusplus
class FDA {
public:

  /// Default constructor
  FDA(FDASettings const& fda_settings = FDASettings());

  /// Returns true if atoms i and j are in fda groups
  gmx_bool atoms_in_groups(int i, int j) const {
	return ((sys_in_g1[i] && sys_in_g2[j]) || (sys_in_g1[j] && sys_in_g2[i]));
  }

  void add_bonded_nocheck(int i, int j, int type, rvec force);

  void add_bonded(int i, int j, int type, rvec force);

  /**
   *  Add a particular type of nonbonded interaction for the kernels where only one type of interaction is computed;
   *  force is passed as scalar along with the distance vector (as dx, dy, dz) from which the vector force is
   *  computed, the same way it's done in the nonbonded kernels
   */
  void add_nonbonded_single(int i, int j, int type, real force, real dx, real dy, real dz);

  /** Add a nonbonded interaction for kernels where both Coulomb and LJ are computed;
   *  this is more efficient than calling the previous one twice because some of the tests are made only once;
   *  forces are passed as scalars along with the distance vector (as dx, dy, dz) from which the vector forces are
   *  computed, the same way it's done in the nonbonded kernels
   */
  void add_nonbonded(int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz);

  /**
   *  The atom virial can be expressed as a 6-real tensor, as it's symmetric.
   *  To avoid defining a new tensor type, the 9-real tensor is used instead.
   */
  void add_virial(int ai, tensor v, real s);

  /**
   *  Origin on j, but for 2 atoms it doesn't matter.
   */
  void add_virial_bond(int ai, int aj, real fbond, real dx, real dy, real dz);

  /**
   *  Translate to origin on the middle (j) atom:
   *  vir = ri*Fi + rj*Fj + rk*Fk
   *      = (ri-rj)*Fi + (rk-rj)*Fk
   *      = r_ij[dim1]*f_i[dim2] + r_kj[dim1]*f_k[dim2]
   */
  void add_virial_angle(int ai, int aj, int ak, rvec r_ij, rvec r_kj, rvec f_i, rvec f_k);

  /**
   *  Translate to origin on the second (j) atom:
   *  vir = ri*Fi + rj*Fj + rk*Fk + rl*Fl
   *      = (ri-rj)*Fi + (rk-rj)*Fk + (rl-rj)*Fl
   *      = (ri-rj)*Fi + (rk-rj)*Fk + ((rl-rk) + (rk-rj))*Fl
   *      = r_ij[dim1]*f_i[dim2] + r_kj[dim1]*f_k[dim2] + (r_kj-r_kl)[dim1]*f_l[dim2]
   */
  void add_virial_dihedral(int i, int j, int k, int l, rvec f_i, rvec f_k, rvec f_l, rvec r_ij, rvec r_kj, rvec r_kl);

  void write_frame(const rvec *x,  gmx_mtop_t *top_global) /* const */;

  void write_frame_detailed(DistributedForces const& forces, FILE* f, int *framenr, const rvec *x, gmx_bool bVector, int Vector2Scalar) /* const */;

  void write_frame_summed(DistributedForces const& forces, FILE* f, int *framenr, const rvec *x, gmx_bool bVector, int Vector2Scalar) /* const */;

  void write_frame_scalar(DistributedForces const& forces, FILE* f, int *framenr) /* const */;

  void write_frame_atoms_scalar_compat(DistributedForces const& forces, FILE *f, int *framenr, gmx_bool ascii) /* const */;

  /**
   * Main function for scalar time averages; saves data and decides when to write it out
   *
   * dealing with residues is more complicated because COMs have to be averaged over time;
   * it's not possible to average the atom positions and calculate COMs only once before
   * writing because for this fda->atoms would need to be initialized (to get the atom
   * number or to get the sys2ps mapping) which only happens when AtomBased is non-zero
   */
  void save_and_write_scalar_time_averages(rvec const *x, gmx_mtop_t *top_global);

  /**
   * Write scalar time averages; this is similar to pf_write_frame, except that time averages are used
   *
   * There are several cases:
   * steps == 0 - there was no frame saved at all or there are no frames saved since the last writing, so no writing needed
   * steps > 0, period == 0 - no frames written so far, but frames saved, so writing one frame
   * steps > 0, period != 0 - frames saved, so writing the last frame
   *
   * if this is called from pf_write_and_save... then steps is certainly > 0
   * period is certainly != 1, otherwise the averages functions are not called
   */
  void write_scalar_time_averages();

  // TODO: Remove following legacy c-functions:

  void atoms_and_residues_init();

  void open();

  void close();

private:

  /// Settings
  FDASettings const& fda_settings;

  /// Distributed forces per atom
  DistributedForces atom_based_forces;

  /// Distributed forces per residue
  DistributedForces residue_based_forces;

public:

  // Helper flag for pairwise forces or punctual stress
  int PFPS;

  // Helper flag for virial stress
  int VS;

  /* Vector2Scalar defines the way the force vector to scalar conversion is done:
   * value | value in fi file | meaning
   * 0     | norm             | takes the norm of the vector
   * 1     | projection       | takes the projection of the force on the direction of the 2 atoms
   * As the conversion affects all scalar writing modes (PF_FILE_OUT*), this is kept as a separate
   * setting rather than creating separates modes for norm and projection.
   */
  int syslen_residues;          /* maximum of residue nr. + 1; residue nr. doesn't have to be continuous, there can be gaps */
  char *sys_in_g1;              /* 0 if atom not in group1, 1 if atom in group1, length of syslen_atoms; always allocated */
  char *sys_in_g2;              /* 0 if atom not in group2, 1 if atom in group2, length of syslen_atoms; always allocated */

  int *atom2residue;        /* stores the residue nr. for each atom; array of length syslen; only initialized if ResidueBased is non-zero */

  int type;                     /* interaction types that are interesting, set based on input file; functions are supposed to test against this before calculating/storing data */

  // File handles for pairwise forces.
  char *ofn_atoms;              /* output file name for atoms if AtomBased is non-zero */
  FILE *of_atoms;               /* output file for atoms if AtomBased is non-zero */
  char *ofn_residues;           /* output file name for residues if ResidueBased is non-zero */
  FILE *of_residues;            /* outpuf file for residues if ResidueBased is non-zero */

  char *groupname;              /* name of group for output in compatibility mode */
  /* the following 2 should be gmx_int64_t, but the file format is defined with int... */
  int nsteps_atoms;             /* nr. of steps for output in compatibility mode; incremented for each step written during run, also used to write the total nr. of steps at the end */
  int nsteps_residues;          /* nr. of steps for output in compatibility mode; incremented for each step written during run, also used to write the total nr. of steps at the end */
  t_pf_time_averages *time_averages;    /* only initialized when time averages are calculated */
  t_pf_per_atom_real *per_atom_real;            /* only initialized if required by user */
  t_pf_per_atom_real_int *per_atom_real_int;    /* only initialized if required by user */
  t_pf_per_atom_real *per_residue_real;          /* only initialized if required by user */
  t_pf_per_atom_real_int *per_residue_real_int; /* only initialized if required by user */
  gmx_bool no_end_zeros;	/* if True, trim the line such that the zeros at the end are not written; if False (default), all per atom/residue data is written */
  tensor *atom_vir;

#else
struct FDA {
#endif
};

} // namespace fda

#endif // GMX_FDA_FDA_H
