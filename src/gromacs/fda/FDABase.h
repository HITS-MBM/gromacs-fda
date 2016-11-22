/*
 * FDABase.h
 *
 *  Created on: Nov 22, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_FDABASE_H_
#define SRC_GROMACS_FDA_FDABASE_H_

#include <iostream>
#include <fstream>
#include <vector>
#include "DistributedForces.h"
#include "gromacs/utility/fatalerror.h"
#include "OnePair.h"
#include "ResultType.h"
#include "Vector2Scalar.h"

namespace fda {

struct Atom {};
struct Residue {};

template <class Base>
class FDABase
{
public:

  FDABase(ResultType result_type, OnePair one_pair, int syslen, std::string const& result_filename, bool no_end_zeros, Vector2Scalar v2s)
   : result_type(result_type),
	 distributed_forces(result_type, one_pair, syslen),
	 total_forces(),
	 result_file(result_filename),
	 no_end_zeros(no_end_zeros),
	 one_pair(one_pair),
	 v2s(v2s)
  {
    if (PF_or_PS_mode()) {
	  if (result_type == ResultType::PUNCTUAL_STRESS) {
        total_forces.resize(syslen, 0.0);
      }
    }
  }

  bool compatibility_mode() const {
    return result_type == ResultType::COMPAT_BIN or result_type == ResultType::COMPAT_ASCII;
  }

  bool stress_mode() const {
    return result_type == ResultType::PUNCTUAL_STRESS or result_type == ResultType::VIRIAL_STRESS or result_type == ResultType::VIRIAL_STRESS_VON_MISES;
  }

  bool PF_or_PS_mode() const {
	return result_type == ResultType::PAIRWISE_FORCES_VECTOR or result_type == ResultType::PAIRWISE_FORCES_SCALAR or result_type == ResultType::PUNCTUAL_STRESS;
  }

  bool VS_mode() const {
	return result_type == ResultType::VIRIAL_STRESS or result_type == ResultType::VIRIAL_STRESS_VON_MISES;
  }

  void write_frame(rvec *x, int nsteps);

  void write_frame_detailed(rvec *x, bool bVector, int nsteps);

  void write_frame_summed(rvec *x, bool bVector, int nsteps);

  void write_frame_scalar(int nsteps);

  void write_total_forces();

  void write_frame_atoms_scalar_compat(DistributedForces const& forces, FILE *f, int *framenr, gmx_bool ascii);

  /**
   * Main function for scalar time averages; saves data and decides when to write it out
   *
   * dealing with residues is more complicated because COMs have to be averaged over time;
   * it's not possible to average the atom positions and calculate COMs only once before
   * writing because for this fda->atoms would need to be initialized (to get the atom
   * number or to get the sys2ps mapping) which only happens when AtomBased is non-zero
   */
  void save_and_write_scalar_time_averages(rvec const *x, gmx_mtop_t *mtop);

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

  void per_atom_real_write_frame(FILE *f, std::vector<real> const& force_per_node, bool no_end_zeros);

  void per_atom_sum(std::vector<real>& force_per_node, DistributedForces& forces, int atoms_len, const rvec *x, Vector2Scalar v2s);

  /// The stress is the negative atom_vir value.
  void write_atom_virial_sum(FILE *f, tensor *atom_vir, int natoms);

  /// For von Mises no negative values are needed, since all items are squared.
  void write_atom_virial_sum_von_mises(FILE *f, tensor *atom_vir, int natoms);

private:

  friend class FDA;

  /// Result type
  ResultType result_type;

  /// Distributed forces
  DistributedForces distributed_forces;

  /// Total force per atom
  std::vector<real> total_forces;

  /// Result file
  std::ofstream result_file;

  /// If True, trim the line such that the zeros at the end are not written.
  /// if False (default), all per atom/residue data is written.
  bool no_end_zeros;

  /// OnePair defines the way the interactions between the same pair of atoms are stored
  OnePair one_pair;

  /// Define conversion from vector to scalar
  Vector2Scalar v2s;

};

template <class Base>
void FDABase<Base>::write_frame(rvec *x, int nsteps)
{
  switch (one_pair) {
	case OnePair::DETAILED:
	  switch (result_type) {
		case ResultType::PAIRWISE_FORCES_VECTOR:
		  write_frame_detailed(x, true, nsteps);
		  break;
		case ResultType::PAIRWISE_FORCES_SCALAR:
		  write_frame_detailed(x, false, nsteps);
		  break;
		case ResultType::VIRIAL_STRESS:
		  gmx_fatal(FARGS, "Per atom detailed output not supported for virial stress.\n");
		  break;
		case ResultType::VIRIAL_STRESS_VON_MISES:
		  gmx_fatal(FARGS, "Per atom detailed output not supported for virial stress von mises.\n");
		  break;
		case ResultType::COMPAT_BIN:
		  break;
		case ResultType::COMPAT_ASCII:
		  gmx_fatal(FARGS, "Can't map detailed pairwise interactions in compatibility mode.\n");
		  break;
		default:
		  gmx_fatal(FARGS, "Per atom detailed output not implemented.\n");
		  break;
	  }
	  break;
	case OnePair::SUMMED:
	  switch (result_type) {
		case ResultType::PAIRWISE_FORCES_VECTOR:
		  write_frame_summed(x, true, nsteps);
		  break;
		case ResultType::PAIRWISE_FORCES_SCALAR:
		  write_frame_summed(x, false, nsteps);
		  break;
		case ResultType::PUNCTUAL_STRESS:
		  per_atom_sum(force_per_atom, atom_based_forces->summed, atom_based_forces->len, x, v2s);
		  per_atom_real_write_frame(of_atoms, force_per_atom, no_end_zeros);
		  break;
		case ResultType::VIRIAL_STRESS:
		  pf_write_atom_virial_sum(of_atoms, atom_vir, syslen_atoms);
		  break;
		case ResultType::VIRIAL_STRESS_VON_MISES:
		  pf_write_atom_virial_sum_von_mises(of_atoms, atom_vir, syslen_atoms);
		  break;
		case ResultType::COMPAT_BIN:
		  break;
		case ResultType::COMPAT_ASCII:
		  pf_write_frame_atoms_summed_compat(atom_based_forces, of_atoms, &nsteps_atoms, x, (AtomBased == FILE_OUT_COMPAT_ASCII), v2s);
		  break;
		default:
		  gmx_fatal(FARGS, "Per atom summed output not implemented.\n");
		  break;
	  }
	  break;
  }
}

template <class Base>
void FDABase<Base>::write_frame_detailed(rvec *x, bool bVector, int nsteps)
{
  result_file << "frame " << nsteps << std::endl;
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
          fprintf(f, "%ld %ld %e %d\n", (long int)ii, (long int)jj, vector2signedscalar(id.force, x[ii], x[jj], v2s), type);
      }
    }
  }
}

template <class Base>
void FDABase<Base>::write_frame_summed(rvec *x, bool bVector, int nsteps)
{
  result_file << "frame " << nsteps << std::endl;
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
        fprintf(f, "%ld %ld %e %d\n", (long int)ii, (long int)jj, vector2signedscalar(is.force, x[ii], x[jj], v2s), is.type);
    }
  }
}

template <class Base>
void FDABase<Base>::write_frame_scalar(int nsteps)
{
  result_file << "frame " << nsteps << std::endl;
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
}

template <class Base>
void FDABase<Base>::write_total_forces()
{
  int j = total_forces.size();
  // Detect the last non-zero item
  if (no_end_zeros) {
    for (; j > 0; --j)
      if (total_forces[j - 1] != 0.0)
        break;
  }

  // j holds the index of first zero item or the length of force
  bool first_on_line = true;
  for (int i = 0; i < j; ++i) {
    if (first_on_line) {
      result_file << total_forces[i];
      first_on_line = false;
    } else {
      result_file << " " << total_forces[i];
    }
  }
  result_file << std::endl;
}

template <class Base>
void FDABase<Base>::save_and_write_scalar_time_averages(const rvec *x, gmx_mtop_t *mtop)
{
  if (fda_settings.time_averaging_period != 1) {
    // First save the data
    if (atom_based_forces.PF_or_PS_mode())
      atom_based_forces.summed_merge_to_scalar(x, v2s);
    if (residue_based_forces.PF_or_PS_mode()) {
  	  rvec *com = pf_residues_com(this, top_global, x);
      residue_based_forces.summed_merge_to_scalar(com, v2s);
      pf_x_inc(residue_based_forces, time_averaging_com, com);
      sfree(com);
    }

    time_averages->steps++;
    if ((time_averages->period != 0) && (time_averages->steps >= time_averages->period))
      this->write_scalar_time_averages();
  } else {
    write_frame(x, top_global);
  }
}

template <class Base>
void FDABase<Base>::write_scalar_time_averages()
{
  if (time_averaging_steps == 0) return;
  if (fda_settings.time_averaging_period == 1) return;

  if (atom_based_forces.PF_or_PS_mode()) {
    atom_based_forces.scalar_real_divide(time_averaging_steps);
    if (atom_based_forces.compatibility_mode())
      this->write_frame_atoms_scalar_compat(atom_based_forces, of_atoms, &nsteps_atoms, (AtomBased == FILE_OUT_COMPAT_ASCII));
    else
      this->write_frame_scalar(atom_based_forces, of_atoms, &nsteps_atoms);
    pf_atoms_scalar_init(atom_based_forces);
  }

  if (residue_based_forces.PF_or_PS_mode()) {
    residue_based_forces.scalar_real_divide(time_averaging_steps);
    pf_x_real_div(time_averaging_com, fda_settings.syslen_residues, time_averaging_steps);
    if (residue_based_forces.compatibility_mode())
      this->write_frame_atoms_scalar_compat(residue_based_forces, of_residues, &nsteps_residues, (ResidueBased == FILE_OUT_COMPAT_ASCII));
    else
      this->write_frame_scalar(residue_based_forces, of_residues, &nsteps_residues);
    pf_atoms_scalar_init(residue_based_forces);
    clear_rvecs(syslen_residues, time_averages->com);
  }

  time_averaging_steps = 0;
}

template <class Base>
void FDABase<Base>::write_frame_atoms_scalar_compat(t_pf_atoms *atoms, FILE *f, int *framenr, gmx_bool ascii)
{
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

template <class Base>
void FDABase<Base>::per_atom_real_write_frame(std::vector<real> const& force_per_node) const
{
  int j = force_per_node.size();
  // Detect the last non-zero item
  if (no_end_zeros) {
    for (; j > 0; --j)
      if (force[j - 1] != 0.0)
        break;
  }

  // j holds the index of first zero item or the length of force
  bool first_on_line = true;
  for (i = 0; i < j; i++) {
    if (first_on_line) {
      ofs << force[i];
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

template <class Base>
void FDABase<Base>::per_atom_sum(std::vector<real>& force_per_node, DistributedForces& forces, int atoms_len, const rvec *x, Vector2Scalar v2s)
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
      switch (v2s) {
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

template <class Base>
void FDABase<Base>::write_atom_virial_sum(FILE *f, tensor *atom_vir, int natoms)
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

template <class Base>
void FDABase<Base>::write_atom_virial_sum_von_mises(FILE *f, tensor *atom_vir, int natoms)
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

} // namespace fda

#endif /* SRC_GROMACS_FDA_FDABASE_H_ */
