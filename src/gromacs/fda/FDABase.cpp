/*
 * FDABase.cpp
 *
 *  Created on: Nov 22, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include "FDABase.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/futil.h"
#include "PureInteractionType.h"
#include "Utilities.h"

namespace fda {

template <class Base>
FDABase<Base>::FDABase(ResultType result_type, int syslen, std::string const& result_filename, FDASettings const& fda_settings)
 : Base(result_type == ResultType::VIRIAL_STRESS or result_type == ResultType::VIRIAL_STRESS_VON_MISES, syslen),
	 result_type(result_type),
	 syslen(syslen),
	 distributed_forces(),
	 total_forces(PF_or_PS_mode() and result_type == ResultType::PUNCTUAL_STRESS ? syslen : 0, 0.0),
	 result_file(result_filename),
	 fda_settings(fda_settings)
{
  if (PF_or_PS_mode()) make_backup(result_filename.c_str());
  write_compat_header(1);
}

template <class Base>
void FDABase<Base>::write_frame(rvec *x, int nsteps)
{
  switch (fda_settings.one_pair) {
	case OnePair::DETAILED:
	  switch (result_type) {
		case ResultType::NO:
		  // do nothing
		  break;
		case ResultType::PAIRWISE_FORCES_VECTOR:
		  write_frame_detailed(x, true, nsteps);
		  break;
		case ResultType::PAIRWISE_FORCES_SCALAR:
		  write_frame_detailed(x, false, nsteps);
		  break;
		case ResultType::PUNCTUAL_STRESS:
		  gmx_fatal(FARGS, "Punctual stress is not supported for detailed output.\n");
		  break;
		case ResultType::VIRIAL_STRESS:
		  gmx_fatal(FARGS, "Virial stress is not supported for detailed output.\n");
		  break;
		case ResultType::VIRIAL_STRESS_VON_MISES:
		  gmx_fatal(FARGS, "Virial stress von Mises is not supported for detailed output.\n");
		  break;
		case ResultType::COMPAT_BIN:
		  gmx_fatal(FARGS, "Compatibility binary mode is not supported for detailed output.\n");
		  break;
		case ResultType::COMPAT_ASCII:
		  gmx_fatal(FARGS, "Compatibility ascii mode is not supported for detailed output.\n");
		  break;
		case ResultType::INVALID:
		  gmx_fatal(FARGS, "ResultType is invalid.\n");
		  break;
	  }
	  break;
	case OnePair::SUMMED:
	  switch (result_type) {
		case ResultType::NO:
		  // do nothing
		  break;
		case ResultType::PAIRWISE_FORCES_VECTOR:
		  write_frame_summed(x, true, nsteps);
		  break;
		case ResultType::PAIRWISE_FORCES_SCALAR:
		  write_frame_summed(x, false, nsteps);
		  break;
		case ResultType::PUNCTUAL_STRESS:
		  sum_total_forces(x);
		  write_total_forces();
		  break;
		case ResultType::VIRIAL_STRESS:
		  write_atom_virial_sum();
		  break;
		case ResultType::VIRIAL_STRESS_VON_MISES:
		  write_atom_virial_sum_von_mises();
		  break;
		case ResultType::COMPAT_BIN:
		  gmx_fatal(FARGS, "Compatibility binary mode is not supported for summed output.\n");
		  break;
		case ResultType::COMPAT_ASCII:
		  write_frame_atoms_summed_compat(x);
		  break;
		case ResultType::INVALID:
		  gmx_fatal(FARGS, "ResultType is invalid.\n");
		  break;
	  }
	  break;
	case OnePair::INVALID:
	  gmx_fatal(FARGS, "OnePair is invalid.\n");
	  break;
  }
}

template <class Base>
void FDABase<Base>::write_frame_detailed(rvec *x, bool print_vector, int nsteps)
{
  result_file << "frame " << nsteps << std::endl;
  for (auto const& v1 : distributed_forces.detailed) {
    int i = v1.first;
    for (auto const& v2 : v1.second) {
      int j = v2.first;
      for (int type = 0; type != to_index(PureInteractionType::NUMBER); ++type) {
    	Vector force = v2.second.force[type];
        if (print_vector) {
          result_file << i << " " << j << " " << force[XX] << " " << force[YY] << " " << force[ZZ] << " " << type << std::endl;
        } else {
          result_file << i << " " << j << " " << vector2signedscalar(force.get_pointer(), x[i], x[j], fda_settings.v2s) << " " << type << std::endl;
        }
      }
    }
  }
}

template <class Base>
void FDABase<Base>::write_frame_summed(rvec *x, bool print_vector, int nsteps)
{
  result_file << "frame " << nsteps << std::endl;
  for (auto const& v1 : distributed_forces.summed) {
	int i = v1.first;
	for (auto const& v2 : v1.second) {
	  int j = v2.first;
	  Vector force = v2.second.force;
	  if (print_vector) {
	    result_file << i << " " << j << " " << force[XX] << " " << force[YY] << " " << force[ZZ] << " " << v2.second.type << std::endl;
	  } else {
		result_file << i << " " << j << " " << vector2signedscalar(force.get_pointer(), x[i], x[j], fda_settings.v2s) << " " << v2.second.type << std::endl;
	  }
	}
  }
}

template <class Base>
void FDABase<Base>::write_frame_scalar(int nsteps)
{
  result_file << "frame " << nsteps << std::endl;
  for (auto const& v1 : distributed_forces.scalar) {
	int i = v1.first;
	for (auto const& v2 : v1.second) {
	  int j = v2.first;
      result_file << i << " " << j << " " << v2.second.force << " " << v2.second.type << std::endl;
	}
  }
}

template <class Base>
void FDABase<Base>::sum_total_forces(rvec *x)
{
  for (auto& v : total_forces) v = 0.0;

  for (auto const& v1 : distributed_forces.summed) {
    int i = v1.first;
    for (auto const& v2 : v1.second) {
      int j = v2.first;
      real scalar_force;
      switch (fda_settings.v2s) {
        case Vector2Scalar::NORM:
          scalar_force = norm(v2.second.force.get_pointer());
          break;
        case Vector2Scalar::PROJECTION:
          scalar_force = vector2unsignedscalar(v2.second.force.get_pointer(), i, j, x);
          break;
        default:
      	  gmx_fatal(FARGS, "Unknown option for Vector2Scalar.\n");
          break;
      }
      total_forces[i] += scalar_force;
      total_forces[j] += scalar_force;
    }
  }
}

template <class Base>
void FDABase<Base>::write_total_forces()
{
  int j = total_forces.size();
  // Detect the last non-zero item
  if (fda_settings.no_end_zeros) {
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
void FDABase<Base>::write_compat_header(int nsteps)
{
  if (!PF_or_PS_mode() or !compatibility_mode()) return;

  result_file << "<begin_block>" << std::endl;
  result_file << "; Forcemat version " << FDASettings::compat_fm_version << std::endl;
  result_file << "; Matrix containing pairwise forces." << std::endl;
  result_file << "; Matrix dimension " << distributed_forces.size() << " x " << distributed_forces.size() << std::endl;
  result_file << "version=" << FDASettings::compat_fm_version << std::endl;
  result_file << "groupname=" << fda_settings.groupname << std::endl;
  result_file << "writefreq=1" << std::endl;
  // Reserve place for up to 8 digits for nsteps
  result_file << "nsteps=" << std::setfill('0') << std::setw(8) << nsteps << std::endl;
  result_file << "sysanr=" << syslen << std::endl;
  result_file << "fmdim=" << distributed_forces.size() << std::endl;
  result_file << "intsize=" << sizeof(int) << std::endl;
  result_file << "realsize=" << sizeof(real) << std::endl;
  result_file << "<end_block>" << std::endl;
}

template <class Base>
void FDABase<Base>::write_frame_atoms_compat(int nsteps)
{
  throw std::runtime_error("Not implemented yet");
  if (result_type == ResultType::COMPAT_ASCII) {
    result_file << "<begin_block>" << std::endl;
    result_file << nsteps << std::endl;
//    result_file << interactions_count << std::endl;
//    pf_write_space_separated_int_array(f, fmatoms, atoms->len);
//    pf_write_space_separated_int_array(f, index, interactions_count);
//    pf_write_space_separated_real_array(f, force, interactions_count);
//    pf_write_space_separated_char_array(f, interaction, interactions_count);
    result_file << "<end_block>" << std::endl;
  } else {
//	fwrite(framenr, sizeof(*framenr), 1, f);
//	fwrite(&interactions_count, sizeof(interactions_count), 1, f);
//	fwrite(fmatoms, sizeof(fmatoms[0]), atoms->len, f);
//	fwrite(index, sizeof(index[0]), interactions_count, f);
//	fwrite(force, sizeof(force[0]), interactions_count, f);
//	fwrite(interaction, sizeof(interaction[0]), interactions_count, f);
//    i = PF_COMPAT_NEW_ENTRY;
//	fwrite(&i, sizeof(i), 1, f);
  }
}

template <class Base>
void FDABase<Base>::write_frame_atoms_scalar_compat()
{
  throw std::runtime_error("Not implemented yet");
//  t_pf_atom_scalar i_atom_scalar;
//  t_pf_interaction_scalar j_atom_scalar;
//  t_pf_interaction_array_scalar *ia_scalar;
//  int i, j, ii, jj, c_ix;
//  gmx_int64_t interactions_count;
//  int *index;
//  real *force;
//  char *interaction;
//  int *fmatoms;
//
//  interactions_count = pf_atoms_scalar_to_fm(atoms, &fmatoms);
//  if (interactions_count == 0)
//    return;
//
//  snew(index, interactions_count);
//  snew(force, interactions_count);
//  snew(interaction, interactions_count);
//
//  /* map between the arrays and a square matrix of size atoms->len;
//   * this assumes that g1 and g2 are identical because in the original PF implementation
//   * there was only one group and the file format was designed for only one group;
//   * this assumption is checked during pf_init
//   */
//  c_ix = 0;
//  for (i = 0; i < atoms->len; i++) {
//    i_atom_scalar = atoms->scalar[i];
//    ia_scalar = &i_atom_scalar.interactions;
//    /* i should be here the same as the atoms->sys2pf[i_atom_scalar.nr], but this might not hold when running in parallel */
//    ii = atoms->sys2pf[i_atom_scalar.nr];
//    for (j = 0; j < ia_scalar->len; j++) {
//      j_atom_scalar = ia_scalar->array[j];
//      jj = atoms->sys2pf[j_atom_scalar.jjnr];
//      //fprintf(stderr, "fm i=%d, j=%d, ii=%d, ii.len=%d, jj=%d, type=%d, fx=%e\n", i_atom_scalar.nr, j_atom_scalar.jjnr, ii, ia_scalar->len, jj, j_atom_scalar.type, j_atom_scalar.force);
//      /* try to simulate the half-matrix storage of the original PF implementation */
//      index[c_ix] = (ii > jj) ? (jj * atoms->len) + ii : (ii * atoms->len) + jj;
//      force[c_ix] = j_atom_scalar.force;
//      interaction[c_ix] = pf_compatibility_map_interaction(j_atom_scalar.type);
//      c_ix++;
//    }
//  }
//
//  pf_write_frame_atoms_compat(atoms, f, framenr, interactions_count, fmatoms, index, force, interaction, ascii);
//
//  sfree(fmatoms);
//  sfree(index);
//  sfree(force);
//  sfree(interaction);
//
//  (*framenr)++;
}

template <class Base>
void FDABase<Base>::write_frame_atoms_summed_compat(rvec *x)
{
  throw std::runtime_error("Not implemented yet");
//  t_pf_atom_summed i_atom_summed;
//  t_pf_interaction_summed j_atom_summed;
//  t_pf_interaction_array_summed *ia_summed;
//  int i, j, ii, jj, c_ix;
//  gmx_int64_t interactions_count;
//  int *index;
//  real *force;
//  char *interaction;
//  int *fmatoms;
//
//  interactions_count = pf_atoms_summed_to_fm(atoms, &fmatoms);
//  if (interactions_count == 0)
//    return;
//
//  snew(index, interactions_count);
//  snew(force, interactions_count);
//  snew(interaction, interactions_count);
//
//  /* map between the arrays and a square matrix of size atoms->len;
//   * this assumes that g1 and g2 are identical because in the original PF implementation
//   * there was only one group and the file format was designed for only one group;
//   * this assumption is checked during pf_init
//   */
//  c_ix = 0;
//  for (i = 0; i < atoms->len; i++) {
//    i_atom_summed = atoms->summed[i];
//    ia_summed = &i_atom_summed.interactions;
//    /* i should be here the same as the atoms->sys2pf[i_atom_summed.nr], but this might not hold when running in parallel */
//    ii = atoms->sys2pf[i_atom_summed.nr];
//    for (j = 0; j < ia_summed->len; j++) {
//      j_atom_summed = ia_summed->array[j];
//      jj = atoms->sys2pf[j_atom_summed.jjnr];
//      //fprintf(stderr, "fm i=%d, j=%d, ii=%d, ii.len=%d, jj=%d, type=%d, fx=%e\n", i_atom_summed.nr, j_atom_summed.jjnr, ii, ia_summed->len, jj, j_atom_summed.type, j_atom_summed.force[0]);
//      /* try to simulate the half-matrix storage of the original PF implementation */
//      index[c_ix] = (ii > jj) ? (jj * atoms->len) + ii : (ii * atoms->len) + jj;
//      force[c_ix] = pf_vector2signedscalar(j_atom_summed.force, x[i_atom_summed.nr], x[j_atom_summed.jjnr], Vector2Scalar);
//      interaction[c_ix] = pf_compatibility_map_interaction(j_atom_summed.type);
//      c_ix++;
//    }
//  }
//
//  pf_write_frame_atoms_compat(atoms, f, framenr, interactions_count, fmatoms, index, force, interaction, ascii);
//
//  sfree(fmatoms);
//  sfree(index);
//  sfree(force);
//  sfree(interaction);
//
//  (*framenr)++;
}

template <>
void FDABase<Atom>::write_atom_virial_sum()
{
  bool first = true;
  for (auto const& v : virial_stress) {
    if (!first) result_file << " ";
    else first = false;
    result_file << -v(XX, XX) << " " << -v(YY, YY) << " " << -v(ZZ, ZZ) << " "
 				<< -v(XX, YY) << " " << -v(XX, ZZ) << " " << -v(YY, ZZ);
  }
  result_file << std::endl;
}

template <>
void FDABase<Residue>::write_atom_virial_sum()
{}

real tensor_to_vonmises(Tensor t)
{
  real txy = t(XX, XX)-t(YY, YY);
  real tyz = t(YY, YY)-t(ZZ, ZZ);
  real tzx = t(ZZ, ZZ)-t(XX, XX);
  return std::sqrt(0.5 * (txy*txy + tyz*tyz + tzx*tzx
		           + 6 * (t(XX, YY)*t(XX, YY) + t(XX, ZZ)*t(XX, ZZ) + t(YY, ZZ)*t(YY, ZZ))));
}

template <>
void FDABase<Atom>::write_atom_virial_sum_von_mises()
{
  bool first = true;
  for (auto const& v : virial_stress) {
	if (!first) result_file << " ";
	else first = false;
	result_file << tensor_to_vonmises(v);
  }
  result_file << std::endl;
}

template <>
void FDABase<Residue>::write_atom_virial_sum_von_mises()
{}

/// template instantiation
template class FDABase<Atom>;
template class FDABase<Residue>;

} // namespace fda
