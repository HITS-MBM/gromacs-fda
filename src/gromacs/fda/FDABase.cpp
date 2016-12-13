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
   distributed_forces(syslen, fda_settings),
   result_file(result_filename),
   fda_settings(fda_settings)
{
  result_file << std::scientific << std::setprecision(6);
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
		  write_total_forces(x);
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
  // Clear arrays for next frame
  distributed_forces.clear();
}

template <class Base>
void FDABase<Base>::write_frame_detailed(rvec *x, bool print_vector, int nsteps)
{
  result_file << "frame " << nsteps << std::endl;
  if (print_vector)
    distributed_forces.write_detailed_vector(result_file);
  else
    distributed_forces.write_detailed_scalar(result_file, x);
}

template <class Base>
void FDABase<Base>::write_frame_summed(rvec *x, bool print_vector, int nsteps)
{
  result_file << "frame " << nsteps << std::endl;
  if (print_vector)
    distributed_forces.write_summed_vector(result_file);
  else
    distributed_forces.write_summed_scalar(result_file, x);
}

template <class Base>
void FDABase<Base>::write_frame_scalar(int nsteps)
{
  result_file << "frame " << nsteps << std::endl;
  distributed_forces.write_scalar(result_file);
}

template <class Base>
void FDABase<Base>::write_total_forces(rvec *x)
{
  distributed_forces.write_total_forces(result_file, x);
}

template <class Base>
void FDABase<Base>::write_compat_header(int nsteps)
{
  if (!PF_or_PS_mode() or !compatibility_mode()) return;

  result_file << "<begin_block>" << std::endl;
  result_file << "; Forcemat version " << FDASettings::compat_fm_version << std::endl;
  result_file << "; Matrix containing pairwise forces." << std::endl;
  result_file << "; Matrix dimension " << syslen << " x " << syslen << std::endl;
  result_file << "version=" << FDASettings::compat_fm_version << std::endl;
  result_file << "groupname=" << fda_settings.groupname << std::endl;
  result_file << "writefreq=1" << std::endl;
  // Reserve place for up to 8 digits for nsteps
  result_file << "nsteps=" << std::setfill('0') << std::setw(8) << nsteps << std::endl;
  result_file << "sysanr=" << syslen << std::endl;
  result_file << "fmdim=" << syslen << std::endl;
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
  result_file << std::scientific << std::setprecision(6);

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
  result_file << std::scientific << std::setprecision(6);

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
