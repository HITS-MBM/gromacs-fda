/*
 * FDABase.cpp
 *
 *  Created on: Nov 22, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <cmath>
#include "FDABase.h"
#include "gromacs/math/vec.h"
#include "Utilities.h"

namespace fda {

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
		case ResultType::INVALID:
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
		  break;
		case ResultType::COMPAT_ASCII:
		  write_frame_atoms_summed_compat(x, true);
		  break;
		case ResultType::INVALID:
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
        if (bVector)
          result_file << i << " " << j << " " << id.force[XX] << " " << id.force[YY] << " " << id.force[ZZ] << " " << type << std::endl;
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
void FDABase<Base>::sum_total_forces(rvec *x)
{
  for (auto& v : total_forces) v = 0.0;

  for (auto const& v1 : distributed_forces.summed) {
    int i = v1.first;
    for (auto const& v2 : v1.second) {
      int j = v2.first;
      real scalar_force;
      switch (v2s) {
        case Vector2Scalar::NORM:
          scalar_force = norm(v2.second.get_rvec());
          break;
        case Vector2Scalar::PROJECTION:
          scalar_force = vector2unsignedscalar(v2.second.get_rvec(), i, j, x);
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

template <>
void FDABase<Atom>::write_frame_atoms_compat(DistributedForces const& forces, FILE *f, int *framenr, int interactions_count, int *fmatoms, int *index, real *force, char *interaction, gmx_bool ascii) {
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

template <>
void FDABase<Atom>::write_frame_atoms_summed_compat(rvec *x, bool ascii)
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

/// template instantiation
template class FDABase<Atom>;
template class FDABase<Residue>;

} // namespace fda
