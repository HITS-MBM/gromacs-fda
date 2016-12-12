/*
 * DistributedForces.cpp
 *
 *  Created on: Nov 3, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "DistributedForces.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "Utilities.h"

namespace fda {

DistributedForces::DistributedForces(int syslen, FDASettings const& fda_settings)
 : syslen(syslen),
   scalar(syslen),
   summed(syslen),
   detailed(syslen),
   fda_settings(fda_settings)
{}

void DistributedForces::add_scalar(int i, int j, real force, InteractionType type)
{
  if (i > j) throw std::runtime_error("Only upper triangle allowed (i < j).");

  auto& scalar_i = scalar[i];

  int pos_j;
  if (lookup.count(j)) {
    pos_j = lookup[j];
  } else {
	pos_j = scalar_i.size();
	scalar_i.resize(scalar_i.size() + 1);
	lookup[j] = pos_j;
	reverse_lookup[pos_j] = j;
  }

  scalar_i[pos_j].force += force;
  scalar_i[pos_j].type |= type;
}

void DistributedForces::add_summed(int i, int j, Vector const& force, InteractionType type)
{
  if (i > j) throw std::runtime_error("Only upper triangle allowed (i < j).");

  auto& summed_i = summed[i];

  int pos_j;
  if (lookup.count(j)) {
    pos_j = lookup[j];
  } else {
	pos_j = summed_i.size();
	summed_i.resize(summed_i.size() + 1);
	lookup[j] = pos_j;
	reverse_lookup[pos_j] = j;
  }

  summed_i[pos_j].force += force;
  summed_i[pos_j].type |= type;
}

void DistributedForces::add_detailed(int i, int j, Vector const& force, PureInteractionType type)
{
  if (i > j) throw std::runtime_error("Only upper triangle allowed (i < j).");

  auto& detailed_i = detailed[i];

  int pos_j;
  if (lookup.count(j)) {
    pos_j = lookup[j];
  } else {
	pos_j = detailed_i.size();
	detailed_i.resize(detailed_i.size() + 1);
	lookup[j] = pos_j;
	reverse_lookup[pos_j] = j;
  }

  detailed_i[pos_j].force[to_index(type)] += force;
  ++detailed_i[pos_j].number[to_index(type)];
}

void DistributedForces::write_detailed_vector(std::ostream& os) const
{
  for (size_t i = 0; i != detailed.size(); ++i) {
	auto detailed_i = detailed[i];
    for (size_t j = 0; j != detailed_i.size(); ++j) {
      auto detailed_j = detailed_i[reverse_lookup.at(j)];
      for (int type = 0; type != static_cast<int>(PureInteractionType::NUMBER); ++type) {
    	Vector force = detailed_j.force[type];
    	if (detailed_j.number[type] == 0) continue;
    	InteractionType interaction_type = from_pure(static_cast<PureInteractionType>(type));
        os << i << " " << j << " " << force[XX] << " " << force[YY] << " " << force[ZZ] << " " << interaction_type << std::endl;
      }
    }
  }
}

void DistributedForces::write_detailed_scalar(std::ostream& os, rvec *x) const
{
  for (size_t i = 0; i != detailed.size(); ++i) {
	auto detailed_i = detailed[i];
    for (size_t j = 0; j < detailed_i.size(); ++j) {
      auto detailed_j = detailed_i[reverse_lookup.at(j)];
      for (int type = 0; type != static_cast<int>(PureInteractionType::NUMBER); ++type) {
    	Vector force = detailed_j.force[type];
    	if (detailed_j.number[type] == 0) continue;
    	InteractionType interaction_type = from_pure(static_cast<PureInteractionType>(type));
        os << i << " " << j << " " << vector2signedscalar(force.get_pointer(), x[i], x[j], fda_settings.v2s) << " " << interaction_type << std::endl;
      }
    }
  }
}

void DistributedForces::write_summed_vector(std::ostream& os) const
{
  for (size_t i = 0; i != summed.size(); ++i) {
	auto summed_i = summed[i];
	for (size_t j = 0; j != summed_i.size(); ++j) {
	  auto summed_j = summed_i[reverse_lookup.at(j)];
      Vector force = summed_j.force;
      os << i << " " << j << " " << force[XX] << " " << force[YY] << " " << force[ZZ] << " " << summed_j.type << std::endl;
	}
  }
}

void DistributedForces::write_summed_scalar(std::ostream& os, rvec *x) const
{
  for (size_t i = 0; i != summed.size(); ++i) {
	auto summed_i = summed[i];
	for (size_t j = 0; j != summed_i.size(); ++j) {
	  auto summed_j = summed_i[reverse_lookup.at(j)];
      Vector force = summed_j.force;
      os << i << " " << j << " " << vector2signedscalar(force.get_pointer(), x[i], x[j], fda_settings.v2s) << " " << summed_j.type << std::endl;
	}
  }
}

void DistributedForces::write_scalar(std::ostream& os) const
{
  for (size_t i = 0; i != scalar.size(); ++i) {
	auto scalar_i = scalar[i];
	for (size_t pos_j = 0; pos_j != scalar_i.size(); ++pos_j) {
	  size_t j = reverse_lookup.at(pos_j);
	  auto scalar_j = scalar_i[j];
	  os << i << " " << j << " " << scalar_j.force << " " << scalar_j.type << std::endl;
	}
  }
}

void DistributedForces::write_total_forces(std::ostream& os, rvec *x) const
{
  std::vector<real> total_forces(syslen, 0.0);
  for (size_t i = 0; i != summed.size(); ++i) {
	auto summed_i = summed[i];
	for (size_t pos_j = 0; pos_j != summed_i.size(); ++pos_j) {
      size_t j = reverse_lookup.at(pos_j);
	  auto summed_j = summed_i[j];
      real scalar_force;
      switch (fda_settings.v2s) {
        case Vector2Scalar::NORM:
          scalar_force = norm(summed_j.force.get_pointer());
          break;
        case Vector2Scalar::PROJECTION:
          scalar_force = vector2unsignedscalar(summed_j.force.get_pointer(), i, j, x);
          break;
        default:
      	  gmx_fatal(FARGS, "Unknown option for Vector2Scalar.\n");
          break;
      }
      total_forces[i] += scalar_force;
      total_forces[j] += scalar_force;
    }
  }

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
      os << total_forces[i];
      first_on_line = false;
    } else {
      os << " " << total_forces[i];
    }
  }
  os << std::endl;
}

void DistributedForces::scalar_real_divide(real divisor)
{
  real inv = 1 / divisor;
  for (auto& i : scalar)
    for (auto& j : i) j.force *= inv;
}

void DistributedForces::summed_merge_to_scalar(const rvec *x)
{
  for (size_t i = 0; i < summed.size(); ++i) {
    auto summed_i = summed[i];
    auto scalar_i = scalar[i];
    for (size_t j = 0; j < summed_i.size(); ++j) {
      auto summed_j = summed_i[reverse_lookup[j]];
      auto scalar_j = scalar_i[reverse_lookup[j]];
      scalar_j.force += vector2signedscalar(summed_j.force.get_pointer(), x[i], x[j], fda_settings.v2s);
      scalar_j.type |= summed_j.type;
    }
  }
}

} // namespace fda
