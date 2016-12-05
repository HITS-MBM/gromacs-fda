/*
 * InteractionType.h
 *
 *  Created on: Nov 18, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_INTERACTIONTYPE_H_
#define SRC_GROMACS_FDA_INTERACTIONTYPE_H_

#include <cstdint>
#include <iostream>
#include <string>
#include <type_traits>

namespace fda {

enum class InteractionType : int
{
  NONE      =      0,
  BOND      = 1 << 0,
  ANGLE     = 1 << 1,
  DIHEDRAL  = 1 << 2,
  POLAR     = 1 << 3,
  COULOMB   = 1 << 4,
  LJ        = 1 << 5,
  NB14      = 1 << 6,
  BONDED    = BOND + ANGLE + DIHEDRAL,
  NONBONDED = COULOMB + LJ + NB14,
  ALL       = BONDED + NONBONDED + POLAR
};

using T = std::underlying_type<InteractionType>::type;

constexpr T to_index(InteractionType e)
{
  return static_cast<T>(e);
}

inline InteractionType operator & (InteractionType lhs, InteractionType rhs)
{
  return (InteractionType)(static_cast<T>(lhs) & static_cast<T>(rhs));
}

inline InteractionType operator | (InteractionType lhs, InteractionType rhs)
{
  return (InteractionType)(static_cast<T>(lhs) | static_cast<T>(rhs));
}

inline InteractionType operator + (InteractionType lhs, InteractionType rhs)
{
  return (InteractionType)(static_cast<T>(lhs) + static_cast<T>(rhs));
}

inline bool operator ! (InteractionType i)
{
  return !static_cast<T>(i);
}

/// Output stream for InteractionType
std::ostream& operator << (std::ostream& os, InteractionType r);

/// Input stream for InteractionType
std::istream& operator >> (std::istream& is, InteractionType& r);

} // namespace fda

#endif /* SRC_GROMACS_FDA_INTERACTIONTYPE_H_ */
