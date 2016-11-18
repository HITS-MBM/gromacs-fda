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

namespace fda {

struct InteractionType
{
  static const int bond      = 1 << 0;
  static const int angle     = 1 << 1;
  static const int dihedral  = 1 << 2;
  static const int polar     = 1 << 3;
  static const int coulomb   = 1 << 4;
  static const int lj        = 1 << 5;
  static const int nb14      = 1 << 6;
  static const int num       =      7;
  static const int bonded    = bond + angle + dihedral;
  static const int nonbonded = coulomb + lj + nb14;
  static const int all       = bonded + nonbonded + polar;
};

#if 0
enum class InteractionType : int
{
  BOND,
  ANGLE,
  DIHEDRAL,
  POLAR,
  COULOMB,
  LJ,
  NB14,
  NUM
};

/// Output stream for InteractionType
std::ostream& operator>>(std::ostream& os, InteractionType r)
{
  switch(r) {
    case InteractionType::BOND:
      return os << "bond";
    case InteractionType::ANGLE:
      return os << "angle";
    case InteractionType::DIHEDRAL:
      return os << "dihedral";
    case InteractionType::POLAR:
      return os << "polar";
    case InteractionType::COULOMB:
      return os << "coulomb";
    case InteractionType::LJ:
      return os << "lj";
    case InteractionType::NB14:
      return os << "nb14";
    default:
      return os << "invalid";
  }
}

/// Input stream for InteractionType
std::istream& operator>>(std::istream& is, InteractionType& r)
{
  std::string s;
  is >> s;
  if (s == "bond")
	r = InteractionType::BOND;
  else if (s == "angle")
	r = InteractionType::ANGLE;
  else if (s == "dihedral")
	r = InteractionType::DIHEDRAL;
  else if (s == "polar")
	r = InteractionType::POLAR;
  else if (s == "coulomb")
	r = InteractionType::COULOMB;
  else if (s == "lj")
	r = InteractionType::LJ;
  else if (s == "nb14")
	r = InteractionType::NB14;
  else
	r = InteractionType::NUM;
  return is;
}

} // namespace fda

#endif
#endif /* SRC_GROMACS_FDA_INTERACTIONTYPE_H_ */
