/*
 * InteractionType.cpp
 *
 *  Created on: Nov 18, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <stdexcept>
#include "InteractionType.h"

namespace fda {

std::ostream& operator << (std::ostream& os, InteractionType r)
{
  switch(r) {
  case InteractionType::NONE:
    return os << "none";
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
    case InteractionType::BONDED:
      return os << "bonded";
    case InteractionType::NONBONDED:
      return os << "nonbonded";
    case InteractionType::ALL:
      return os << "all";
    default:
      return os << "invalid";
  }
}

std::istream& operator >> (std::istream& is, InteractionType& r)
{
  std::string s;
  is >> s;
  if (s == "none")
	r = InteractionType::NONE;
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
  else if (s == "bonded")
	r = InteractionType::BONDED;
  else if (s == "nonbonded")
	r = InteractionType::NONBONDED;
  else if (s == "all")
	r = InteractionType::ALL;
  else
	throw std::runtime_error("InteractionType " + s + "unknown");
  return is;
}

} // namespace fda
