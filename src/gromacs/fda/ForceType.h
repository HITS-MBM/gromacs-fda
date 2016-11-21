/*
 * ForceType.h
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_FORCETYPE_H_
#define SRC_GROMACS_FDA_FORCETYPE_H_

#include <cstdint>
#include <iostream>
#include <string>

namespace fda {

enum class ForceType : std::int8_t
{
  ATOMS,
  RESIDUES,
  INVALID
};

/// Output stream for ForceType
std::ostream& operator<<(std::ostream& os, ForceType r)
{
  switch(r) {
    case ForceType::ATOMS:
      return os << "atoms";
    case ForceType::RESIDUES:
      return os << "residues";
    default:
      return os << "invalid";
  }
}

/// Input stream for ForceType
std::istream& operator>>(std::istream& is, ForceType& r)
{
  std::string s;
  is >> s;
  if (s == "atoms")
	r = ForceType::ATOMS;
  else if (s == "residues")
	r = ForceType::RESIDUES;
  else
	r = ForceType::INVALID;
  return is;
}

} // namespace fda

#endif /* SRC_GROMACS_FDA_FORCETYPE_H_ */
