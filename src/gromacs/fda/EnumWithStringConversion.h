/*
 * EnumWithStringConversion.h
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_ENUMWITHSTRINGCONVERSION_H_
#define SRC_GROMACS_FDA_ENUMWITHSTRINGCONVERSION_H_

#include <cstdint>
#include <iostream>
#include <string>

#define ENUM_WITH_STRING_CONVERSION(name, type, ...) \
  enum class name : type \
  { \
    __VA_ARGS__ \
  };

/// Output stream for ForceType
std::ostream& operator>>(std::ostream& os, ForceType r)
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

#endif /* SRC_GROMACS_FDA_ENUMWITHSTRINGCONVERSION_H_ */
