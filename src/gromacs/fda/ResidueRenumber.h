/*
 * ResiduesRenumber.h
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_RESIDUESRENUMBER_H_
#define SRC_GROMACS_FDA_RESIDUESRENUMBER_H_

#include <cstdint>
#include <iostream>
#include <string>

namespace fda {

enum class ResiduesRenumber : std::int8_t
{
  AUTO,
  DO,
  DONT,
  INVALID
};

/// Output stream for ResiduesRenumber
std::ostream& operator>>(std::ostream& os, ResiduesRenumber r)
{
  switch(r) {
    case ResiduesRenumber::AUTO:
      return os << "auto";
    case ResiduesRenumber::DO:
      return os << "do";
    case ResiduesRenumber::DONT:
      return os << "dont";
    default:
      return os << "invalid";
  }
}

/// Input stream for ResiduesRenumber
std::istream& operator>>(std::istream& is, ResiduesRenumber& r)
{
  std::string s;
  is >> s;
  if (s == "auto")
	r = ResiduesRenumber::AUTO;
  else if (s == "do")
	r = ResiduesRenumber::DO;
  else if (s == "dont")
	r = ResiduesRenumber::DONT;
  else
	r = ResiduesRenumber::INVALID;
  return is;
}

} // namespace fda

#endif /* SRC_GROMACS_FDA_RESIDUESRENUMBER_H_ */
