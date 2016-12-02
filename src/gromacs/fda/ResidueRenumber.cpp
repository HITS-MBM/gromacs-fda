/*
 * ResiduesRenumber.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "ResidueRenumber.h"

namespace fda {

std::ostream& operator << (std::ostream& os, ResiduesRenumber r)
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

std::istream& operator >> (std::istream& is, ResiduesRenumber& r)
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
