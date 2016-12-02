/*
 * OnePair.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "OnePair.h"

namespace fda {

std::ostream& operator << (std::ostream& os, OnePair r)
{
  switch(r) {
    case OnePair::DETAILED:
      return os << "detailed";
    case OnePair::SUMMED:
      return os << "summed";
    default:
      return os << "invalid";
  }
}

std::istream& operator >> (std::istream& is, OnePair& r)
{
  std::string s;
  is >> s;
  if (s == "detailed")
	r = OnePair::DETAILED;
  else if (s == "summed")
	r = OnePair::SUMMED;
  else
	r = OnePair::INVALID;
  return is;
}

} // namespace fda
