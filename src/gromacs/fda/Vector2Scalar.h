/*
 * Vector2Scalar.h
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_VECTOR2SCALAR_H_
#define SRC_GROMACS_FDA_VECTOR2SCALAR_H_

#include <cstdint>
#include <iostream>
#include <string>

namespace fda {

enum class Vector2Scalar : std::int8_t
{
  NORM,
  PROJECTION,
  INVALID
};

/// Output stream for Vector2Scalar
std::ostream& operator>>(std::ostream& os, Vector2Scalar r)
{
  switch(r) {
    case Vector2Scalar::NORM:
      return os << "norm";
    case Vector2Scalar::PROJECTION:
      return os << "projection";
    default:
      return os << "invalid";
  }
}

/// Input stream for Vector2Scalar
std::istream& operator>>(std::istream& is, Vector2Scalar& r)
{
  std::string s;
  is >> s;
  if (s == "norm")
	r = Vector2Scalar::NORM;
  else if (s == "projection")
	r = Vector2Scalar::PROJECTION;
  else
	r = Vector2Scalar::INVALID;
  return is;
}

} // namespace fda

#endif /* SRC_GROMACS_FDA_VECTOR2SCALAR_H_ */
