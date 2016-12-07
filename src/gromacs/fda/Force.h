/*
 * Force.h
 *
 *  Created on: Dec 7, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_FORCE_H_
#define SRC_GROMACS_FDA_FORCE_H_

#include <iostream>
#include "gromacs/utility/real.h"
#include "InteractionType.h"

namespace fda {

template <typename T>
struct Force
{
  Force(T force = 0.0, InteractionType type = 0)
   : force(force), type(type)
  {}

  /// Scalar force
  T force;

  /// Interaction type
  InteractionType type;
};

/// Output stream for ResultType
template <typename T>
std::ostream& operator << (std::ostream& os, Force<T> const& f)
{
  return os << f.force << f.type;
}

/// Input stream for ResultType
template <typename T>
std::istream& operator >> (std::istream& is, Force<T> & f)
{
  return is >> f.force >> f.type;
}

} // namespace fda

#endif /* SRC_GROMACS_FDA_FORCE_H_ */
