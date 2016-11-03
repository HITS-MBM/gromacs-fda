/*
 * ResultType.h
 *
 *  Created on: Nov 3, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_RESULTTYPE_H_
#define SRC_GROMACS_FDA_RESULTTYPE_H_

#include <cstdint>
#include <iostream>
#include <string>

namespace fda {

enum class ResultType : std::int8_t
{
  NO,                       // no storing (default)
  PAIRWISE_FORCES_VECTOR,
  PAIRWISE_FORCES_SCALAR,
  PUNCTUAL_STRESS,
  VIRIAL_STRESS,
  VIRIAL_STRESS_VON_MISES,
  COMPAT_BIN,               // DEPRICATED! compatibility mode (signed scalars) in binary
  COMPAT_ASCII,             // DEPRICATED! compatibility mode (signed scalars) in ascii
  INVALID                   // unknown type from string conversion
};

std::string to_string(ResultType const& r)
{
  switch(r) {
    case ResultType::NO:
      return "none";
    case ResultType::PAIRWISE_FORCES_VECTOR:
      return "pairwise_forces_vector";
    case ResultType::PAIRWISE_FORCES_SCALAR:
      return "pairwise_forces_scalar";
    case ResultType::PUNCTUAL_STRESS:
      return "punctual_stress";
    case ResultType::VIRIAL_STRESS:
      return "virial_stress";
    case ResultType::VIRIAL_STRESS_VON_MISES:
      return "virial_stress_von_mises";
    case ResultType::COMPAT_BIN:
      return "compat_bin";
    case ResultType::COMPAT_ASCII:
      return "compat_ascii";
    default:
      return "invalid";
  }
}

/// Input operator for ResultType
std::istream& operator>>(std::istream& is, ResultType& r)
{
  std::string s;
  is >> s;
  if (s == "none")
    r = ResultType::NO;
  else if (s == "pairwise_forces_vector")
    r = ResultType::PAIRWISE_FORCES_VECTOR;
  else if (s == "pairwise_forces_scalar")
    r = ResultType::PAIRWISE_FORCES_SCALAR;
  else if (s == "punctual_stress")
    r = ResultType::PUNCTUAL_STRESS;
  else if (s == "virial_stress")
    r = ResultType::VIRIAL_STRESS;
  else if (s == "virial_stress_von_mises")
    r = ResultType::VIRIAL_STRESS_VON_MISES;
  else if (s == "compat_bin")
    r = ResultType::COMPAT_BIN;
  else if (s == "compat_ascii")
    r = ResultType::COMPAT_ASCII;
  else
    r = ResultType::INVALID;
  return is;
}

} // namespace fda

#endif /* SRC_GROMACS_FDA_RESULTTYPE_H_ */
