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
  INVALID                   // unknown type by string conversion
};

/// Conversion to string
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

/// Conversion from string
ResultType from_string(std::string const& s)
{
  if (s == "none")
	return ResultType::NO;
  else if (s == "pairwise_forces_vector")
	return ResultType::PAIRWISE_FORCES_VECTOR;
  else if (s == "pairwise_forces_scalar")
	return ResultType::PAIRWISE_FORCES_SCALAR;
  else if (s == "punctual_stress")
	return ResultType::PUNCTUAL_STRESS;
  else if (s == "virial_stress")
	return ResultType::VIRIAL_STRESS;
  else if (s == "virial_stress_von_mises")
	return ResultType::VIRIAL_STRESS_VON_MISES;
  else if (s == "compat_bin")
	return ResultType::COMPAT_BIN;
  else if (s == "compat_ascii")
	return ResultType::COMPAT_ASCII;
  else
	return ResultType::INVALID;
}

/// Input operator for ResultType
std::istream& operator>>(std::istream& is, ResultType& r)
{
  std::string s;
  is >> s;
  r = from_string(s);
  return is;
}

} // namespace fda

#endif /* SRC_GROMACS_FDA_RESULTTYPE_H_ */
