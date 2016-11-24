/*
 * PureInteractionType.h
 *
 *  Created on: Nov 24, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_PUREINTERACTIONTYPE_H_
#define SRC_GROMACS_FDA_PUREINTERACTIONTYPE_H_

#include "InteractionType.h"

namespace fda {

/// Used for array indexes!  Don't change the numbers!
enum class PureInteractionType : int
{
  BOND      = 0,
  ANGLE     = 1,
  DIHEDRAL  = 2,
  POLAR     = 3,
  COULOMB   = 4,
  LJ        = 5,
  NB14      = 6,
  NUMBER    = 7
};

PureInteractionType to_pure(InteractionType i)
{
  switch(i) {
	case InteractionType::BOND:
	  return PureInteractionType::BOND;
	case InteractionType::ANGLE:
	  return PureInteractionType::ANGLE;
	case InteractionType::DIHEDRAL:
	  return PureInteractionType::DIHEDRAL;
	case InteractionType::POLAR:
	  return PureInteractionType::POLAR;
	case InteractionType::COULOMB:
	  return PureInteractionType::COULOMB;
	case InteractionType::LJ:
	  return PureInteractionType::LJ;
	case InteractionType::NB14:
	  return PureInteractionType::NB14;
	default:
	  throw std::runtime_error("Is not a pure interaction");
  }
}

InteractionType from_pure(PureInteractionType i)
{
  switch(i) {
	case PureInteractionType::BOND:
	  return InteractionType::BOND;
	case PureInteractionType::ANGLE:
	  return InteractionType::ANGLE;
	case PureInteractionType::DIHEDRAL:
	  return InteractionType::DIHEDRAL;
	case PureInteractionType::POLAR:
	  return InteractionType::POLAR;
	case PureInteractionType::COULOMB:
	  return InteractionType::COULOMB;
	case PureInteractionType::LJ:
	  return InteractionType::LJ;
	case PureInteractionType::NB14:
	  return InteractionType::NB14;
	default:
	  throw std::runtime_error("Is not a pure interaction");
  }
}

} // namespace fda

#endif /* SRC_GROMACS_FDA_PUREINTERACTIONTYPE_H_ */
