/*
 * Vector.h
 *
 *  Created on: Nov 21, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_VECTOR_H_
#define SRC_GROMACS_FDA_VECTOR_H_

#include "gromacs/math/vectypes.h"

namespace fda {

/// As currently no RAII algebra types are available in GROMACS
/// a wrapper class is used around rvec.
class Vector
{
public:

	Vector()
	{ v[0] = 0.0; v[1] = 0.0; v[2] = 0.0; }

	void operator+=(rvec other)
	{
	    v[0] += other[0]; v[1] += other[1]; v[2] += other[2];
	}

	real* get_rvec() { return v; }
	real const* get_rvec() const { return v; }

private:

	rvec v;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_VECTOR_H_ */
