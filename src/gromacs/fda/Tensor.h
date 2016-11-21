/*
 * Tensor.h
 *
 *  Created on: Nov 21, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_TENSOR_H_
#define SRC_GROMACS_FDA_TENSOR_H_

#include "gromacs/math/vectypes.h"

namespace fda {

/// As currently no RAII algebra types are availbale in GROMACS
/// a wrapper class is used around tensor.
class Tensor
{
public:

	Tensor()
	{
		t[0][0] = 0.0; t[0][1] = 0.0; t[0][2] = 0.0;
		t[1][0] = 0.0; t[1][1] = 0.0; t[1][2] = 0.0;
		t[2][0] = 0.0; t[2][1] = 0.0; t[2][2] = 0.0;
	}

private:

	tensor t;
};

} // namespace fda

#endif /* SRC_GROMACS_FDA_TENSOR_H_ */
