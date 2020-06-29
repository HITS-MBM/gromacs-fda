/*
 * Utilities.cpp
 *
 *  Created on: Nov 22, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "gromacs/math/vec.h"
#include "Utilities.h"

namespace fda {

real vector2signedscalar(const rvec v, const rvec xi, const rvec xj, const matrix box, Vector2Scalar v2s)
{
    rvec r;
    rvec_sub(xj, xi, r);

    // periodic boundary correction of bonds crossing the border
    if      (std::abs(r[0] + box[0][0]) < std::abs(r[0])) r[0] += box[0][0];
    else if (std::abs(r[0] - box[0][0]) < std::abs(r[0])) r[0] -= box[0][0];
    if      (std::abs(r[1] + box[1][1]) < std::abs(r[1])) r[1] += box[1][1];
    else if (std::abs(r[1] - box[1][1]) < std::abs(r[1])) r[1] -= box[1][1];
    if      (std::abs(r[2] + box[2][2]) < std::abs(r[2])) r[2] += box[2][2];
    else if (std::abs(r[2] - box[2][2]) < std::abs(r[2])) r[2] -= box[2][2];

    real c = cos_angle(r, v);
    switch (v2s) {
        case Vector2Scalar::NORM:
            // based on the cos of the angle between vectors:
            // if it's positive, the vector v is oriented within -90 to 90 degress with respect to r - considered to be same direction = attractive = negative;
            // if it's negative, the vector v is oriented within 90 to 270 degrees with respect to r - considered to be in opposite direction = repulsvive = positive;
            // if it's zero, the vector v is oriented at exactly 90 or 270 degrees with respect to r - the force can be considered neither attractive nor repulsvive, so it's set to zero
            if (c > 0) {
                return -norm(v);
            } else if (c < 0) {
                return norm(v);
            } else {
                return 0.0;
            }
            break;
        case Vector2Scalar::PROJECTION:
            // it's minus c to go along the norm based calculation above: positive cos = in the same direction = attractive = negative
            return -c * norm(v);
            break;
    }
    return GMX_FLOAT_MAX;
}

real vector2unsignedscalar(const rvec v, int i, int j, gmx::PaddedHostVector<gmx::RVec> const& x)
{
    rvec r;
    rvec_sub(x[j], x[i], r);
    real p = cos_angle(r, v) * norm(v);
    return std::abs(p);
}

} // namespace fda
