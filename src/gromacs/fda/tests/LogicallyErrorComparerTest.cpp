/*
 * LogicallyErrorComparerTest.cpp
 *
 *  Created on: Jun 26, 2017
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <gtest/gtest.h>
#include "testutils/LogicallyErrorComparer.h"

namespace fda
{

TEST(LogicallyErrorComparerTest, Test1)
{
    LogicallyEqualComparer<true, true> comparer(1e2);

    EXPECT_TRUE(comparer(1.516308e-06, -1.140742e-06));
}

} // namespace fda
