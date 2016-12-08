/*
 * PairwiseForces.cpp
 *
 *  Created on: Dec 8, 2016
 *      Author: Bernd Doser, HITS gGmbH
 */

#include <gtest/gtest.h>
#include "gromacs/fda/PairwiseForces.h"
#include "testutils/integrationtests.h"

namespace fda {

//! Test fixture for PairwiseForces
class PairwiseForcesTest : public ::gmx::test::IntegrationTestFixture
{};

TEST_F(PairwiseForcesTest, DefaultConstructor)
{
  PairwiseForces<real> pf;
  EXPECT_TRUE(pf.all_pairwise_forces.empty());
}

TEST_F(PairwiseForcesTest, ReadFile)
{
  std::string data_path = std::string(fileManager_.getInputDataDirectory()) + "/data";

  PairwiseForces<Force<real>> pf(data_path + "/test1.pfa");

  EXPECT_EQ(0, pf.all_pairwise_forces.size());
}

} // namespace fda
