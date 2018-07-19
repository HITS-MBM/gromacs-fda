/*
 * PairwiseForces.cpp
 *
 *  Created on: Dec 8, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <iostream>
#include <gtest/gtest.h>
#include "gromacs/fda/PairwiseForces.h"
#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace fda
{

//! Test fixture for PairwiseForces
class PairwiseForcesTest : public gmx::test::CommandLineTestBase
{};

TEST_F(PairwiseForcesTest, ReadFile1)
{
    std::string data_path = std::string(fileManager().getInputDataDirectory()) + "/data";

    PairwiseForces<Force<real>> pf(data_path + "/test1.pfa");
    auto&& pf_all = pf.get_all_pairwise_forces(true);
	std::cout << pf << std::endl;

    double abs_error_float = 1e-6;

    EXPECT_EQ(1, pf_all.size());
    EXPECT_EQ(3, pf_all[0].size());
    EXPECT_EQ(0, pf_all[0][0].i);
    EXPECT_EQ(14, pf_all[0][0].j);
    EXPECT_NEAR(3.780676, pf_all[0][0].force.force, 3.780676 * abs_error_float);
    EXPECT_EQ(16, pf_all[0][0].force.type);
    EXPECT_EQ(11, pf_all[0][2].i);
    EXPECT_EQ(14, pf_all[0][2].j);
    EXPECT_NEAR(8.904991e+02, pf_all[0][2].force.force, 8.904991e+02 * abs_error_float);
    EXPECT_EQ(64, pf_all[0][2].force.type);
}

TEST_F(PairwiseForcesTest, ReadFile2)
{
    std::string data_path = std::string(fileManager().getInputDataDirectory()) + "/data";

    PairwiseForces<Force<real>> pf(data_path + "/test2.pfa");
    auto pf_all = pf.get_all_pairwise_forces(true);

    EXPECT_EQ(11, pf_all.size());
    EXPECT_EQ(142, pf_all[0].size());
    EXPECT_EQ(0, pf_all[0][0].i);
    EXPECT_EQ(12, pf_all[0][0].j);
}

TEST_F(PairwiseForcesTest, ReadFile3)
{
    std::string data_path = std::string(fileManager().getInputDataDirectory()) + "/data";

    PairwiseForces<Force<Vector>> pf(data_path + "/test3.pfa");
    auto pf_all = pf.get_all_pairwise_forces(true);

    EXPECT_EQ(11, pf_all.size());
    EXPECT_EQ(0, pf_all[8][5].i);
    EXPECT_EQ(17, pf_all[8][5].j);
}

} // namespace fda
