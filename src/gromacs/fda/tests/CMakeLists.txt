set(testname "FDATests")
set(exename "fda-tests")

gmx_add_gtest_executable(
    ${exename}
    PairwiseForcesTest.cpp
)

gmx_register_unit_test(
    ${testname}
    ${exename}
)
