set(testname "FDATests")
set(exename "fda-test")

gmx_add_gtest_executable(
    ${exename}
    # files with code for tests
    FDATest.cpp
    PairwiseForces.cpp
    # pseudo-library for code for mdrun
    $<TARGET_OBJECTS:mdrun_objlib>
)
    
gmx_register_integration_test(
    ${testname}
    ${exename}
)
