set(testname "FDAIntegrationTests")
set(exename "fda-integration-tests")

gmx_add_gtest_executable(
    ${exename}
    # files with code for tests
    FDATest.cpp
    # pseudo-library for code for mdrun
    $<TARGET_OBJECTS:mdrun_objlib>
)
    
gmx_register_gtest_test(
    ${testname}
    ${exename}
    INTEGRATION_TEST
)
