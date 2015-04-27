set(testname "FDATests")
set(exename "fda-test")

gmx_build_unit_test(
    ${testname}
    ${exename}
    # files with code for tests
    FDATest.cpp
    # pseudo-library for code for mdrun
    $<TARGET_OBJECTS:mdrun_objlib>
)
    
gmx_register_integration_test(
    ${testname}
    ${exename}
)
