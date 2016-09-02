node {
    stage 'Checkout'
    git url: "https://github.com/HITS-MBM/gromacs-fda.git", branch: 'fda'
    
    stage 'Build'
    sh 'cmake -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DGMX_BUILD_MDRUN_ONLY=OFF -DGMX_BUILD_FDA=ON -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_pf2 -DGMX_SIMD=NONE -DGMX_BUILD_UNITTESTS=ON -DGMX_BUILD_OWN_FFTW=ON .'
    sh 'make'
    
    stage 'Test'
    sh 'make check'
}
