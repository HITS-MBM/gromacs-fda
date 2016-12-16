#!groovy

#!groovy

stage('Checkout') {
  node('master') {
    git url: "https://github.com/HITS-MBM/gromacs-fda.git", branch: 'release-2016-fda'
  }
}

stage('Build') {
  node('master') {
    sh 'cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DGMX_BUILD_MDRUN_ONLY=OFF -DGMX_BUILD_FDA=ON -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_fda -DGMX_SIMD=NONE -DGMX_BUILD_UNITTESTS=ON -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF .'
    sh 'make'
  }
}

stage('Test') {
  node('master') {
    sh 'make check'
  }
}

stage('Promotion') {
  try {
    timeout(time:5, unit:'MINUTES') {
      input message:'Approve deployment?', submitter: 'admin'
    }
  } catch (Exception ex) {
    annotate('currentStage.status', 'UNSTABLE');
  }
}

stage('Deployment') {
  node('master') {
    sh 'make install'
  }
}

