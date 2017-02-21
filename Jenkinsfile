#!groovy

pipeline {

  agent {
    label "docker-nodes"
  }
  
  stages {

    stage("Build") {
      parallel {
        "gcc-4.9" : {
          agent {
            docker "bernddoser/docker-devel-cpp:ubuntu-gcc-4.9"
          }
          steps {
            sh 'mkdir -p build-gcc-4.9'
            sh 'cd build-gcc-4.9 && cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DGMX_BUILD_MDRUN_ONLY=OFF -DGMX_BUILD_FDA=ON -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_fda -DGMX_SIMD=NONE -DGMX_BUILD_UNITTESTS=ON -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF ..'
            sh 'cd build-gcc-4.9 && make'
          }
        }
        "clang-3.9" : {
          agent {
            docker "bernddoser/docker-devel-cpp:ubuntu-clang-3.9"
          }
          steps {
            sh 'mkdir -p build-clang-3.9'
            sh 'cd build-clang-3.9 && cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DGMX_BUILD_MDRUN_ONLY=OFF -DGMX_BUILD_FDA=ON -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_fda -DGMX_SIMD=NONE -DGMX_BUILD_UNITTESTS=ON -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF ..'
            sh 'cd build-clang-3.9 && make'
          }
        }
      }
    }
    
    stage("Test") {
      steps {
        sh 'cd build && make check'
      }
    }
  }

//  post {
//    always {
//      deleteDir()
//    }
//    failure {
//      mail to:"bernd.doser@h-its.org", subject:"FAILURE: ${currentBuild.fullDisplayName}", body: "Boo, we failed."
//    }
//  }

}
