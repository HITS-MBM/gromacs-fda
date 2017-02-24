#!groovy

pipeline {
  agent {
    label "docker-nodes"
  }
  stages {
    stage('Build') {
      steps {
        parallel(
          "gcc-4.9": {
            agent {
              docker 'bernddoser/docker-devel-cpp:ubuntu-16.04-gcc-4.9'
              label "docker-nodes"
            }
            steps {
              sh 'mkdir -p build'
              sh 'cd build && cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DGMX_BUILD_MDRUN_ONLY=OFF -DGMX_BUILD_FDA=ON -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_fda -DGMX_SIMD=NONE -DGMX_BUILD_UNITTESTS=ON -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF ..'
              sh 'cd build && make'
            }
          },
          "clang-3.9": {
            agent {
              docker 'bernddoser/docker-devel-cpp:ubuntu-16.04-clang-3.9'
              label "docker-nodes"
            }
            steps {
              sh 'mkdir -p build'
              sh 'cd build && cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++ -DGMX_BUILD_MDRUN_ONLY=OFF -DGMX_BUILD_FDA=ON -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_fda -DGMX_SIMD=NONE -DGMX_BUILD_UNITTESTS=ON -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF ..'
              sh 'cd build && make'
            }
          }
        )
      }
    }
    stage('Test') {
      steps {
        parallel(
          "gcc-4.9": {
            sh 'cd build && make check'
            
          },
          "clang-3.9": {
            sh 'cd build && make check'
            
          }
        )
      }
    }
    stage('Dokumentation') {
      steps {
        sh 'cd build && make doc'
      }
    }
    stage('Deploy') {
      steps {
        archiveArtifacts '*.deb'
      }
    }
  }
  post {
    always {
      deleteDir()
    }
    failure {
      mail to:"bernd.doser@h-its.org", subject:"FAILURE: ${currentBuild.fullDisplayName}", body: "Boo, we failed."
    }
  }
}
