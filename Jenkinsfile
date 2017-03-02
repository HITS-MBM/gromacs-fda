#!groovy

pipeline {
  agent {
    label 'docker-nodes'
  }
  stages {
    stage('Build') {
      agent {
        docker 'bernddoser/docker-devel-cpp:ubuntu-16.04-cmake-3.7.2-gcc-4.9-gtest-1.8.0'
        label 'docker-nodes'
      }
      steps {
        sh 'mkdir -p build'
        sh 'cd build && cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DGMX_BUILD_MDRUN_ONLY=OFF -DGMX_BUILD_FDA=ON -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_fda -DGMX_SIMD=NONE -DGMX_BUILD_UNITTESTS=ON -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF ..'
        sh 'cd build && make'
      }
    }
    stage('Test') {
      steps {
        sh 'cd build && make check'
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
      mail to:'bernd.doser@h-its.org', subject:"FAILURE: ${currentBuild.fullDisplayName}", body: "Boo, we failed."
    }
  }
}
