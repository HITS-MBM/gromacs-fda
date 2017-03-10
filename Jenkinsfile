#!groovy

pipeline {
  agent {
    docker {
      image 'bernddoser/docker-devel-cpp:ubuntu-16.04-gcc-4.9-tools-1'
      label 'docker-nodes'
    }
  }
  stages {
    stage('Build') {
      steps {
        sh 'mkdir -p build'
        dir('build') {
          sh 'cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DGMX_BUILD_MDRUN_ONLY=OFF -DGMX_BUILD_FDA=ON -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_fda -DGMX_SIMD=NONE -DGMX_BUILD_UNITTESTS=ON -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=OFF ..'
          sh 'make'
        }
      }
    }
    stage('Test') {
      steps {
        dir('build') {
          sh 'make check'
        }
      }
    }
    stage('Doxygen') {
      steps {
        dir('build') {
          sh 'make doxygen-all'
        }
      }
    }
  }
  post {
    always {
      step([$class: 'XUnitBuilder',
        thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
        tools: [[$class: 'GoogleTestType', pattern: 'build/Testing/Temporary/*.xml']]])
        
      publishHTML target: [$class: 'HtmlPublisherTarget',
        reportName: 'Doxygen', reportDir: 'docs/html/doxygen/html-full', reportFiles: 'index.html']
   
      deleteDir()
    }
    success {
      archiveArtifacts artifacts: 'build/bin/gmx_fda', fingerprint: true
    }
//    failure {
//      mail to:'bernd.doser@h-its.org', subject:"FAILURE: ${currentBuild.fullDisplayName}", body: "Boo, we failed."
//    }
  }
}
