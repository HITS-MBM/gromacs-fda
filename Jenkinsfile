#!groovy

pipeline {

  agent {
    dockerfile {
      label 'docker-gpu-host'
      filename 'devel/Dockerfile'
    }
  }

  options {
    timeout(time: 2, unit: 'HOURS')
    disableConcurrentBuilds()
  }

  stages {
    stage('Build') {
      steps {
        sh '''
          mkdir -p build
          cd build
          cmake -DGMX_BUILD_MDRUN_ONLY=OFF \
                -DGMX_BUILD_FDA=ON \
                -DGMX_DEFAULT_SUFFIX=OFF \
                -DGMX_BINARY_SUFFIX=_fda \
                -DGMX_SIMD=NONE \
                -DGMX_BUILD_UNITTESTS=ON \
                -DGMX_BUILD_OWN_FFTW=ON \
                ..
          make 2>&1 |tee make.out
        '''
      }
      post {
        always {
          recordIssues enabledForFailure: true, aggregatingResults: false,
            tool: gcc(id: 'gcc', pattern: 'build/make.out')
        }
      }
    }
    stage('Test') {
      steps {
        sh '''
          cd build
          make check
        '''
      }
      post {
        always {
          step([
            $class: 'XUnitPublisher',
            thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
            tools: [[$class: 'GoogleTestType', pattern: 'build/Testing/Temporary/*.xml']]
          ])
        }
      }
    }
    stage('Doxygen') {
      steps {
        sh '''
          cd build
          make doxygen-all
        '''
        publishHTML( target: [
          allowMissing: false,
          alwaysLinkToLastBuild: false,
          keepAll: true,
          reportName: 'Doxygen',
          reportDir: 'build/docs/html/doxygen/html-full',
          reportFiles: 'index.xhtml'
        ])
      }
    }
  }
  post {
    success {
      archiveArtifacts artifacts: 'build/bin/gmx_fda', fingerprint: true
    }
    failure {
      mail to: 'bernd.doser@h-its.org', subject: "FAILURE: ${currentBuild.fullDisplayName}", body: "Failed."
    }
  }
}
