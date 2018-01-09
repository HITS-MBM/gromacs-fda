#!groovy

pipeline {
  agent {
    dockerfile {
      filename 'Dockerfile-gcc-4.9'
      label 'docker-nodes'
    }
  }
  stages {
    stage('Build') {
      steps {
        sh '''
          mkdir -p build
          cd build
          cmake -DCMAKE_BUILD_TYPE=release \
                -DCMAKE_C_COMPILER=/usr/bin/gcc \
                -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
                -DGMX_BUILD_MDRUN_ONLY=OFF \
                -DGMX_BUILD_FDA=ON \
                -DGMX_DEFAULT_SUFFIX=OFF \
                -DGMX_BINARY_SUFFIX=_fda \
                -DGMX_SIMD=NONE \
                -DGMX_BUILD_UNITTESTS=ON \
                -DGMX_BUILD_OWN_FFTW=ON \
                ..
          make
        '''
      }
    }
    stage('Test') {
      steps {
        script {
          try {
            sh '''
              cd build
              make check
            '''
          } catch (err) {
            echo "Failed: ${err}"
          } finally {
            step([
              $class: 'XUnitBuilder',
              thresholds: [[$class: 'FailedThreshold', unstableThreshold: '1']],
              tools: [[$class: 'GoogleTestType', pattern: 'build/Testing/Temporary/*.xml']]
            ])
          }
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
    stage('Deploy') {
      when {
        branch "*/tags/*"
      }
      steps {
        sh '''
         ## Extract version and project name from CMakeLists.txt
         #version="$(grep project CMakeLists.txt | cut -d " " -f3)"
         #project="$(grep project CMakeLists.txt | cut -d " " -f1 | cut -d "(" -f2)"
         #modulefile_path=/hits/fast/jenkins/modules
         #installation_path=/hits/fast/jenkins/software

         ## Write modulefile
         #mkdir -p &{modulefile_path}/${project}
         #cat > ${installation_path}/${project}/${version} << EOF
         ##%Module
         #module-whatis "Gromacs FDA ${version}";
         #prepend-path PATH ${installation_path}/${version}
         #EOF
        '''
      }
    }
  }
  }
  post {
    success {
      archiveArtifacts artifacts: 'build/bin/gmx_fda', fingerprint: true
      mail to: 'bernd.doser@h-its.org', subject: "SUCCESS: ${currentBuild.fullDisplayName}", body: "All fine."
    }
    failure {
      mail to: 'bernd.doser@h-its.org', subject: "FAILURE: ${currentBuild.fullDisplayName}", body: "Failed."
    }
  }
}
