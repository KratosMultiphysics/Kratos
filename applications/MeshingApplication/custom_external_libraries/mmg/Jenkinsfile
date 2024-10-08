#!groovy

// Jenkinsfile for compiling, testing, and packaging MMG


pipeline {
    agent {
        dockerfile true
    }
    stages {
        stage('Checkout') {
            steps {
                checkout scm
                sh 'mkdir -p build local'
            }
        }
        stage('Compile') {
            steps {
                sh '''
                	cd build &&
                	cmake -D CMAKE_BUILD_TYPE=Debug -D BUILD_TESTING=ON .. &&
                	make
                '''
            }
        }
        stage('Test') {
            steps {
                sh '''
                	cd build &&
                	ctest
                '''
            }
        }
        stage('Package') {
            steps {
                sh '''
                	cd build &&
                	cmake -D CMAKE_BUILD_TYPE=Release -D BUILD_TESTING=OFF -D CMAKE_INSTALL_PREFIX:PATH=$PWD/../local ..
                	make
                	make install
                '''
								archiveArtifacts artifacts: 'local/**', fingerprint: true, onlyIfSuccessful: true
            }
        }
    }
}
