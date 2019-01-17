#!/bin/bash

# Performs an analysis of SpM source code:
# - we consider to be in SpM's source code root
# - we consider having the coverage file spm.lcov in the root directory
# - we consider having cppcheck, rats, sonar-scanner programs available in the environment

# filter sources:
# - consider generated files in ${BUILDDIR}
# - exclude base *z* files to avoid duplication
# - exclude cblas.h and lapacke-.h because not really part of spm and make cppcheck analysis too long

if [ $# -gt 0 ]
then
    BUILDDIR=$1
fi
BUILDDIR=${BUILDDIR:-build}

./.gitlab-ci-filelist.sh $BUILDDIR

# Generate coverage xml output
lcov_cobertura.py spm.lcov --output spm-coverage.xml

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"

# to get it displayed and captured by gitlab to expose the badge on the main page
cat ./spm-gcov.log

# run cppcheck analysis
cppcheck -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingInclude ${UNDEFINITIONS} --file-list=./filelist.txt 2> spm-cppcheck.xml

# run rats analysis
rats -w 3 --xml  `cat filelist.txt` > spm-rats.xml

# Set the default for the project key
SONARQUBE_PROJECTKEY=${SONARQUBE_PROJECTKEY:-hiepacs:spm:gitlab:dev}

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.bordeaux.inria.fr/sonarqube
sonar.login=$SONARQUBE_LOGIN

sonar.links.homepage=https://gitlab.inria.fr/solverstack/spm
sonar.links.scm=https://gitlab.inria.fr/solverstack/spm.git
sonar.links.ci=https://gitlab.inria.fr/solverstack/spm/pipelines
sonar.links.issue=https://gitlab.inria.fr/solverstack/spm/issues

sonar.projectKey=$SONARQUBE_PROJECTKEY
sonar.projectDescription=Parallel Sparse direct Solver
sonar.projectVersion=master

sonar.language=c
sonar.sources=$BUILDDIR/src, $BUILDDIR/tests, include, src, tests
sonar.inclusions=`cat filelist.txt | xargs echo | sed 's/ /, /g'`
sonar.sourceEncoding=UTF-8
sonar.c.errorRecoveryEnabled=true
sonar.c.compiler.charset=UTF-8
sonar.c.compiler.parser=GCC
sonar.c.compiler.regex=^(.*):(\\d+):\\d+: warning: (.*)\\[(.*)\\]$
sonar.c.compiler.reportPath=spm-build.log
sonar.c.coverage.reportPath=spm-coverage.xml
sonar.c.cppcheck.reportPath=spm-cppcheck.xml
sonar.c.rats.reportPath=spm-rats.xml
sonar.c.jsonCompilationDatabase=build/compile_commands.json
EOF

# run sonar analysis + publish on sonarqube-dev
sonar-scanner -X > sonar.log
