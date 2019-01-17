#!/bin/bash

# Performs an analysis of PaStiX source code:
# - we consider to be in PaStiX's source code root
# - we consider having the coverage file pastix.lcov in the root directory
# - we consider having cppcheck, rats, sonar-scanner programs available in the environment

# filter sources:
# - consider generated files in ${BUILDDIR}
# - exclude base *z* files to avoid duplication
# - exclude cblas.h and lapacke-.h because not really part of pastix and make cppcheck analysis too long

if [ $# -gt 0 ]
then
    BUILDDIR=$1
fi
BUILDDIR=${BUILDDIR:-build}

./.gitlab-ci-filelist.sh $BUILDDIR

# Generate coverage xml output
lcov_cobertura.py pastix.lcov --output pastix-coverage.xml

# Undefine this because not relevant in our configuration
export UNDEFINITIONS="-UWIN32 -UWIN64 -U_MSC_EXTENSIONS -U_MSC_VER -U__SUNPRO_C -U__SUNPRO_CC -U__sun -Usun -U__cplusplus"
export UNDEFINITIONS="$UNDEFINITIONS -UPARSEC_PROF_DRY_BODY -UPARSEC_PROF_DRY_RUN -UPARSEC_PROF_TRACE -UPARSEC_PROF_GRAPHER -UPARSEC_SIM -DPINS_ENABLE -UBUILD_PARSEC"
export UNDEFINITIONS="$UNDEFINITIONS -UPARSEC_DEBUG_NOISIER -UPARSEC_DEBUG_PARANOID -UPARSEC_DEBUG_HISTORY -UPARSEC_C_PARSEC_HAVE_VISIBILITY"
export UNDEFINITIONS="$UNDEFINITIONS -UBUILDING_STRAPU"
export UNDEFINITIONS="$UNDEFINITIONS -UNAPA_SOPALIN -UPASTIX_WITH_STARPU_DIST"

# to get it displayed and captured by gitlab to expose the badge on the main page
cat ./pastix-gcov.log

# run cppcheck analysis
cppcheck -v -f --language=c --platform=unix64 --enable=all --xml --xml-version=2 --suppress=missingInclude ${UNDEFINITIONS} --file-list=./filelist-c.txt 2> pastix-cppcheck.xml

# run rats analysis
rats -w 3 --xml  `cat filelist.txt` > pastix-rats.xml

# Set the default for the project key
SONARQUBE_PROJECTKEY=${SONARQUBE_PROJECTKEY:-hiepacs:pastix:gitlab:dev}

# create the sonarqube config file
cat > sonar-project.properties << EOF
sonar.host.url=https://sonarqube.bordeaux.inria.fr/sonarqube
sonar.login=$SONARQUBE_LOGIN

sonar.links.homepage=https://gitlab.inria.fr/solverstack/pastix
sonar.links.scm=https://gitlab.inria.fr/solverstack/pastix.git
sonar.links.ci=https://gitlab.inria.fr/solverstack/pastix/pipelines
sonar.links.issue=https://gitlab.inria.fr/solverstack/pastix/issues

sonar.projectKey=$SONARQUBE_PROJECTKEY
sonar.projectDescription=Parallel Sparse direct Solver
sonar.projectVersion=master

sonar.language=c
sonar.sources=$BUILDDIR, bcsc, blend, common, example, graph, include, kernels, order, refinement, sopalin, symbol, test
sonar.inclusions=`cat filelist.txt | grep -v spm | xargs echo | sed 's/ /, /g'`
sonar.sourceEncoding=UTF-8
sonar.c.errorRecoveryEnabled=true
sonar.c.compiler.charset=UTF-8
sonar.c.compiler.parser=GCC
sonar.c.compiler.regex=^(.*):(\\d+):\\d+: warning: (.*)\\[(.*)\\]$
sonar.c.compiler.reportPath=pastix-build.log
sonar.c.coverage.reportPath=pastix-coverage.xml
sonar.c.cppcheck.reportPath=pastix-cppcheck.xml
sonar.c.rats.reportPath=pastix-rats.xml
sonar.c.jsonCompilationDatabase=build/compile_commands.json
EOF

# run sonar analysis + publish on sonarqube-dev
sonar-scanner -X > sonar.log
