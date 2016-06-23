#!/bin/bash

## Step0: Define
LOG_DIR=/tmp
MAIL_GCC=${HOME}/mail_gcc
MAIL_CLANG=${HOME}/mail_clang
MAIL_TO=${MAIL_TO_ADDRESS}

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ubuntu/Kratos/libs:/home/ubuntu/CompiledLibs/clang-3.8.0-16.04-prebuilt/lib
export PYTHONPATH=$PYTHONPATH:/home/ubuntu/Kratos

## Step1: Prepare
mkdir Kratos

cd ${HOME}/Kratos
svn co https://svn.cimne.upc.edu/p/kratos/kratos .

mkdir -p cmake_gcc
mkdir -p cmake_clang

## Copy the already prepared configure files
cp ${HOME}/Kratos/scripts/build/nightly/configure_gcc.sh ${HOME}/Kratos/cmake_gcc/configure.sh
cp ${HOME}/Kratos/scripts/build/nightly/configure_clang.sh ${HOME}/Kratos/cmake_clang/configure.sh

## Step2: Gcc

# Build
cd ${HOME}/Kratos/cmake_gcc
sh configure.sh > ${LOG_DIR}/configure_gcc.log 2> ${LOG_DIR}/configure_gcc.log
make install -j2 -k > ${LOG_DIR}/compile_gcc.log 2> ${LOG_DIR}/compile_gcc.log

# UnitTesting
cd ${HOME}/Kratos/kratos/python_scripts
python3 run_tests.py > ${LOG_DIR}/unittest_gcc.log 2> ${LOG_DIR}/unittest_gcc.log

# Benchmarking
cd ${HOME}/Kratos/benchmarking
python3 run_all_benchmarks.py > ${LOG_DIR}/benchmarking_gcc.log 2> ${LOG_DIR}/benchmarking_gcc.log

echo "This is the Kratos Nightly report for GCC" > ${MAIL_GCC}
echo "=========================================" >> ${MAIL_GCC}
echo "\n" >> ${MAIL_GCC}
echo -n "Failed Applications : "  >> ${MAIL_GCC};
cat ${LOG_DIR}/compile_gcc.log | grep "CMake.*all' failed" | awk '{print $5}' | awk -F'/' '{print $2}' | wc -l >> ${MAIL_GCC};
cat ${LOG_DIR}/compile_gcc.log | grep "CMake.*all' failed" | awk '{print $5}' | awk -F'/' '{print $2 " --> " $4}' >> ${MAIL_GCC};
echo -n "Compiling Errors    : "  >> ${MAIL_GCC};
cat ${LOG_DIR}/compile_gcc.log | grep " error:" | wc -l >> ${MAIL_GCC};
cat ${LOG_DIR}/compile_gcc.log | grep " error:" >> ${MAIL_GCC};
echo -n "Compiling Warnings  : "  >> ${MAIL_GCC};
cat ${LOG_DIR}/compile_gcc.log | grep " warning:" | wc -l >> ${MAIL_GCC};
cat ${LOG_DIR}/compile_gcc.log | grep " warning:" >> ${MAIL_GCC};
echo -n "Linking Error       : "  >> ${MAIL_GCC};
grep -c "cannot find" ${LOG_DIR}/compile_gcc.log >> ${MAIL_GCC};
grep "cannot find" ${LOG_DIR}/compile_gcc.log >> ${MAIL_GCC};
echo "\n" >> ${MAIL_GCC}
echo "UnitTests:    \n" >> ${MAIL_GCC}
echo "============= \n" >> ${MAIL_GCC}
cat ${LOG_DIR}/unittest_gcc.log >> ${MAIL_GCC};
echo "\n" >> ${MAIL_GCC}
echo "Benchmarking: \n" >> ${MAIL_GCC}
echo "============= \n" >> ${MAIL_GCC}
cat ${LOG_DIR}/benchmarking_gcc.log >> ${MAIL_GCC};

tar -zcvf ${LOG_DIR}/logs_gcc.tar.gz ${LOG_DIR}/configure_gcc.log ${LOG_DIR}/compile_gcc.log ${LOG_DIR}/benchmarking_gcc.log ${LOG_DIR}/unittest_gcc.log
sendEmail -s ${KRATOS_MAIL_SERVER} -xu ${KRATOS_MAIL_USER} -xp ${KRATOS_MAIL_PASSWD} -f "Kratos@no-reply.com" -t "${MAIL_TO}" -o message-file=${MAIL_GCC} -u "Kratos Nightly Report GCC" -a ${LOG_DIR}/logs_gcc.tar.gz

## ------------------------------------------------------ ##

## Cleanup previous installation

cd ${HOME}/Kratos
rm -rf libs

## ------------------------------------------------------ ##

## Step3: Clang

cd ${HOME}/Kratos/cmake_clang
sh configure.sh > ${LOG_DIR}/configure_clang.log 2> ${LOG_DIR}/configure_clang.log
make install -j2 -k > ${LOG_DIR}/compile_clang.log 2> ${LOG_DIR}/compile_clang.log

# UnitTesting
cd ${HOME}/Kratos/kratos/python_scripts
python3 run_tests.py > ${LOG_DIR}/unittest_clang.log 2> ${LOG_DIR}/unittest_clang.log

# Benchmarking
cd ${HOME}/Kratos/benchmarking
python3 run_all_benchmarks.py > ${LOG_DIR}/benchmarking_clang.log 2> ${LOG_DIR}/benchmarking_clang.log

echo "This is the Kratos Nightly report for CLANG" > ${MAIL_CLANG}
echo "===========================================" >> ${MAIL_CLANG}
echo "\n" >> ${MAIL_CLANG}
echo -n "Failed Applications : "  >> ${MAIL_CLANG};
cat ${LOG_DIR}/compile_clang.log | grep "CMake.*all' failed" | awk '{print $5}' | awk -F'/' '{print $2}' | wc -l >> ${MAIL_CLANG};
cat ${LOG_DIR}/compile_clang.log | grep "CMake.*all' failed" | awk '{print $5}' | awk -F'/' '{print $2 " --> " $4}' >> ${MAIL_CLANG};
echo -n "Compiling Errors    : "  >> ${MAIL_CLANG};
cat ${LOG_DIR}/compile_clang.log | grep " error:" | wc -l >> ${MAIL_CLANG};
cat ${LOG_DIR}/compile_clang.log | grep " error:" >> ${MAIL_CLANG};
echo -n "Compiling Warnings  : "  >> ${MAIL_CLANG};
cat ${LOG_DIR}/compile_clang.log | grep " warning:" | wc -l >> ${MAIL_CLANG};
cat ${LOG_DIR}/compile_clang.log | grep " warning:" >> ${MAIL_CLANG};
echo -n "Linking Error       : "  >> ${MAIL_CLANG};
grep -c "cannot find" ${LOG_DIR}/compile_clang.log >> ${MAIL_CLANG};
grep "cannot find" ${LOG_DIR}/compile_clang.log >> ${MAIL_CLANG};
echo "\n" >> ${MAIL_CLANG}
echo "UnitTests:    \n" >> ${MAIL_CLANG}
echo "============= \n" >> ${MAIL_CLANG}
cat ${LOG_DIR}/unittest_clang.log >> ${MAIL_CLANG};
echo "\n" >> ${MAIL_CLANG}
echo "Benchmarking: \n" >> ${MAIL_CLANG}
echo "============= \n" >> ${MAIL_CLANG}
cat ${LOG_DIR}/benchmarking_clang.log >> ${MAIL_CLANG};

tar -zcvf ${LOG_DIR}/logs_clang.tar.gz ${LOG_DIR}/configure_clang.log ${LOG_DIR}/compile_clang.log ${LOG_DIR}/benchmarking_clang.log ${LOG_DIR}/unittest_clang.log
sendEmail -s ${KRATOS_MAIL_SERVER} -xu ${KRATOS_MAIL_USER} -xp ${KRATOS_MAIL_PASSWD} -f "Kratos@no-reply.com" -t "${MAIL_TO}" -o message-file=${MAIL_CLANG} -u "Kratos Nightly Report Clang" -a ${LOG_DIR}/logs_clang.tar.gz

## Step4: Finish

## Give the email some time to be queued and delivered
sleep 150 # 2.5 minutes

## Stop the machine
shutdown -h now
