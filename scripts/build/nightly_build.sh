#!/bin/bash

## Step0: Define
cd /tmp
mkdir logs

LOG_DIR=/tmp/logs
MAIL_GCC=${HOME}/mail_gcc
MAIL_CLANG=${HOME}/mail_clang
MAIL_TO=${MAIL_TO_ADDRESS}

# Indicate that we want the stacktraces of crashes.
export LIBC_FATAL_STDERR_=1

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ubuntu/Kratos/bin/Linux/Custom/libs:/home/ubuntu/CompiledLibs/clang-3.8.0-16.04-prebuilt/lib
export PYTHONPATH=$PYTHONPATH:/home/ubuntu/Kratos/bin/Linux/Custom

## Step0: Prevent apt-get from triggering before the unasisted updates have finished
while fuser /var/lib/dpkg/lock >/dev/null 2>&1 ; do
  echo "Waiting for apt-get to finish"
done

# Download additional dependencies
sudo apt-get install -y unzip python3-h5py libhdf5-dev libio-socket-ssl-perl  libdigest-hmac-perl  libterm-readkey-perl libmime-lite-perl libfile-libmagic-perl libio-socket-inet6-perl python3-numpy python3-scipy

# We move to home directory
cd ${HOME}
# mmg library
wget https://drive.google.com/uc\?export\=download\&id\=1RM-SpxbynS6_zby8TjaJvIEC9NA-IxRZ -O mmg.zip
unzip mmg.zip
# Eigen library
wget https://bitbucket.org/eigen/eigen/get/dbed8786ceed.tar.gz -O eigen.tar.gz
tar xzf eigen.tar.gz
mv ${HOME}/eigen-eigen-dbed8786ceed ${HOME}/eigen
# ANurbs library
# a specific commit is specified, this has to be tested before updating
ANUBS_COMMIT_HASH=aa59b9ea2ff2c8e2fec321807e3eaf43c4394070
wget https://github.com/oberbichler/ANurbs/archive/${ANUBS_COMMIT_HASH}.tar.gz -O AnurbsLibrary.tar.gz
tar xzf AnurbsLibrary.tar.gz
mv ${HOME}/ANurbs-${ANUBS_COMMIT_HASH} ${HOME}/ANurbs

## Step1: Prepare
wget http://www.logix.cz/michal/devel/smtp-cli/smtp-cli
chmod 777 smtp-cli

wget https://github.com/KratosMultiphysics/Kratos/archive/master.tar.gz -O KratosMaster.tar.gz
tar xzf KratosMaster.tar.gz
mv ${HOME}/Kratos-master ${HOME}/Kratos
cd ${HOME}/Kratos

mkdir -p cmake_gcc
mkdir -p cmake_clang

## Copy the already prepared configure files
cp ${HOME}/Kratos/scripts/build/nightly/configure_gcc.sh ${HOME}/Kratos/configure_gcc.sh
cp ${HOME}/Kratos/scripts/build/nightly/configure_clang.sh ${HOME}/Kratos/configure_clang.sh

## Step2: Gcc

# Build
cd ${HOME}/Kratos
sh configure_gcc.sh > ${LOG_DIR}/configure_gcc.log 2>&1

# UnitTesting
cd ${HOME}/Kratos/bin/Release/KratosMultiphysics
python3 run_tests.py -l nightly -c python3 > ${LOG_DIR}/unittest_gcc.log 2>&1

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
cat ${LOG_DIR}/configure_gcc.log >> ${MAIL_GCC};
echo "\n" >> ${MAIL_GCC}
echo "UnitTests:    \n" >> ${MAIL_GCC}
echo "============= \n" >> ${MAIL_GCC}
cat ${LOG_DIR}/unittest_gcc.log >> ${MAIL_GCC};

cd ${HOME}
tar -zcvf /tmp/logs_gcc.tar.gz ${LOG_DIR}/*
./smtp-cli --host=${KRATOS_MAIL_SERVER} --enable-auth --user ${KRATOS_MAIL_USER} --password ${KRATOS_MAIL_PASSWD}  --to ${MAIL_TO} --body-plain ${MAIL_GCC} --attach /tmp/logs_gcc.tar.gz --subject "Kratos Nightly Report" --mail-from "kratosmultiphysics@gmail.com" --from "Kratos Nightly Report GCC"

## Give the email some time to be queued and delivered
sleep 300 # 5 minutes

## ------------------------------------------------------ ##

## Cleanup previous installation

rm ${LOG_DIR}/*
cd ${HOME}/Kratos
rm -rf libs

## ------------------------------------------------------ ##

## Step3: Clang

cd ${HOME}/Kratos
sh configure_clang.sh > ${LOG_DIR}/configure_clang.log 2>&1

# UnitTesting
cd ${HOME}/Kratos/bin/Release/KratosMultiphysics
python3 run_tests.py -l nightly -c python3 > ${LOG_DIR}/unittest_clang.log 2>&1

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
cat ${LOG_DIR}/configure_clang.log >> ${MAIL_CLANG};
echo "\n" >> ${MAIL_CLANG}
echo "UnitTests:    \n" >> ${MAIL_CLANG}
echo "============= \n" >> ${MAIL_CLANG}
cat ${LOG_DIR}/unittest_clang.log >> ${MAIL_CLANG};

cd ${HOME}
tar -zcvf /tmp/logs_clang.tar.gz ${LOG_DIR}/*
./smtp-cli --host=${KRATOS_MAIL_SERVER} --enable-auth --user ${KRATOS_MAIL_USER} --password ${KRATOS_MAIL_PASSWD}  --to ${MAIL_TO} --body-plain ${MAIL_CLANG} --attach /tmp/logs_clang.tar.gz --subject "Kratos Nightly Report" --mail-from "kratosmultiphysics@gmail.com" --from "Kratos Nightly Report CLANG"

## Step4: Finish

## Give the email some time to be queued and delivered
sleep 300 # 5 minutes

## Stop the machine
shutdown -h now
