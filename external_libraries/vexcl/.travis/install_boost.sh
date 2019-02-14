#!/bin/bash

export BOOST_ROOT=${HOME}/${BOOST_BASENAME}

if [ ! -e ${BOOST_ROOT}/boost/config.hpp ]; then
    pushd ${HOME}
    BOOST_VERSION=$(echo $BOOST_BASENAME | awk -F '_' '{print $2 "." $3 "." $4 }')
    wget https://dl.bintray.com/boostorg/release/$BOOST_VERSION/source/$BOOST_BASENAME.tar.bz2
    rm -rf $BOOST_BASENAME
    tar xf ${BOOST_BASENAME}.tar.bz2
    cd ${BOOST_BASENAME}
    ./bootstrap.sh --with-libraries=atomic,chrono,date_time,filesystem,program_options,system,test,thread
    ./b2
    popd
fi
