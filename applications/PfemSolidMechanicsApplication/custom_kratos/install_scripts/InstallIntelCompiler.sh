#!/bin/bash

LIBC5_DEB=libstdc++5_3.3.6-20_amd64.deb
INTEL_COM=l_ccompxe_intel64_2011.10.319
INTEL_TAR=$INTEL_COM.tgz

echo -n "Install linux packages (y [yes], n [no] ) [default= y]: "
read -e LINUX_PACKAGES

if [ -z $LINUX_PACKAGES ]; then
    LINUX_PACKAGES="y"
fi

if [ $LINUX_PACKAGES == "y" ]; then

    sudo apt-get install -y gcc
    sudo apt-get install -y build-essential
    sudo apt-get install -y gcc-multilib
    sudo apt-get install -y g++
    sudo apt-get install -y rpm
    sudo apt-get install -y openjdk-6-jre-headless
    sudo apt-get install -y lib32stdc++6
    sudo apt-get install -y libc6-dev-i386
    sudo apt-get install -y g++-multilib

fi

PROGRAM_DIR=$PWD 

echo -n "DECOMPRESS libstdc++5 in a ./lib folder (y [yes], n [no] ) [default= y]: "
read -e DECOMPRESS_LIBC5

if [ -z $DECOMPRESS_LIBC5 ]; then
    DECOMPRESS_LIBC5="y"
fi

if [ $DECOMPRESS_LIBC5 == "y" ]; then

    dpkg --extract $LIBC5_DEB ./lib
    cd usr/lib
    sudo cp libstdc++.so.5.0.7 /usr/lib32
    cd /usr/lib32
    ln -s libstdc++.so.5.0.7 libstdc++.so.5

fi

echo "Decompression of $LIBC5_DEB Completed"

echo -n "SET libstdc++5 in a ./lib folder (y [yes], n [no] ) [default= y]: "
read -e SET_LIBC5

if [ -z $SET_LIBC5 ]; then
    SET_LIBC5="y"
fi

if [ $SET_LIBC5 == "y" ]; then

    cd lib/usr/lib
    sudo cp libstdc++.so.5.0.7 /usr/lib32
    cd /usr/lib32
    ln -s libstdc++.so.5.0.7 libstdc++.so.5

fi

echo "Set of $LIBC5_DEB Completed"

echo -n "DECOMPRESS intel compiler in a ./lib folder (y [yes], n [no] ) [default= y]: "
read -e DECOMPRESS_INTEL

if [ -z $DECOMPRESS_INTEL ]; then
    DECOMPRESS_INTEL="y"
fi

if [ $DECOMPRESS_INTEL == "y" ]; then

    tar -C lib/ -zxvf $INTEL_TAR >/dev/null 2>&1

fi

echo "Decompression of $INTEL_TAR Completed"

echo -n "INSTALL intel compiler (y [yes], n [no] ) [default= y]: "
read -e INTEL_INSTALL

if [ -z $INTEL_INSTALL ]; then
    INTEL_INSTALL="y"
fi

if [ $INTEL_INSTALL == "y" ]; then

    cd    lib/
    cd    $INTEL_COM
    sudo ./install.sh

fi

echo "Installation of $INTEL_COM Completed:"
echo "add this enviroment variables to .bashrc"
echo "source /opt/intel/bin/compilervars.sh intel64"
echo "source /opt/intel/mkl/bin/intel64/mklvars_intel64.sh intel64"
echo "source /opt/intel/mkl/bin/mklvars.sh intel64"

echo -n "DO IT AUTOMATICALLY (y [yes], n [no] ) [default= n]: "
read -e MODIFY_BASHRC

if [ -z $MODIFY_BASHRC ]; then
    MODIFY_BASHRC="n"
fi

if [ $MODIFY_BASHRC == "y" ]; then
    echo "source /opt/intel/bin/compilervars.sh intel64" >> $HOME/.bashrc
    echo "source /opt/intel/mkl/bin/intel64/mklvars_intel64.sh intel64" >> $HOME/.bashrc
    echo "source /opt/intel/mkl/bin/mklvars.sh intel64" >> $HOME/.bashrc
fi


echo
echo "********NOTE*********"
echo
echo "When compiling KRATOS: change mkl_solvers_application/CMakeLists.txt"
echo "Add: "
echo "set(USE_INTEL_GREATER_THAN_12 TRUE)"
echo 
echo "and cross your fingers when compiling!!"
echo
echo "*********************"
echo
APP_DIR=$HOME/kratos/applications/mkl_solvers_application/
emacs $APP_DIR/CMakeLists.txt &

