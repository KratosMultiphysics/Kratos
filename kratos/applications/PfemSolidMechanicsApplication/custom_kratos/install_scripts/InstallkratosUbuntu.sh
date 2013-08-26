#!/bin/bash

echo -n "Install linux packages (y [yes], n [no] ) [default= y]: "
read -e LINUX_PACKAGES

if [ -z $LINUX_PACKAGES ]; then
    LINUX_PACKAGES="y"
fi

if [ $LINUX_PACKAGES == "y" ]; then

    sudo apt-get install -y build-essential
    sudo apt-get install -y subversion
    sudo apt-get install -y rapidsvn
    sudo apt-get install -y libblas3gf
    sudo apt-get install -y liblapack3gf
    sudo apt-get install -y libblas-dev
    sudo apt-get install -y liblapack-dev
    sudo apt-get install -y python-dev
    sudo apt-get install -y gfortran
    sudo apt-get install -y csh
    sudo apt-get install -y swig
    sudo apt-get install -y python-matplotlib
    sudo apt-get install -y python-scipy
    sudo apt-get install -y python-numpy
    sudo apt-get install -y python-scientific
    sudo apt-get install -y g++
    sudo apt-get install -y python
    sudo apt-get install libboost-all-dev
    sudo apt-get install libboost-dbg
    sudo apt-get install cmake

fi

PROGRAM_DIR=$HOME
PROGRAM_FOLDER=kratos

echo -n "INSTALL Kratos (y [yes], n [no] ) [default= y]: "
read -e KRATOS_INSTALL

if [ -z $KRATOS_INSTALL ]; then
    KRATOS_INSTALL="y"
fi

if [ $KRATOS_INSTALL == "y" ]; then

    cd cmake_build

    echo -n "INITIAL INSTALL (y [yes], n [no] ) [default= n]: "
    read -e INITIAL_INSTALL

    if [ -z $INITIAL_INSTALL ]; then
	INITIAL_INSTALL="n"
    fi
    
    if [ $INITIAL_INSTALL == "y" ]; then
	cp example_configure.sh.do_not_touch configure.sh
    fi

    sh configure.sh

fi

echo
echo "Installation of KRATOS Completed:"
echo "add this enviroment variables to .bashrc"
echo "export PYTHONPATH=$HOME/kratos/libs:$PYTHONPATH"
echo "export PYTHONPATH=$HOME/kratos:$PYTHONPATH"
echo

echo -n "DO IT AUTOMATICALLY (y [yes], n [no] ) [default= n]: "
read -e MODIFY_BASHRC

if [ -z $MODIFY_BASHRC ]; then
    MODIFY_BASHRC="n"
fi

if [ $MODIFY_BASHRC == "y" ]; then
    echo "export PYTHONPATH=\$HOME/kratos/libs:\$PYTHONPATH" >> $HOME/.bashrc
    echo "export PYTHONPATH=\$HOME/kratos:\$PYTHONPATH" >> $HOME/.bashrc
fi

