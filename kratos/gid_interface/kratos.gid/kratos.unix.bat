#!/bin/bash
# OutputFile: "$2/$1.info"
# ErrorFile: "$2/$1.err"
#delete previous result file 
rm -f "$2/$1*.post.bin" 
rm -f "$2/$1*.post.res" 
rm -f "$2/$1*.post.msh" 
rm -f "$2/$1.info"
rm -f "$2/$1.err"
rm -f "$2/$1.flavia.dat"

# include .bashrc if it exists
if [ -f "$HOME/.bashrc" ]; then
    . "$HOME/.bashrc"
fi

# gid redefines LD_LIBRARY_PATH to its own libs directory
# and maintains OLD_LD_LIBRARY_PATH with previous settings
# therefore, we use the OLD_LD_LIBRARY_PATH and prepend the path to the kratos libs
if [ "$OLD_LD_LIBRARY_PATH" != "" ]; then
    export LD_LIBRARY_PATH="$3/kratos":"$3/kratos/libs":$OLD_LD_LIBRARY_PATH
else
    # do not add the ':'
    export LD_LIBRARY_PATH="$3/kratos":"$3/kratos/libs"
fi

export PYTHONPATH="$3/kratos/python27.zip":"$3/kratos"

# Set the number of threads for OpenMP
export OMP_NUM_THREADS=$5

# Run Python using the script KratosOpenMP.py
"$3/kratos/runkratos" KratosOpenMP.py > "$2/$1.info" 2> "$2/$1.err"
