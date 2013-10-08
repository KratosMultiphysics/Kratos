#!/bin/bash -i
#  echo hola
#  echo par1 : $1
#  echo "1. param: $1" > /tmp/kk.txt
#  echo "2. param: $2" >> /tmp/kk.txt
#  echo "3. param: $3" >> /tmp/kk.txt
#  echo "4. param: $4" >> /tmp/kk.txt
#  echo "5. param: $5" >> /tmp/kk.txt
#    OutputFile: "$2/$1.info"
#    ErrorFile: "$2/$1.err"
#delete previous result file
rm -f "$2/$1"*.post.bin
rm -f "$2/$1"*.post.res 
rm -f "$2/$1"*.post.msh
rm -f "$2/$1.info"
rm -f "$2/$1.err"


# Set the number of threads for OpenMP
export OMP_NUM_THREADS=1

export PYTHONPATH="$3/kratos/python27.zip":"$3/kratos":$PYTHONPATH

# Set mpi paths
export OPAL_PREFIX="$3/openmpi"
export PATH="$3/openmpi/bin":$PATH
export LD_LIBRARY_PATH="$3/openmpi/lib":$LD_LIBRARY_PATH

"$3/openmpi/bin/mpirun" -np $5 "$3/kratos/runkratos" KratosStructuralOpenMP.py >"$2/$1.info" 2>"$2/$1.err"
