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
rm -f "$2/$1.post.bin"
rm -f "$2/$1.info"
rm -f "$2/$1.err"
rm -f "$2/$1.flavia.dat"


# gid redefined LD_LIBRARY_PATH to its own libs directory
# and maintains OLD_LD_LIBRARY_PATH with previous settings
# therefore, we use the OLD_LD_LIBRARY_PATH and prepend the path to the kratos libs
export LD_LIBRARY_PATH="$3/kratos/libs":$OLD_LD_LIBRARY_PATH

# Set the number of threads for OpenMP
export OMP_NUM_THREADS=1

mpirun -np $5 "$3/kratos/runkratos" KratosMPI.py >"$2/$1.info" 2>"$2/$1.err"
