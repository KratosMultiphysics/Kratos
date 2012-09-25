#!/bin/bash -i
echo hola
echo par1 : $1
echo "1. param: $1" > /tmp/kk.txt
echo "2. param: $2" >> /tmp/kk.txt
echo "3. param: $3" >> /tmp/kk.txt
echo "4. param: $4" >> /tmp/kk.txt
echo "5. param: $5" >> /tmp/kk.txt
#    OutputFile: $2/$1.info
#    ErrorFile: $2/$1.err
#delete previous result file
rm -f "$2/$1.post.bin"
rm -f "$2/$1.info"
rm -f "$2/$1.err"
rm -f "$2/$1.flavia.dat"

mpirun -np $5 /usr/bin/python KratosMPI.py >"$2/$1.info" 2>"$2/$1.err"
