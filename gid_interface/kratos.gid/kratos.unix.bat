#!/bin/bash -i
#    OutputFile: $2/$1.info
#    ErrorFile: $2/$1.err
#delete previous result file 
rm -f $2/$1.post.bin 
rm -f $2/$1.info
rm -f $2/$1.err
rm -f $2/$1.flavia.dat

python KratosOpenMP.py > $2/$1.info 2> $2/$1.err