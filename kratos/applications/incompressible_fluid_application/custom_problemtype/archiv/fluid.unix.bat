#!/bin/bash -f
#    OutputFile: $2/$1.info
#    ErrorFile: $2/$1.err
#delete previous result file 
rm -f  $2/$1.flavia.res 
rm $2/$1.info
rm $2/$1.flavia.dat
rm $2/$1.err
rm $2/$1.node
rm $2/$1.prop
rm $2/$1.elem
rm $2/$1.cond
rm $2/$1.init
mv $2/$1-1.dat $2/$1.prop
mv $2/$1-2.dat $2/$1.elem
mv $2/$1-3.dat $2/$1.cond
mv $2/$1-4.dat $2/$1.init
mv $2/$1-5.dat $2/$1.node
