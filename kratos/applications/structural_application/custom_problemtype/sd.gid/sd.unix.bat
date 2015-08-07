#!/bin/sh -f
#    $1: name of the current project
#    $2: path of the current project
#    $3: path of the problem type
#    OutputFile: $2/$1.info
#    ErrorFile: $2/$1.err
#    delete previous result file 
rm -f $2/*.post.res 
rm -f $2/*.post.msh
rm -f $2/*.post.bin
# renaming Kratos input files
python $3/clean_mdpa.py $2/$1.dat $2/$1.bak
mv $2/$1.bak $2/$1.mdpa
rm $2/$1.dat
mv $2/$1-1.dat $2/${1}.py
mv $2/$1-2.dat $2/${1}_distributed_include.py
mv $2/$1-3.dat $2/${1}_layers.py
mv $2/$1-4.dat $2/${1}_shared_include.py
#check if ess file exist
if [ -f $2/$1.ess ] ; then pass ; else touch $2/$1.ess ; fi
#append ess to python script
cat $2/$1.ess >> $2/$1.py
#replace the string represent project name
sed s/rEpLaCeMeNtStRiNg/$1/g < $2/$1.py > $2/$1.py_changed
mv $2/$1.py_changed $2/$1.py
cd $2
cd ..
