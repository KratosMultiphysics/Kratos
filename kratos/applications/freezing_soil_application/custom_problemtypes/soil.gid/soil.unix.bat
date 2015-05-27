#!/bin/csh -f
#    OutputFile: $2/$1.info
#    ErrorFile: $2/$1.err
#    delete previous result file 
rm -f  $2/*.post.res 
rm -f $2/*.post.msh
rm -f $2/*.post.bin
# renaming Kratos input files
mv $2/$1.dat $2/$1.mdpa
mv $2/$1-1.dat $2/${1}.py
mv $2/$1-2.dat $2/${1}_include.py
touch $2/$1.ess
#touch $2/set_material_data.py
cat $2/$1.ess >> $2/$1.py
sed s/rEpLaCeMeNtStRiNg/$1/g < $2/$1.py > $2/$1.py_changed
mv $2/$1.py_changed $2/$1.py
cd $2
setenv OMP_NUM_THREADS 2
python $1.py $1 > $1.info
cd ..
