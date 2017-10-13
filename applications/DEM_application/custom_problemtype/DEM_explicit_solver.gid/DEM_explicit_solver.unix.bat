#!/bin/bash
#    OutputFile: $2/$1.info
#    ErrorFile: $2/$1.err
#delete previous result file

#rm -f "$2/${1}_*.post.bin"


rm -f "$2/$1_Post_Files/${1}_*.post.bin"
rm -f "$2/$1_Neigh_Lists/neigh_list_*"
rm -f "$2/$1.info"
rm -f "$2/$1.mdpa"

mv "$2/$1.dat" "$2/$1.mdpa"
mv "$2/$1-1.dat" "$2/${1}_aux.unix.bat"
rm -f "$2/$1-2.dat"



chmod 700 "$2/${1}_aux.unix.bat"

# The following line is required when python and kratos
# are linked to different versions of libstdc++.so
export LD_PRELOAD=libstdc++.so.6

./"${1}_aux.unix.bat" $1 $2 $3

