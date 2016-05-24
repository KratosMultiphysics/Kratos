#!/bin/sh

# USED BY GiD
# OutputFile: "$1.log"
# ErrorFile: "$1.err"

# Information
#basename = $1
#currentdirectory = $2
#problemtypedirectory = $3

# Remove old files
rm -f "$1.post.bin"
rm -f "$1.post.res"
rm -f "$1.post.msh"
rm -f "$1.log"
rm -f "$1.err"

# Update input files
mv "$1.dat" "$1.mdpa"
mv "$1-1.dat" "ProjectParameters.py"
cp "$3/../../python_scripts/main_script.py" "$2/"

# Setting paths. WARNING: one must properly set them before running this file

#Linux
#Setting PATHs for kratos
export PYTHONPATH="/path/to/kratos:$PYTHONPATH"
export LD_LIBRARY_PATH="/path/to/kratos/libs:/path/to/boost/stage/lib:$LD_LIBRARY_PATH"
#Setting PATH for Python 3
export PATH="/usr/bin/python3:$PATH"

#OS X
#Setting PATHs for kratos
#export PYTHONPATH="/path/to/kratos:$PYTHONPATH"
#export DYLD_LIBRARY_PATH="/path/to/kratos/libs:/path/to/boost/stage/lib:$DYLD_LIBRARY_PATH"
#Setting PATH for Python 3
#export PATH="/Library/Frameworks/Python.framework/Versions/3.5/bin:${PATH}"

# Execute the program
"python3" "main_script.py" > "$1.log" 2> "$1.err"
