#!/bin/bash
#    OutputFile: $2/$1.info
#    ErrorFile: $2/$1.err

#delete previous data and settings 
rm -f "${2}/${1}.info"
rm -f "${2}/${1}.err"
rm -f "${2}/${1}.mdpa"
rm -f "${2}/ProjectParameters.py"

mv "${2}/${1}.dat" "${2}/${1}.mdpa"
mv "${2}/${1}-1.dat" "${2}/ProjectParameters.py"

cp "${3}/../../python_scripts/main_script.py" "${2}/"


# Warning: one must properly set the following paths before running this file

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
#Setting PATH for Python 3.5
#export PATH="/Library/Frameworks/Python.framework/Versions/3.5/bin:${PATH}"


# Launch the script
if [ -f "${2}/main_script.py" ]; then
 python3 "${2}/main_script.py" > "${2}/${1}.info" 2> "${2}/${1}.err"
fi
