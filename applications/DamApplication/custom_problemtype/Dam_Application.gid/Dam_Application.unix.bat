#!/bin/bash
# OutputFile: $2/$1.info
# ErrorFile: $2/$1.err

# Identification for arguments
# basename               = ${1}
# Project directory      = ${2}
# Problem Type directory = ${3}


rm -f "${2}/${1}.info"
rm -f "${2}/${1}.err"

mv "${2}/${1}.dat" "${2}/${1}.mdpa"
mv "${2}/${1}-1.dat" "${2}/ProjectParameters.py"
cp "${3}/../../python_scripts/thermo_mechanic_script.py" "${2}/"

# WARNING: one must properly set the following paths before running this file

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


python3 "${2}/thermo_mechanic_script.py" > "${2}/${1}.info" 2> "${2}/${1}.err"
