#!/bin/sh

# USED BY GiD
# OutputFile: "$1.info"
# ErrorFile: "$1.err"

# Information
# basename = $1
# currentdirectory = $2
# problemtypedirectory = $3

# Setting paths. WARNING: one must properly set them before running this file

# Setting PATHs for Kratos
export PYTHONPATH="$HOME/Kratos:$PYTHONPATH"

## Linux
# Setting PATHs for Kratos libs
export LD_LIBRARY_PATH="$HOME/Kratos/libs:$HOME/boost_1_62_0/stage/lib:$LD_LIBRARY_PATH"
# Setting PATH for Python 3
export PATH="/usr/bin/python3:$PATH"

## OS X
#~ # Setting PATHs for Kratos
#~ export DYLD_LIBRARY_PATH="$HOME/Kratos/libs:$HOME/boost_1_62_0/stage/lib:$DYLD_LIBRARY_PATH"
#~ # Setting PATH for Python 3
#~ export PATH="/Library/Frameworks/Python.framework/Versions/3.5/bin:${PATH}"

# Execute the program
"python3" "MainKratos.py" > "$1.info" 2> "$1.err"
