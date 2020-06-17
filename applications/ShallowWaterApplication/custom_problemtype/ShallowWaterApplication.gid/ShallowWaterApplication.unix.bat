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
export LD_LIBRARY_PATH="$HOME/Kratos/bin/Release/libs:$LD_LIBRARY_PATH"

# Execute the program
"python3" "MainKratos.py" > "$1.info" 2> "$1.err"
