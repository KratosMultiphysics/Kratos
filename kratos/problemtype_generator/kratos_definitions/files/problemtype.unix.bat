#!/bin/bash
# OutputFile: $2/$1.info
# ErrorFile: $2/$1.err

if [ -f ${HOME}/.bashrc ]; then
 source ${HOME}/.bashrc
fi

#delete previous data and settings 
rm -f "${2}/${1}.info"
rm -f "${2}/${1}.err"
rm -f "${2}/${1}.mdpa"
rm -f "${2}/${1}.kpt"
rm -f "${2}/problem_settings.py"

mv "${2}/${1}.dat" "${2}/${1}.mdpa"
mv "${2}/${1}-1.dat" "$2/${1}.kpt"
rm "${2}/${1}-2.dat"
mv "${2}/${1}-3.dat" "${2}/problem_settings.py"

# Read additional settings (kpt file)
while read name value; do
 if [ ${name} = "CUSTOMFILE" ]; then
   script_type=${value}
 elif [ ${name} = "FILEPATH" ]; then
   script_path=${value}
 fi;
done < "${2}/${1}.kpt"

if [ ${script_type} == "Use_Default" ]; then
 cp "${3}/script.py" "${2}/"
# cp "${3}/run_example_trilinos.py" "${2}/"
elif [ $script_type == "Copy_From" ]; then
 cp "$script_path" "${2}/script.py"
fi

# Launch the script
if [ -f "${2}/script.py" ]; then
 python "${2}/script.py" > "${2}/${1}.info" 2> "${2}/${1}.err"
fi


