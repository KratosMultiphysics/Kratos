#!/bin/bash -i

write_python_file=*GenData(Write_a_new_Python_script_file_for_Kratos,Int)
problemtype_name=*Tcl(GiD_Info Project ProblemType)

mv $2/$1-3.dat $2/${problemtype_name}_var.py

echo "problem_name=\"${1}\"" >> ${problemtype_name}_var.py
echo "problem_path=\"${2}\"" >> ${problemtype_name}_var.py
echo "kratos_path=\"${KRATOS_PATH}\"" >> ${problemtype_name}_var.py

if [ $write_python_file = 1 ]
then
 cp $3/**.py $2/
fi
