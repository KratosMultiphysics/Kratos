#!/bin/bash -i

write_python_file="*GenData(Python_script_file)"
file_location="*GenData(Python_file)"
problemtype_name=*Tcl(GiD_Info Project ProblemType)

mv $2/$1-3.dat $2/${problemtype_name}_var.py

echo "problem_name=\"${1}\"" >> ${problemtype_name}_var.py
echo "problem_path=\"${2}\"" >> ${problemtype_name}_var.py
echo "kratos_path=\"${KRATOS_PATH}\"" >> ${problemtype_name}_var.py

if [ $write_python_file = "Use_Default" ]
then
 cp $3/script.py $2/
 cp $3/run_example_trilinos.py $2/
elif [ $write_python_file = "Copy_From" ]
then
 cp $file_location $2/script.py
fi

if [ $write_python_file = 1 ]
then
 cp $3/**.py $2/
fi

if [ -f script.py ]
then
 python $2/script.py >& $2/$1.info
fi
