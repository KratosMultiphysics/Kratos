#!/bin/bash -i

write_python_file="Use_Default"
file_location="script.py"
problemtype_name="DEM_explicit_solver"

export OMP_NUM_THREADS=1

mv "$2/$1-3.dat" "$2/${problemtype_name}_var.py"

if [ $write_python_file = "Use_Default" ]
then
 cp "$3/script.py" "$2/"
# cp "$3/run_example_trilinos.py" "$2/"
elif [ $write_python_file = "Copy_From" ]
then
 cp "$file_location" "$2/script.py"
fi

cp "$3/script_mpi.py" "$2/"
cp "$3/DEM_procedures.py" "$2/"
cp "$3/pressure_script.py" "$2/"

echo "Running on: " >& "$2/$1.info"
echo $OMP_NUM_THREADS >& "$2/$1.info"
echo "processors" >& "$2/$1.info"

if [ -f script.py ]
then
 python "$2/script.py" >& "$2/$1.info"
fi
