#!/bin/bash

source $HOME/.bashrc

write_python_file="*GenData(Python_script_file)"
file_location="*GenData(Python_file)"
problemtype_name="*Tcl(GiD_Info Project ProblemType)"

# gid redefines LD_LIBRARY_PATH to its own libs directory
# and maintains OLD_LD_LIBRARY_PATH with previous settings
# therefore, we use the OLD_LD_LIBRARY_PATH and prepend the path to the kratos libs
#if [ "$OLD_LD_LIBRARY_PATH" != "" ]; then
#    export LD_LIBRARY_PATH="$3/kratos":"$3/kratos/libs":$OLD_LD_LIBRARY_PATH
#else
#    # do not add the ':'
#    export LD_LIBRARY_PATH="$3/kratos":"$3/kratos/libs"
#fi

export OMP_NUM_THREADS=*GenData(number_of_processors)

mv "$2/$1-3.dat" "$2/${problemtype_name}_var.py"
mv "$2/$1-4.dat" "$2/$1DEM_FEM_boundary.mdpa"

if [ $write_python_file = "Use_Default" ]
then
 cp "$3/script.py" "$2/"
 cp "$3/spheric_particle_script.py" "$2/"
 cp "$3/continuum_spheric_particle_script.py" "$2/"
 cp "$3/continuum_spheric_particle_script_mpi.py" "$2/"
elif [ $write_python_file = "Copy_From" ]
then
 cp "$file_location" "$2/script.py"
 cp "$file_location" "$2/spheric_particle_script.py"
 cp "$file_location" "$2/continuum_spheric_particle_script.py"
 cp "$file_location" "$2/continuum_spheric_particle_script_mpi.py"

fi

cp "$3/load_graf.dem" "$2/"

echo "Running on: " >& "$2/$1.info"
echo $OMP_NUM_THREADS >& "$2/$1.info"
echo "processors" >& "$2/$1.info"

if [ -f script.py ]
then
python "$2/script.py" >& "$2/$1.info"
fi
