 #!/bin/bash
 
if [ -z "$1" ]
then
    echo staring build for branch master
    BRANCH=master
else
    echo staring build for branch $1
    BRANCH=$1
fi

## Currently does nothing, as the build script calls to -j$(nproc) 
if [ -z "$2" ]
then
    echo using 4 cpus
    CPUS=4
else
    echo using $2 cpus
    CPUS=$2
fi

cd /workspace/kratos
git clone --depth 1 --single-branch -b $BRANCH https://github.com/KratosMultiphysics/Kratos.git

cd /workspace/kratos/Kratos/scripts/wheels/linux/
chmod +x build_mpi.sh
./build_mpi.sh