 #!/bin/bash
 
if [ -z "$1" ];
then
    echo staring build for branch master
    BRANCH=master
else
    echo staring build for branch $1
    BRANCH=$1
fi

cd /workspace/kratos
git clone --depth 1 --single-branch -b $BRANCH https://github.com/KratosMultiphysics/Kratos.git

cd /workspace/kratos/Kratos/scripts/wheels/linux/
chmod +x build.sh
./build.sh
 

