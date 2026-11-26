NVARCH=`uname -s`_`uname -m`; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/23.5/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/23.5/compilers/bin:$PATH; export PATH
echo "NVARCH = $NVARCH"
echo "NVCOMPILERS = $NVCOMPILERS"
echo "PATH = $PATH"
