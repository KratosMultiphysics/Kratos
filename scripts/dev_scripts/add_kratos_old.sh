# add kratos to paths

KRATOS_ROOT="/home/fmeister/GitkrakenRepos/Kratos_old/bin/Release"

NUM_OF_OMP_THREADS=4

if [[ ":$PYTHONPATH:" != *"Kratos"* ]]; then
    export PYTHONPATH="$PYTHONPATH:$KRATOS_ROOT"
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$KRATOS_ROOT/libs"
    echo "Kratos added to paths"
else
    echo "Kratos Path already exists in Path."
    echo $PYTHONPATH
    echo $LD_LIBRARY_PATH
    echo "Should PYTHONPATH and LD_LIBRARY_PATH be cleared?"
    read -p "y/n: " YES_NO
    if [ "$YES_NO" = "y" ]; then
        export PYTHONPATH="$KRATOS_ROOT"
        export LD_LIBRARY_PATH="$KRATOS_ROOT/libs"
    else
        echo "Paths not overriden"
    fi
fi

echo $PYTHONPATH
echo $LD_LIBRARY_PATH

# set number of threads
export OMP_NUM_THREADS=$NUM_OF_OMP_THREADS
echo "Set number of threads to Variable: NUM_OF_OMP_THREADS = ${NUM_OF_OMP_THREADS}"

