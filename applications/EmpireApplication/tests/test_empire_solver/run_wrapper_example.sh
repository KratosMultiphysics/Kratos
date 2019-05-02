#!/bin/bash
# This script is used to execute the example

# ATTENTION make sure to set the EMPIRE-environment before running this script (e.g. startEMPIRE)
if [ -z ${EMPIRE_API_LIBSO_ON_MACHINE+x} ]; then
    echo "EMPIRE Environment in not set!";
    return 1
fi

# Options for how to run the example
# use_intel_mpi:
# true: If standard OpenMPI is used
# false: If Intel MPI or OpenMPI compiled with the "–disable-dlopen" option is used
use_intel_mpi=false

# run_in_xterm
# true: spawn three xterm terminals for each process (1 for EMPIRE and 2 for the Kratos Clients)
# false: writes output into files
run_in_xterm=true

# ompi-server path: this is the ABSOLUTE path to the ompi-server-file.
# Unused if Intel MPI or OpenMPI compiled with the "–disable-dlopen" option is used
# To start a server: "ompi-server -r ~/ompi_server.port"
ompi_server_file_path="/home/philippb/.ompi_server.port"

if $run_in_xterm; then
    if $use_intel_mpi; then
        echo "Using Intel MPI or OpenMPI compiled with the \"–disable-dlopen\" option and running with XTerm"
        xterm -hold -e mpiexec -np 1 Emperor emperorInput.xml & EPID=$!
        sleep 0.8s
        xterm -hold -e mpiexec -np 1 python3 fluid_dynamics_analysis_with_empire.py &
        sleep 0.8s
        xterm -hold -e mpiexec -np 1 python3 MainCoSim.py &
        wait $EPID
    else
        echo "Using OpenMPI and running with XTerm"
        xterm -hold -e mpiexec -np 1 --ompi-server file:$ompi_server_file_path Emperor emperorInput.xml & EPID=$!
        sleep 0.8s
        xterm -hold -e mpiexec -np 1 --ompi-server file:$ompi_server_file_path python3 fluid_dynamics_analysis_with_empire.py &
        sleep 0.8s
        xterm -hold -e mpiexec -np 1 --ompi-server file:$ompi_server_file_path python3 MainCoSim.py &
        wait $EPID
    fi
else
    if $use_intel_mpi; then
        echo "Using Intel MPI or OpenMPI compiled with the \"–disable-dlopen\" option and writting output to logfiles"
        mpiexec -np 1 Emperor emperorInput.xml > Emperor.log 2> Emperor.err & EPID=$!
        sleep 0.4s
        mpiexec -np 1 python3 fluid_dynamics_analysis_with_empire.py > fluid_dynamics_analysis_with_empire.log 2>fluid_dynamics_analysis_with_empire.err &
        sleep 0.4s
        mpiexec -np 1 python3 MainCoSim.py > MainCoSim.log 2>MainCoSim.err &
        wait $EPID
    else
        echo "Using OpenMPI and writting output to logfiles"
        mpiexec -np 1 --ompi-server file:$ompi_server_file_path Emperor emperorInput.xml > Emperor.log 2> Emperor.err & EPID=$!
        sleep 0.4s
        mpiexec -np 1 --ompi-server file:$ompi_server_file_path python3 fluid_dynamics_analysis_with_empire.py > fluid_dynamics_analysis_with_empire.log 2>fluid_dynamics_analysis_with_empire.err &
        sleep 0.4s
        mpiexec -np 1 --ompi-server file:$ompi_server_file_path python3 MainCoSim.py > MainCoSim.log 2>MainCoSim.err &
        wait $EPID
    fi
fi

