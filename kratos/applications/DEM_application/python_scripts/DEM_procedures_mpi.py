from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.MetisApplication import *

from DEM_explicit_solver_var import *
from pressure_script import *

from numpy import *

from KratosMultiphysics.mpi import *


def AddMpiVariables(model_part):

    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
    model_part.AddNodalSolutionStepVariable(PARTITION_MASK)


def PerformInitialPartition(model_part, model_part_io_solid, input_file_name):

    domain_size = 3

    print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "before performing the division")
    number_of_partitions = mpi.size  # we set it equal to the number of processors
    if mpi.rank == 0:
        print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "start partition process")
        partitioner = MetisDivideNodalInputToPartitionsProcess(model_part_io_solid, number_of_partitions, domain_size);
        partitioner.Execute()

    print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "division performed")
    mpi.world.barrier()

    MPICommSetup = SetMPICommunicatorProcess(model_part)
    MPICommSetup.Execute()

    print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "Comunicator Set")

    print("(" + str(mpi.rank) + "," + str(mpi.size) + ")" + "Reading: "+input_file_name+"_"+str(mpi.rank))

    my_input_filename = input_file_name + "_" + str(mpi.rank)
    model_part_io_solid = ModelPartIO(my_input_filename)

    return [model_part_io_solid, model_part]
