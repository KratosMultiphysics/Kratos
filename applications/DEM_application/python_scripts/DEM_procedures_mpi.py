from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.mpi import *

def KRATOS_MPI_WORLD_BARRIER():
    mpi.world.barrier()

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
    model_part_io_solid = ModelPartIO(my_input_filename, True)

    return [model_part_io_solid, model_part]

def CreateDirectories(main_path,problem_name):

    KRATOS_MPI_WORLD_BARRIER()

    root             = main_path + '/' + problem_name

    post_path        = root + '_Post_Files'
    list_path        = root + '_Post_Lists'
    data_and_results = root + '_Results_and_Data'
    graphs_path      = root + '_Graphs'
    MPI_results      = root + '_MPI_results'

    if mpi.rank == 0:
        for directory in [post_path, list_path, data_and_results, graphs_path, MPI_results]:
            if not os.path.isdir(directory):
                os.makedirs(str(directory))

    KRATOS_MPI_WORLD_BARRIER()

    return [post_path,list_path,data_and_results,graphs_path,MPI_results]

def KRATOSprint(message):
    if (mpi.rank == 0):
        print(message)    
