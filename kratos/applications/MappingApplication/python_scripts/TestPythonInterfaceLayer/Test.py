# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

from KratosMultiphysics import *
from KratosMultiphysics.MappingApplication import *

import Mapper
# import KratosMultiphysics.MappingApplication as hgfdsa
from KratosMultiphysics.mpi import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.TrilinosApplication import *

def partition_and_read_model_part(model_part_name, model_part_input_file, size_domain, variable_list, number_of_partitions = mpi.size):
    model_part = ModelPart(model_part_name)
    for variable in variable_list:
        model_part.AddNodalSolutionStepVariable(variable)

    # number of partitions is by default equal to mpi.size
    if mpi.size > 1:
        if mpi.rank == 0:
            model_part_io = ReorderConsecutiveModelPartIO(model_part_input_file)

            partitioner = MetisDivideHeterogeneousInputProcess(
                model_part_io,
                number_of_partitions,
                size_domain,
                0, # verbosity, set to 1 for more detailed output
                True)

            partitioner.Execute()

        mpi.world.barrier()
        model_part_input_file = model_part_input_file + "_" + str(mpi.rank)

    model_part_io = ModelPartIO(model_part_input_file)
    model_part_io.ReadModelPart(model_part)

    MPICommSetup = SetMPICommunicatorProcess(model_part)
    MPICommSetup.Execute()

    ParallelFillComm = ParallelFillCommunicator(model_part.GetRootModelPart())
    ParallelFillComm.Execute()

    model_part.ProcessInfo.SetValue(DOMAIN_SIZE, size_domain)
    model_part.SetBufferSize(1)

    return model_part


input_file_structure = "mdpa_files/FSI_Example4Mapper_1_Structural"
input_file_fluid     = "mdpa_files/FSI_Example4Mapper_1_Fluid"

parameter_file = open("TestPythonInterfaceLayer.json",'r')
ProjectParameters = Parameters( parameter_file.read())

# Create and Partition Model Parts
variable_list = [PRESSURE, VELOCITY, PARTITION_INDEX]
model_part_origin  = partition_and_read_model_part("ModelPartNameOrigin", input_file_fluid, 3, variable_list)
model_part_destination = partition_and_read_model_part("ModelPartNameDestination", input_file_structure, 3, variable_list)

mapper = Mapper.Mapper(model_part_origin, model_part_destination, ProjectParameters)

mapper.Map(PRESSURE, PRESSURE)
mapper.InverseMap(PRESSURE, PRESSURE)
