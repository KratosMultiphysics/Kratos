from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

kratos_path = '/home/philipp/software/KRATOS/kratos'
import sys
sys.path.append(kratos_path)

from KratosMultiphysics import *
from KratosMultiphysics.MappingApplication import *
import Mapper

from KratosMultiphysics.mpi import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.TrilinosApplication import *

#*******************************************************************************
##### START AUXILLARY FUNCTIONS #####
#*******************************************************************************

### Function to partition a ModelPart ###
def partition_and_read_model_part(model_part_name, model_part_input_file, size_domain, number_of_partitions = mpi.size):
    model_part = ModelPart(model_part_name)
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(VELOCITY)

    # number of partitions is by default equal to mpi.size
    if mpi.size > 1:
        if mpi.rank == 0:
            print("PARTITIONING", model_part_input_file)
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

### Function to set values on nodes ###
def set_values_on_nodes_id(model_part, variable, value, on_ghost_nodes, sub_model_part=None):
    if sub_model_part:
        for node in model_part.GetSubModelPart(sub_model_part).Nodes:
            value = node.Id
            if on_ghost_nodes:
                node.SetSolutionStepValue(variable, value)
            else:
                if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
                    node.SetSolutionStepValue(variable, value)
    else:
        for node in model_part.Nodes:
            value = node.Id
            if on_ghost_nodes:
                node.SetSolutionStepValue(variable, value)
            else:
                if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
                    node.SetSolutionStepValue(variable, value)
###
def set_values_on_nodes_pos(model_part, variable, value, on_ghost_nodes, sub_model_part=None):
    if sub_model_part:
        for node in model_part.GetSubModelPart(sub_model_part).Nodes:
            value = node.X**2 + node.Y**2
            if on_ghost_nodes:
                node.SetSolutionStepValue(variable, value)
            else:
                if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
                    node.SetSolutionStepValue(variable, value)
    else:
        for node in model_part.Nodes:
            value = node.X**2 + node.Y**2
            if on_ghost_nodes:
                node.SetSolutionStepValue(variable, value)
            else:
                if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
                    node.SetSolutionStepValue(variable, value)
###
def set_values_on_nodes(model_part, variable, value, on_ghost_nodes, sub_model_part=None):
    if sub_model_part:
        for node in model_part.GetSubModelPart(sub_model_part).Nodes:
            if on_ghost_nodes:
                node.SetSolutionStepValue(variable, value)
            else:
                if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
                    node.SetSolutionStepValue(variable, value)
    else:
        for node in model_part.Nodes:
            if on_ghost_nodes:
                node.SetSolutionStepValue(variable, value)
            else:
                if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
                    node.SetSolutionStepValue(variable, value)

### Function to set "INTERFACE"-Flag on nodes ###
def set_interface_flag(model_part, sub_model_part, on_ghost_nodes):
    for node in model_part.GetSubModelPart(sub_model_part).Nodes:
            if on_ghost_nodes:
                node.Set(INTERFACE)
            else:
                if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
                    node.Set(INTERFACE)

### Function to write meshes ###
def write_mesh(model_part, gid_io):
    gid_io.InitializeMesh(0)
    gid_io.WriteMesh(model_part.GetMesh())
    gid_io.FinalizeMesh()

def WriteNodalResultsCustom(gid_io, variable, model_part, output_time):
    if mpi.size > 1:
        gid_io.WriteNodalResults(variable, model_part.GetCommunicator().LocalMesh().Nodes, output_time, 0)
    else:
        gid_io.WriteNodalResults(variable, model_part.Nodes, output_time, 0)

#*******************************************************************************
##### END AUXILLARY FUNCTIONS #####
#*******************************************************************************

# Mdpa Input files
input_file_structure = "FSI_Example4Mapper_1_Structural"
input_file_fluid     = "FSI_Example4Mapper_1_Fluid"

parameter_file = open("MainKratos_Iterative_Mortar_Mapper_Parallel.json",'r')
ProjectParameters = Parameters( parameter_file.read())

# Create and Partition Model Parts
model_part_master  = partition_and_read_model_part("ModelPartNameMaster", input_file_fluid, 3)
model_part_slave = partition_and_read_model_part("ModelPartNameSlave", input_file_structure, 3)

# Initialize GidIO
output_file_master = "out_mapTest_master_r" + str(mpi.rank)
output_file_slave  = "out_mapTest_slave_r" + str(mpi.rank)

gid_mode = GiDPostMode.GiD_PostAscii
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions

gid_io_master = GidIO(output_file_master, gid_mode, multifile, deformed_mesh_flag, write_conditions)
gid_io_slave  = GidIO(output_file_slave,  gid_mode, multifile, deformed_mesh_flag, write_conditions)

# Initialize Results Output
gid_io_master.InitializeResults(0, model_part_master.GetMesh())
gid_io_slave.InitializeResults( 0, model_part_slave.GetMesh())

# Print original meshes
write_mesh(model_part_master, gid_io_master)
write_mesh(model_part_slave,  gid_io_slave)

# Write Initial values
write_time = 0.0
initial_nodal_values_master = 10.0
set_values_on_nodes(model_part_master, PRESSURE,   initial_nodal_values_master, False)
set_values_on_nodes(model_part_master, VELOCITY_X, initial_nodal_values_master, False)
set_values_on_nodes(model_part_master, VELOCITY_Y, initial_nodal_values_master, False)
set_values_on_nodes(model_part_master, VELOCITY_Z, initial_nodal_values_master, False)
WriteNodalResultsCustom(gid_io_master, PRESSURE, model_part_master, write_time)
WriteNodalResultsCustom(gid_io_master, VELOCITY, model_part_master, write_time)


initial_nodal_values_slave =  15.0
set_values_on_nodes(model_part_slave, PRESSURE,   initial_nodal_values_slave, False)
set_values_on_nodes(model_part_slave, VELOCITY_X, initial_nodal_values_slave, False)
set_values_on_nodes(model_part_slave, VELOCITY_Y, initial_nodal_values_slave, False)
set_values_on_nodes(model_part_slave, VELOCITY_Z, initial_nodal_values_slave, False)
WriteNodalResultsCustom(gid_io_slave, PRESSURE, model_part_slave, write_time)
WriteNodalResultsCustom(gid_io_slave, VELOCITY, model_part_slave, write_time)



##### Mapper stuff #####
nearestNeighborMapper = Mapper.Mapper(model_part_master, model_part_slave, ProjectParameters)

# Time step 1.0: Map; constant values
write_time = 1.0
initial_nodal_values_master = 11.0
set_values_on_nodes(model_part_master, PRESSURE,   initial_nodal_values_master, False)
set_values_on_nodes(model_part_master, VELOCITY_X, initial_nodal_values_master, False)
set_values_on_nodes(model_part_master, VELOCITY_Y, initial_nodal_values_master, False)
set_values_on_nodes(model_part_master, VELOCITY_Z, initial_nodal_values_master, False)
WriteNodalResultsCustom(gid_io_master, PRESSURE, model_part_master, write_time)
WriteNodalResultsCustom(gid_io_master, VELOCITY, model_part_master, write_time)

nearestNeighborMapper.Map(PRESSURE,PRESSURE)
nearestNeighborMapper.Map(VELOCITY,VELOCITY)

WriteNodalResultsCustom(gid_io_slave, PRESSURE, model_part_slave, write_time)
WriteNodalResultsCustom(gid_io_slave, VELOCITY, model_part_slave, write_time)

# Time step 2.0: InverseMap; constant values
write_time = 2.0
initial_nodal_values_slave =  100.0
set_values_on_nodes(model_part_slave, PRESSURE,   initial_nodal_values_slave, False)
set_values_on_nodes(model_part_slave, VELOCITY_X, initial_nodal_values_slave, False)
set_values_on_nodes(model_part_slave, VELOCITY_Y, initial_nodal_values_slave, False)
set_values_on_nodes(model_part_slave, VELOCITY_Z, initial_nodal_values_slave, False)
WriteNodalResultsCustom(gid_io_slave, PRESSURE, model_part_slave, write_time)
WriteNodalResultsCustom(gid_io_slave, VELOCITY, model_part_slave, write_time)

nearestNeighborMapper.InverseMap(PRESSURE,PRESSURE)
nearestNeighborMapper.InverseMap(VELOCITY,VELOCITY)

WriteNodalResultsCustom(gid_io_master, PRESSURE, model_part_master, write_time)
WriteNodalResultsCustom(gid_io_master, VELOCITY, model_part_master, write_time)

# Time step 3.0: Map; non-constant values
write_time = 3.0
set_values_on_nodes_id(model_part_master, PRESSURE,   initial_nodal_values_master, False)
set_values_on_nodes_pos(model_part_master, VELOCITY_X, initial_nodal_values_master, False)
set_values_on_nodes_pos(model_part_master, VELOCITY_Y, initial_nodal_values_master, False)
set_values_on_nodes_pos(model_part_master, VELOCITY_Z, initial_nodal_values_master, False)
WriteNodalResultsCustom(gid_io_master, PRESSURE, model_part_master, write_time)
WriteNodalResultsCustom(gid_io_master, VELOCITY, model_part_master, write_time)

nearestNeighborMapper.Map(PRESSURE,PRESSURE)
nearestNeighborMapper.Map(VELOCITY,VELOCITY)

WriteNodalResultsCustom(gid_io_slave, PRESSURE, model_part_slave, write_time)
WriteNodalResultsCustom(gid_io_slave, VELOCITY, model_part_slave, write_time)

# Time step 4.0: InverseMap; non-constant values
write_time = 4.0
set_values_on_nodes_pos(model_part_slave, PRESSURE,   initial_nodal_values_slave, False)
set_values_on_nodes_id(model_part_slave, VELOCITY_X, initial_nodal_values_slave, False)
set_values_on_nodes_id(model_part_slave, VELOCITY_Y, initial_nodal_values_slave, False)
set_values_on_nodes_id(model_part_slave, VELOCITY_Z, initial_nodal_values_slave, False)
WriteNodalResultsCustom(gid_io_slave, PRESSURE, model_part_slave, write_time)
WriteNodalResultsCustom(gid_io_slave, VELOCITY, model_part_slave, write_time)
#
nearestNeighborMapper.InverseMap(PRESSURE,PRESSURE)
nearestNeighborMapper.InverseMap(VELOCITY,VELOCITY)

WriteNodalResultsCustom(gid_io_master, PRESSURE, model_part_master, write_time)
WriteNodalResultsCustom(gid_io_master, VELOCITY, model_part_master, write_time)



# Finalize Result Output
gid_io_master.FinalizeResults()
gid_io_slave.FinalizeResults()
