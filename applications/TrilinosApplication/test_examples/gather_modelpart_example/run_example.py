from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MetisApplication import *
import sys, math

if mpi.size != 2:
    print("Terminating due to mpi.size != 2")
    sys.exit()

model = Model()
model_part = model.CreateModelPart("LocalPart",1)
input_filename = "gather_model_part"
number_of_partitions = 2
domain_size = 2

# PRESSURE variable is used to demonstrate
model_part.AddNodalSolutionStepVariable(PRESSURE)
model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
if mpi.rank == 0:
    model_part_io = ModelPartIO(input_filename)
    partitioner = MetisDivideHeterogeneousInputProcess(
        model_part_io,
        number_of_partitions,
        domain_size,
        1,
        True)
    partitioner.Execute()
mpi.world.barrier()

MPICommSetup = SetMPICommunicatorProcess(model_part)
MPICommSetup.Execute()
my_input_filename = input_filename + "_" + str(mpi.rank)
model_part_io = ModelPartIO(my_input_filename)
model_part_io.ReadModelPart(model_part)
Comm = CreateCommunicator()
model_part.SetBufferSize(1)

gather_model_part = model.CreateModelPart("GatherPart",1)
# gather mesh 0 to process 0
gather_model_part_util = GatherModelPartUtility(0, model_part, 0, gather_model_part)

def PrintPressureOnNodes(model_part):
    for node in model_part.Nodes:
        if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
            print("Local Node[", node.Id, "]: PRESSURE = ", node.GetSolutionStepValue(PRESSURE))
        else:
            print("Ghost Node[", node.Id, "]: PRESSURE = ", node.GetSolutionStepValue(PRESSURE))
    sys.stdout.flush()

# first we assign some values on the local model part
for node in model_part.Nodes:
    if node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank:
        node.SetSolutionStepValue(PRESSURE, math.exp(-float(node.Id)))

if mpi.rank == 0:
    print("\nNew pressure values written to local model part")
    print("==========================================")
    print("Local model part on process 0")
    print("==========================================")
    PrintPressureOnNodes(model_part)
mpi.world.barrier()
if mpi.rank == 1:
    print("==========================================")
    print("Local model part on process 1")
    print("==========================================")
    PrintPressureOnNodes(model_part)
mpi.world.barrier()
if mpi.rank == 0:
    print("==========================================")
    print("Gather model part on process 0")
    print("==========================================")
    PrintPressureOnNodes(gather_model_part)
    print("Calling GatherOnMaster(PRESSURE) ...")
mpi.world.barrier()
gather_model_part_util.GatherOnMaster(PRESSURE)
if mpi.rank == 0:
    print("==========================================")
    print("Gather model part on process 0")
    print("==========================================")
    PrintPressureOnNodes(gather_model_part)
mpi.world.barrier()


if mpi.rank == 0:
    for node in gather_model_part.Nodes:
        node.SetSolutionStepValue(PRESSURE, math.exp(float(node.Id)))

if mpi.rank == 0:
    print("\nNew pressure values written to gather model part")
    print("==========================================")
    print("Gather model part on process 0")
    print("==========================================")
    PrintPressureOnNodes(gather_model_part)
mpi.world.barrier()
if mpi.rank == 0:
    print("==========================================")
    print("Local model part on process 0")
    print("==========================================")
    PrintPressureOnNodes(model_part)
mpi.world.barrier()
if mpi.rank == 1:
    print("==========================================")
    print("Local model part on process 1")
    print("==========================================")
    PrintPressureOnNodes(model_part)
    print("Calling ScatterFromMaster(PRESSURE) ...")
mpi.world.barrier()
gather_model_part_util.ScatterFromMaster(PRESSURE)
if mpi.rank == 0:
    print("==========================================")
    print("Local model part on process 0")
    print("==========================================")
    PrintPressureOnNodes(model_part)
mpi.world.barrier()
if mpi.rank == 1:
    print("==========================================")
    print("Local model part on process 1")
    print("==========================================")
    PrintPressureOnNodes(model_part)
mpi.world.barrier()
