try:
  import boost.mpi as mpi
except ImportError:
  import mpi
  
##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################
##################################################################
## ATTENTION: here the order is important

## from now on the order is not anymore crucial
##################################################################
###``lllll
kratos_path = '../../../..'
#kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking

import sys
sys.path.append(kratos_path)
sys.path.append(kratos_benchmarking_path)
import benchmarking

##################################################################
##################################################################
 # importing kratos
from KratosMultiphysics import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MetisApplication import *


#defining a model part
model_part = ModelPart("FluidPart");  


## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosThermoMechanicalApplication import *
from KratosTrilinosApplication import *
from KratosMetisApplication import *
#KKKKKKKK
##########################################################
thermal_settings = ConvectionDiffusionSettings()
thermal_settings.SetDensityVariable(DENSITY)
thermal_settings.SetDiffusionVariable(CONDUCTIVITY)
thermal_settings.SetUnknownVariable(TEMPERATURE)
thermal_settings.SetVolumeSourceVariable(HEAT_FLUX)
thermal_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
thermal_settings.SetConvectionVariable(VELOCITY)
thermal_settings.SetMeshVelocityVariable(MESH_VELOCITY)
##########################################################

#importing the solver files
import trilinos_thermal_solver
trilinos_thermal_solver.AddVariables(model_part,thermal_settings)




#adding of Variables to Model Part should be here when the "very fix container will be ready"
#reading a model
input_file_name="square_convection"
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)

number_of_partitions = mpi.size #we set it equal to the number of processors
print "number_of_partitions", number_of_partitions
partitioner = MetisPartitioningProcess(model_part, gid_io, number_of_partitions, domain_size);
partitioner.Execute()
print "GetRank()",mpi.rank

#model_part_io_fluid = ModelPartIO(input_file_name)

#print "before performing the division"
#number_of_partitions = mpi.size #we set it equal to the number of processors
#if mpi.rank == 0 :
    #partitioner = MetisDivideInputToPartitionsProcess(model_part_io_fluid, number_of_partitions, domain_size);
    #partitioner.Execute()

#print "division performed"

#mpi.world.barrier()

#MPICommSetup = SetMPICommunicatorProcess(fluid_model_part)
#MPICommSetup.Execute()

#my_input_filename = input_file_name + "_" + str(mpi.rank)
#model_part_io_fluid = ModelPartIO(my_input_filename)
#model_part_io_fluid.ReadModelPart(fluid_model_part)

## This is to correct the condition partition
#ccc= ParallelFillCommunicator(fluid_model_part)
#ccc.Execute()

#normal_calculator = NormalCalculationUtils()
#normal_calculator.CalculateOnSimplex(fluid_model_part.Conditions,domain_size)
#fluid_model_part.GetCommunicator().AssembleCurrentData(NORMAL)
#mpi.world.barrier()


##write down the mesh
mesh_name = mpi.rank
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetCommunicator().LocalMesh());
gid_io.FinalizeMesh()



#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#importing the solver files
trilinos_thermal_solver.AddDofs(model_part,thermal_settings)
print "111111111111111111111111111111111111111111"

#creating a fluid solver object
solver = trilinos_thermal_solver.Solver(model_part,domain_size,thermal_settings)
solver.Initialize()

#assigning the fluid properties
conductivity = 1.0;
density = 900.0;
specific_heat = 2400.0;
for node in model_part.Nodes:
    node.SetSolutionStepValue(thermal_settings.GetDiffusionVariable(),0,conductivity);
    node.SetSolutionStepValue(DENSITY,0,density);
    node.SetSolutionStepValue(SPECIFIC_HEAT,0,specific_heat);

#assigning a rotational velocity field
for node in model_part.Nodes:
    if(node.Y > 0.499 ):
         node.SetSolutionStepValue(TEMPERATURE,0,500.0);
         node.Fix(TEMPERATURE)

for node in model_part.Nodes:
    node.SetSolutionStepValue(VELOCITY_Y,0,-0.10);
    node.Fix(VELOCITY_Y)  
print "222222222222222222222222222222222222222222222222222222"

#settings to be changed
nsteps = 5000
output_step = 1

Dt = 1.0;

out = 0

gid_io.InitializeResults(mesh_name,(model_part).GetCommunicator().LocalMesh())

for step in range(0,nsteps):
    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]    #solving the fluid problem
    if(step > 3):
        solver.Solve()

    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(SPECIFIC_HEAT,model_part.Nodes,time,0) 
        gid_io.WriteNodalResults(DENSITY,model_part.Nodes,time,0) 
        gid_io.WriteNodalResults(CONDUCTIVITY,model_part.Nodes,time,0)      
        gid_io.WriteNodalResults(MESH_VELOCITY,model_part.Nodes,time,0)            
        out = 0
    out = out + 1

gid_io.FinalizeResults()          
        

