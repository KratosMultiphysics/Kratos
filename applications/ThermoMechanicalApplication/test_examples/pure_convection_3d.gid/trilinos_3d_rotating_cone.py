##################################################################
##################################################################
#setting the domain size for the problem to be solved
try:
  import boost.mpi as mpi
except ImportError:
  import mpi
  
domain_size = 3

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_ThermoMechanicalApplication = True
applications_interface.Import_KratosTrilinosApplication = True
applications_interface.Import_KratosMetisApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)


#defining a model part
model_part = ModelPart("FluidPart");  


## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosThermoMechanicalApplication import *
from KratosTrilinosApplication import *
from KratosMetisApplication import *

##########################################################
convection_settings = ConvectionDiffusionSettings()
convection_settings.SetUnknownVariable(TEMPERATURE)
convection_settings.SetConvectionVariable(VELOCITY)
convection_settings.SetMeshVelocityVariable(MESH_VELOCITY)
##########################################################

#importing the solver files
import trilinos_pureconvection_solver
trilinos_pureconvection_solver.AddVariables(model_part,convection_settings)




#adding of Variables to Model Part should be here when the "very fix container will be ready"
#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("pure_convection_3d",gid_mode,multifile,deformed_mesh_flag, write_conditions)
#gid_io.ReadModelPart(model_part)

#read with old input format
number_of_partitions = mpi.size #we set it equal to the number of processors
print "number_of_partitions", number_of_partitions
partitioner = MetisPartitioningProcess(model_part, gid_io, number_of_partitions, domain_size);
partitioner.Execute()
print "GetRank()",mpi.rank

model_part.Conditions.clear()

#read with new input format
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

mesh_name = mpi.rank
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
#print model_part


#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#importing the solver files
trilinos_pureconvection_solver.AddDofs(model_part,convection_settings)

#creating a fluid solver object
solver = trilinos_pureconvection_solver.Solver(model_part,domain_size,convection_settings)
#pDiagPrecond = DiagonalPreconditioner()
#solver.linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
solver.Initialize()
#solver.SetEchoLevel(3)

#assigning a rotational velocity field
vel = Vector(3);
for node in model_part.Nodes:
      vel[0] = 0.0
      vel[1] = -node.Z
      vel[2] = node.Y
  #else:
      #vel[0] = 0.0
      #vel[1] = 0.0
      #vel[2] = 0.0    
##   if(node.X**2 + node.Y**2 > 0.24):
##        vel[0] = 0.0
##        vel[1] = 0.0        
      node.SetSolutionStepValue(convection_settings.GetConvectionVariable(),0,vel);

#assigning a cone shaped temperature distribution
xc = 0.0
yc = 0.0
zc = 1.00/6.00
sigma = 0.2
import math
for node in model_part.Nodes:
    X1 = (node.X - xc) / sigma
    X2 = (node.Y - yc) / sigma
    X3 = (node.Z - zc) / sigma    
    if( (X1**2 + X2**2 + X3**2) <= 1.00  ):
        temp = 0.25*(1.00+math.cos(math.pi*X1) )*(1.00+math.cos(math.pi*X2) )*(1.00+math.cos(math.pi*X3) )
        node.SetSolutionStepValue(convection_settings.GetUnknownVariable(),0,temp)
    else:
        node.SetSolutionStepValue(convection_settings.GetUnknownVariable(),0,0.0)
    if(node.Z == 0.5 or node.Z == -0.5 or node.Y== 0.5 or node.Y == -0.5):
        node.SetSolutionStepValue(convection_settings.GetUnknownVariable(),0,0.0)
        node.Fix(convection_settings.GetUnknownVariable())
 
print "222222222222222222222222222222222222222222222222222222"

#settings to be changed
nsteps = 500
output_step = 1

Dt = 2.00*math.pi/200.0;

out = 0

gid_io.InitializeResults(mesh_name,(model_part).GetCommunicator().LocalMesh())

for step in range(0,nsteps):

    time = Dt*step
    model_part.CloneTimeStep(time)

    if(mpi.rank == 0):
      print "time = ",time, " Dt = ",Dt
    #print model_part.ProcessInfo()[TIME]    #solving the fluid problem
    if(step > 3):
        solver.Solve()

    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(convection_settings.GetUnknownVariable(),model_part.Nodes,time,0)
        gid_io.WriteNodalResults(convection_settings.GetConvectionVariable(),model_part.Nodes,time,0)
        gid_io.WriteNodalResults(convection_settings.GetMeshVelocityVariable(),model_part.Nodes,time,0)   
        out = 0
    out = out + 1

gid_io.FinalizeResults()          
          
        

