import mpi
##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications
kratos_benchmarking_path = '../../../../benchmarking' ##kratos_root/benchmarking
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(kratos_benchmarking_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_KratosTrilinosApplication = True
applications_interface.Import_KratosMetisApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosIncompressibleFluidApplication import *
from KratosTrilinosApplication import *
from KratosMetisApplication import *
import benchmarking


#defining a model part
model_part = ModelPart("FluidPart");  

#providing the variable list to the model part
import trilinos_fs_fluid_solver
trilinos_fs_fluid_solver.AddVariables(model_part)
model_part.AddNodalSolutionStepVariable(NORMAL)

#reading a model
gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly

input_file_name = "cavity2D"
gid_io = GidIO(input_file_name, gid_mode_flag, use_multifile, deformed_print_flag, write_conditions)


model_part_io_fluid = ModelPartIO(input_file_name)

print "before performing the division"
number_of_partitions = mpi.size #we set it equal to the number of processors
if mpi.rank == 0 :
    partitioner = MetisDivideInputToPartitionsProcess(model_part_io_fluid, number_of_partitions, domain_size);
    partitioner.Execute()

"print division performed"

mpi.world.barrier()

MPICommSetup = SetMPICommunicatorProcess(model_part)
MPICommSetup.Execute()

my_input_filename = input_file_name + "_" + str(mpi.rank)
model_part_io_fluid = ModelPartIO(my_input_filename)
model_part_io_fluid.ReadModelPart(model_part)

##number_of_partitions = mpi.size #we set it equal to the number of processors
##print "number_of_partitions", number_of_partitions
##partitioner = MetisPartitioningProcess(model_part, gid_io, number_of_partitions, domain_size);
##partitioner.Execute()
##print "GetRank()",GetRank()




print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#adding dofs
trilinos_fs_fluid_solver.AddDofs(model_part)
for node in model_part.Nodes:
    node.AddDof(DISPLACEMENT_X)
    node.AddDof(DISPLACEMENT_Y)
    node.AddDof(DISPLACEMENT_Z)


#creating a fluid solver object
fluid_solver = trilinos_fs_fluid_solver.IncompressibleFluidSolver(model_part,domain_size)
#fluid_solver = flexible_trilinos_fs_fluid_solver.IncompressibleFluidSolver(model_part,domain_size)
fluid_solver.laplacian_form =1;
fluid_solver.predictor_corrector = True;
fluid_solver.vel_toll = 1e-3
fluid_solver.time_order = 2
fluid_solver.ReformDofAtEachIteration = False
fluid_solver.echo_level = 0

##pILUPrecond = ILU0Preconditioner() 
##pDiagPrecond = DiagonalPreconditioner() 
##fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)
##fluid_solver.pressure_linear_solver = SkylineLUFactorizationSolver();

##for node in model_part.Nodes:
##    if(node.X < 0.001 and node.Y<0.001):
##        node.Fix(PRESSURE)

import PressureMultiLevelSolver
fluid_solver.pressure_linear_solver =  PressureMultiLevelSolver.MultilevelLinearSolver(1e-4,1000)


fluid_solver.Initialize()

#settings to be changed
Re = 100.0
nsteps = 100
output_step = 2

Dt = 1000.0/Re
if(Dt > 0.1):
    Dt = 0.1
out = 0
for node in model_part.Nodes:
    node.SetSolutionStepValue(DENSITY,0,1.0)
    node.SetSolutionStepValue(VISCOSITY,0,1.0/Re)



for elem in model_part.Elements:
    if (elem.Id < 200):
        elem.SetValue(SPLIT_ELEMENT,True)
#for elem in model_part.Elements:
#    elem.SetValue(SPLIT_ELEMENT,True)

ccc= ParallelFillCommunicator(model_part)
ccc.Execute()

print "*******************************************************************"

Comm = CreateCommunicator()
mesh_utility = TrilinosRefineMesh(model_part,Comm)
refine_on_reference = False
interpolate_internal_variables = False
mesh_utility.Local_Refine_Mesh(refine_on_reference,interpolate_internal_variables,domain_size)


print "-----------------------------------------------------------------"
ccc= ParallelFillCommunicator(model_part)
ccc.Execute()
##ccc.PrintDebugInfo()



mesh_utility.PrintDebugInfo()
aaa = BodyNormalCalculationUtils()
aaa.CalculateBodyNormals(model_part,2)

mpi.world.barrier()

print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
#Comm = CreateCommunicator()
#mesh_utility = TrilinosRefineMesh(model_part,Comm)
#refine_on_reference = False
#interpolate_internal_variables = False
#mesh_utility.Local_Refine_Mesh(refine_on_reference,interpolate_internal_variables)



##zero = Vector(3);
##zero[0] = 0.0;
##zero[1] = 0.0;
##zero[2] = 0.0;

mesh_name = mpi.rank
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
gid_io.Flush()


gid_io.InitializeResults(mesh_name , model_part.GetMesh())
for step in range(0,nsteps):

    time = Dt*step
    model_part.CloneTimeStep(time)

    if(mpi.rank == 0):
        print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 4):
        fluid_solver.Solve()

##        for node in model_part.Nodes:
##            node.SetSolutionStepValue(PRESS_PROJ,0,zero);


    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(PARTITION_INDEX,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(NORMAL,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        out = 0
    out = out + 1
    
gid_io.FinalizeResults()

          
        

