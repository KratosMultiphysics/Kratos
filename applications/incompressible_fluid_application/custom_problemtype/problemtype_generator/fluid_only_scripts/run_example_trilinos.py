import mpi
import fluid_only_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = fluid_only_var.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
import sys
sys.path.append(fluid_only_var.kratos_path)


#importing Kratos main library
from KratosMultiphysics import *

## from now on the order is not anymore crucial
##################################################################
##################################################################
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MetisApplication import *

#defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart");  

#############################################
##importing the solvers needed
SolverType = fluid_only_var.SolverType
if(SolverType == "FractionalStep"):
    import trilinos_fs_fluid_solver
    trilinos_fs_fluid_solver.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    import trilinos_monolithic_solver_eulerian
    trilinos_monolithic_solver_eulerian.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian_compressible"):
    import trilinos_monolithic_solver_eulerian_compressible
    trilinos_monolithic_solver_eulerian_compressible.AddVariables(fluid_model_part)
else:
    raise "solver type not supported: options are FractionalStep - Monolithic"

#introducing input file name
input_file_name = fluid_only_var.problem_name

#reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)

model_part_io_fluid = ModelPartIO(input_file_name)

#do the reading through the partitioner
number_of_partitions = mpi.size #we set it equal to the number of processors
partitioner = MetisPartitioningProcess(fluid_model_part, model_part_io_fluid, number_of_partitions, domain_size);
partitioner.Execute()

#write down the mesh
mesh_name = mpi.rank
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((fluid_model_part).GetMesh());
gid_io.FinalizeMesh()


##adding dofs
if(SolverType == "FractionalStep"):
    #setting up the buffer size: SHOULD BE DONE AFTER READING!!!
    fluid_model_part.SetBufferSize(3)
    trilinos_fs_fluid_solver.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    #setting up the buffer size: SHOULD BE DONE AFTER READING!!!
    fluid_model_part.SetBufferSize(2)
    trilinos_monolithic_solver_eulerian.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian_compressible"):
    #setting up the buffer size: SHOULD BE DONE AFTER READING!!!
    fluid_model_part.SetBufferSize(2)
    trilinos_monolithic_solver_eulerian_compressible.AddDofs(fluid_model_part)

#########select here the laplacian form!!!!!!!!!!!!!!!!!
laplacian_form = fluid_only_var.laplacian_form 
if(laplacian_form >= 2):
    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

##check to ensure that no node has zero density or pressure
for node in fluid_model_part.Nodes:
    if(node.GetSolutionStepValue(DENSITY) == 0.0):
        print "node ",node.Id," has zero density!"
        raise 'node with zero density found'
    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print "node ",node.Id," has zero viscosity!"
        raise 'node with zero VISCOSITY found'

#creating the solvers
#fluid solver
if(SolverType == "FractionalStep"):
    fluid_solver = trilinos_fs_fluid_solver.IncompressibleFluidSolver(fluid_model_part,domain_size)
    fluid_solver.laplacian_form = laplacian_form; #standard laplacian form
    fluid_solver.predictor_corrector = fluid_only_var.predictor_corrector
    fluid_solver.max_press_its = fluid_only_var.max_press_its
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian"): 
    fluid_solver = trilinos_monolithic_solver_eulerian.MonolithicSolver(fluid_model_part,domain_size)
    oss_swith = fluid_only_var.use_oss
    dynamic_tau = fluid_only_var.dynamic_tau
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_swith);				
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian_compressible"): 
    fluid_solver = trilinos_monolithic_solver_eulerian_compressible.MonolithicSolver(fluid_model_part,domain_size)
    oss_swith = fluid_only_var.use_oss
    dynamic_tau = fluid_only_var.dynamic_tau
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_swith);				
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.Initialize()


print "fluid solver created"

#settings to be changed
Dt = fluid_only_var.Dt 
full_Dt = Dt 
initial_Dt = 0.001 * full_Dt #0.05 #0.01
final_time = fluid_only_var.max_time
output_step = fluid_only_var.output_step

out = 0




gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh())


time = 0.0
step = 0
while(time < final_time):

    if(step < 5):
        Dt = initial_Dt
    else:
        Dt = full_Dt
        
    time = time + Dt
    fluid_model_part.CloneTimeStep(time)

    if(step >= 3):
        fluid_solver.Solve()

    if(out == output_step):
        gid_io.WriteNodalResults(PARTITION_INDEX,fluid_model_part.Nodes,time,0)
        if(SolverType == "FractionalStep"):
            gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
        else:
            gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(AIR_PRESSURE,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(WATER_PRESSURE,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(MESH_VELOCITY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(IS_STRUCTURE,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(IS_BOUNDARY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(IS_POROUS,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(IS_FREE_SURFACE,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(ADVPROJ,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(DIVPROJ,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(DENSITY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(DENSITY_AIR,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(VISCOSITY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(SOUND_VELOCITY,fluid_model_part.Nodes,time,0)
            gid_io.WriteNodalResults(AIR_SOUND_VELOCITY,fluid_model_part.Nodes,time,0)

        out = 0

    out = out + 1
    step = step + 1
      
gid_io.FinalizeResults()
          
        

