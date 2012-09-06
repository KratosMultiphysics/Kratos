#################################################################
##################################################################
#import the configuration data as read from the GiD
import ProjectParameters
import define_output

def PrintResults(nodes):
    if(mpi.rank == 0):
      print "Writing results. Please run Gid for viewing results of analysis."
    for variable_name in ProjectParameters.nodal_results:
        gid_io.WriteNodalResults(variables_dictionary[variable_name],nodes,time,0)
    for variable_name in ProjectParameters.gauss_points_results:
        gid_io.PrintOnGaussPoints(variables_dictionary[variable_name],model_part,time)
    gid_io.Flush()
    mpi.world.barrier()


##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

##################################################################
##################################################################
import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MetisApplication import *

## defining variables to be used

variables_dictionary = {"PRESSURE" : PRESSURE,
                        "VELOCITY" : VELOCITY,
                        "REACTION" : REACTION,
                        "DISTANCE" : DISTANCE,}

#defining a model part for the fluid 
fluid_model_part = ModelPart("FluidPart");

if "REACTION" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(REACTION)
if "DISTANCE" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

#############################################
##importing the solvers needed
SolverType = ProjectParameters.SolverType
if(SolverType == "FractionalStep"):
    import trilinos_vms_fs_fluid_solver as solver
    solver.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    import vms_monolithic_solver as solver
    solver.AddVariables(fluid_model_part)
else:
    raise NameError("solver type not supported: options are FractionalStep  - monolithic_solver_eulerian")

#introducing input file name
input_file_name = ProjectParameters.problem_name

#reading the fluid part

# initialize GiD  I/O
if ProjectParameters.GiDPostMode == "Binary":
    gid_mode = GiDPostMode.GiD_PostBinary
elif ProjectParameters.GiDPostMode == "Ascii":
    gid_mode = GiDPostMode.GiD_PostAscii
elif ProjectParameters.GiDPostMode == "AsciiZipped":
    gid_mode = GiDPostMode.GiD_PostAsciiZipped
else:
    print "Unknown GiD post mode, assuming Binary"
    gid_mode = GiDPostMode.GiD_PostBinary

if ProjectParameters.GiDWriteMeshFlag == True:
    deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
else:
    deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed

if(ProjectParameters.VolumeOutput == True):
    if ProjectParameters.GiDWriteConditionsFlag == True:
	write_conditions = WriteConditionsFlag.WriteConditions
    else:
	write_conditions = WriteConditionsFlag.WriteElementsOnly
else:
    write_conditions = WriteConditionsFlag.WriteConditions

multifile = MultiFileFlag.MultipleFiles
    
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)

########################## do parallel reading ######################
number_of_partitions = mpi.size #we set it equal to the number of processors
if mpi.rank == 0 :
    partitioner = MetisDivideInputToPartitionsProcess(model_part_io_fluid, number_of_partitions, domain_size);
    partitioner.Execute()

mpi.world.barrier()

MPICommSetup = SetMPICommunicatorProcess(fluid_model_part)
MPICommSetup.Execute()

my_input_filename = input_file_name + "_" + str(mpi.rank)
model_part_io_fluid = ModelPartIO(my_input_filename)
model_part_io_fluid.ReadModelPart(fluid_model_part)
########################################################################

Comm = CreateCommunicator()
       

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
if SolverType == "FractionalStep":
    fluid_model_part.SetBufferSize(3)
else:
    fluid_model_part.SetBufferSize(2)

##adding dofs
if(SolverType == "FractionalStep"):
    solver.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    solver.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian_compressible"):
    solver.AddDofs(fluid_model_part)

# If Lalplacian form = 2, free all pressure Dofs
laplacian_form = ProjectParameters.laplacian_form 
if(laplacian_form >= 2):
    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)



dynamic_tau = ProjectParameters.use_dt_in_stabilization
oss_switch = ProjectParameters.use_orthogonal_subscales
#creating the solvers
#fluid solver
if(SolverType == "FractionalStep"):
    fluid_solver = solver.IncompressibleFluidSolver(fluid_model_part,domain_size)
    fluid_solver.max_val_its = ProjectParameters.max_vel_its
    fluid_solver.max_press_its = ProjectParameters.max_press_its
    fluid_solver.predictor_corrector = ProjectParameters.predictor_corrector            
    fluid_solver.vel_toll = ProjectParameters.velocity_relative_tolerance
    fluid_solver.press_toll = ProjectParameters.pressure_relative_tolerance
    fluid_solver.dynamic_tau = float(dynamic_tau)
    fluid_solver.compute_reactions = ProjectParameters.Calculate_reactions
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian"): 
    fluid_solver = solver.MonolithicSolver(fluid_model_part,domain_size)
    fluid_solver.oss_switch = int(oss_switch)
    fluid_solver.dynamic_tau = float(dynamic_tau)
    fluid_solver.rel_vel_tol = ProjectParameters.velocity_relative_tolerance
    fluid_solver.abs_vel_tol = ProjectParameters.velocity_absolute_tolerance
    fluid_solver.rel_pres_tol = ProjectParameters.pressure_relative_tolerance
    fluid_solver.abs_pres_tol = ProjectParameters.pressure_absolute_tolerance
    fluid_solver.max_iter = ProjectParameters.max_iterations
    fluid_solver.compute_reactions = ProjectParameters.Calculate_reactions
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian_compressible"): 
    print "monolithic_solver_eulerian_compressible is not available in mpi"
    err


print "fluid solver created"



cut_model_part = ModelPart("CutPart");
if(ProjectParameters.VolumeOutput == True):
    mesh_name = mpi.rank
    gid_io.InitializeMesh( mesh_name )
    gid_io.WriteMesh( fluid_model_part.GetMesh() )
    gid_io.FinalizeMesh()

    gid_io.InitializeResults(mesh_name,(fluid_model_part).GetMesh())
    gid_io.Flush()
    Multifile = True
    
    # Write .post.list file (GiD postprocess list)
    if(mpi.rank == 0):
	f = open(ProjectParameters.problem_name+'.post.lst','w')
	f.write('Merge\n')
	if ProjectParameters.GiDPostMode == "Binary":
	    nproc = mpi.size
	    for i in range(0,nproc):
	      f.write(ProjectParameters.problem_name+'_'+str(i)+'.post.bin\n')  
	elif ProjectParameters.GiDPostMode == "Ascii":
	    nproc = mpi.size
	    for i in range(0,nproc):
	      f.write(ProjectParameters.problem_name+'_'+str(i)+'.post.msh\n')
	f.close()
	
else:
    #generate the cuts
    Cut_App = TrilinosCuttingApplication(Comm);
    Cut_App.FindSmallestEdge(fluid_model_part)
    
    cut_number = 1

    cut_list = define_output.DefineCutPlanes()
    print cut_list
    print "***************** i am rank ",mpi.rank
    
    for item in cut_list:
       print item
       n = Vector( item[0] )
       p = Vector( item[1] )
       
       Cut_App.GenerateCut(fluid_model_part,cut_model_part,n,p, cut_number , 0.01)
       cut_number = cut_number + 1
       print "generated cut number =",cut_number
      
    if(len(cut_model_part.Conditions) != 0):
	Cut_App.AddSkinConditions(fluid_model_part,cut_model_part, cut_number)
	cut_number += 1      
    
    ###mesh to be printed (single mesh case)
    mesh_name = mpi.rank
    gid_io.InitializeMesh( mesh_name )
    gid_io.WriteMesh( cut_model_part.GetMesh() )
    gid_io.FinalizeMesh()

    gid_io.InitializeResults(mesh_name,(cut_model_part).GetMesh())
    gid_io.Flush()
    Multifile = True
    
    # Write .post.list file (GiD postprocess list)
    if(mpi.rank == 0):
	f = open(ProjectParameters.problem_name+'.post.lst','w')
	f.write('Merge\n')
	if ProjectParameters.GiDPostMode == "Binary":
	    nproc = mpi.size
	    for i in range(0,nproc):
	      f.write(ProjectParameters.problem_name+'_'+str(i)+'.post.bin\n')  
	elif ProjectParameters.GiDPostMode == "Ascii":
	    nproc = mpi.size
	    for i in range(0,nproc):
	      f.write(ProjectParameters.problem_name+'_'+str(i)+'.post.msh\n')
	f.close()  
    
    
    
#######################################33
#preparing output of point graphs
import point_graph_printer

output_nodes_list = define_output.DefineOutputPoints()
graph_printer = point_graph_printer.PrintGraphPrinter(output_nodes_list, fluid_model_part, variables_dictionary, domain_size)


# Stepping and time settings
Dt = ProjectParameters.Dt 
full_Dt = Dt 
initial_Dt = 0.001 * full_Dt #0.05 #0.01
Nsteps  = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0

while(time <= final_time):

    if(step < 3):
        Dt = initial_Dt
    else:
        Dt = full_Dt
        
    time = time + Dt
    step = step + 1
    fluid_model_part.CloneTimeStep(time)

    print "STEP = ", step
    print "TIME = ", time

    if(step >= 3):
        fluid_solver.Solve()
        
        graph_printer.PrintGraphs(time)

    if(output_time <= out):
	if(ProjectParameters.VolumeOutput == True):
	    
	    local_nodes = fluid_model_part.GetCommunicator().LocalMesh().Nodes
	    PrintResults(local_nodes)
	    gid_io.Flush()
	    out = 0
	else:
	    cut_model_part.CloneTimeStep(time)
	    Cut_App.UpdateCutData(cut_model_part,fluid_model_part)
	    	    
	    PrintResults(cut_model_part.Nodes)
	    
	    out = 0	    

    out = out + Dt

if Multifile:
    f.close()
else:
    gid_io.FinalizeResults()
    
          
        
