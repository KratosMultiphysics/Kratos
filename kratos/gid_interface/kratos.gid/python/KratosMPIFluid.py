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
                
    if(ProjectParameters.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras"):
	gid_io.WriteNodalResults(VISCOSITY,nodes,time,0)
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
from KratosMultiphysics.MeshingApplication import *

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

##################################################################
##################################################################
##importing the solvers needed
SolverSettings = ProjectParameters.FluidSolverConfiguration
solver_constructor = __import__(SolverSettings.solver_type)

##################################################################
##################################################################
##importing variables
solver_constructor.AddVariables( fluid_model_part, SolverSettings)
  
  
  
if(ProjectParameters.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras"):
  fluid_model_part.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY);
  fluid_model_part.AddNodalSolutionStepVariable(MOLECULAR_VISCOSITY);
  fluid_model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)

  fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)
  
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
    partitioner = MetisDivideHeterogeneousInputProcess(model_part_io_fluid, number_of_partitions , domain_size, 1)
    partitioner.Execute()

mpi.world.barrier()

MPICommSetup = SetMPICommunicatorProcess(fluid_model_part)
MPICommSetup.Execute()

(ParallelFillCommunicator(fluid_model_part)).Execute()

my_input_filename = input_file_name + "_" + str(mpi.rank)
model_part_io_fluid = ModelPartIO(my_input_filename)
model_part_io_fluid.ReadModelPart(fluid_model_part)
########################################################################

Comm = CreateCommunicator()
       

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)


##adding dofs
solver_constructor.AddDofs( fluid_model_part, SolverSettings)
    
if(ProjectParameters.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras"):
    for node in fluid_model_part.Nodes:
       node.AddDof(TURBULENT_VISCOSITY)
       
# If Lalplacian form = 2, free all pressure Dofs
#laplacian_form = ProjectParameters.FluidSolverConfiguration.laplacian_form 
#if(laplacian_form >= 2):
    #for node in fluid_model_part.Nodes:
        #node.Free(PRESSURE)

#copy Y_WALL
for node in fluid_model_part.Nodes:
    y = node.GetSolutionStepValue(Y_WALL,0)
    node.SetValue(Y_WALL,y)

##################################################################
##################################################################
##Creating the fluid solver
fluid_solver = solver_constructor.CreateSolver( fluid_model_part, SolverSettings)

##activate turbulence model
if(ProjectParameters.FluidSolverConfiguration.TurbulenceModel == "Smagorinsky-Lilly"):
    fluid_solver.ActivateSmagorinsky(ProjectParameters.SmagorinskyConstant)
elif(ProjectParameters.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras"):
    ##apply the initial turbulent viscosity on all of the nodes
    turb_visc = ProjectParameters.TurbulentViscosity
    for node in fluid_model_part.Nodes:
      node.SetSolutionStepValue(TURBULENT_VISCOSITY,0,turb_visc);
      visc = node.GetSolutionStepValue(VISCOSITY)
      node.SetSolutionStepValue(MOLECULAR_VISCOSITY,0,visc);
      if(node.IsFixed(VELOCITY_X)):
	  node.Fix(TURBULENT_VISCOSITY)
	  
	  
    ##select nodes on the wall
    wall_nodes = []
    for i in ProjectParameters.SA_wall_group_ids:
       nodes = fluid_model_part.GetNodes(i) ##get the nodes of the wall for SA.
       for node in nodes:
	  wall_nodes.append(node)
	  node.SetSolutionStepValue(TURBULENT_VISCOSITY,0,0.0);
	  node.Fix(TURBULENT_VISCOSITY)
	  
    DES = False
    fluid_solver.ActivateSpalartAllmaras(wall_nodes,DES)

fluid_solver.Initialize()
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
#######################################33
#define the drag computation list   
drag_list = define_output.DefineDragList()
drag_file_output_list = []

if(mpi.rank == 0): 
  for it in drag_list:
      f = open(it[1],'w')
      drag_file_output_list.append(f)
      tmp = "#Drag for group " + it[1] + "\n"
      f.write(tmp)
      tmp = "time RX RY RZ"
      f.write(tmp)
      f.flush()
    
print drag_file_output_list
    
def PrintDrag(drag_list,drag_file_output_list,fluid_model_part,time):
    i = 0
    for it in drag_list:
      print it[0]
      nodes = fluid_model_part.GetNodes(it[0])
      dx = 0.0;
      dy = 0.0;
      dz = 0.0;
      
      for node in nodes:
	  if(node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank):
	      reaction = node.GetSolutionStepValue(REACTION,0)
	      dx += reaction[0]
	      dy += reaction[1]
	      dz += reaction[2]
	  
      auxx = mpi.gather(mpi.world,dx,0)
      auxy = mpi.gather(mpi.world,dy,0)
      auxz = mpi.gather(mpi.world,dz,0)
      print auxx
      
      rx = 0.0;
      ry = 0.0;
      rz = 0.0;
      for k in auxx:
	rx += k
      for k in auxy:
	ry += k
      for k in auxz:
	rz += k
	
      if(mpi.rank == 0):
	  output = str(time) + " " + str(rx) +  " " + str(ry) +  " " + str(rz) + "\n"
	  #print drag_file_output_list[i]
	  #print output
	  drag_file_output_list[i].write(output)
	  drag_file_output_list[i].flush()
	  
      i = i+1    
    
#######################################33
#preparing output of point graphs
import point_graph_printer

output_nodes_list = define_output.DefineOutputPoints()
graph_printer = point_graph_printer.PrintGraphPrinter(output_nodes_list, fluid_model_part, variables_dictionary, domain_size)


# Stepping and time settings
Dt = ProjectParameters.Dt 
Nsteps  = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0
zero_vector = Vector(3)
zero_vector[0] = 0.0
zero_vector[1] = 0.0
zero_vector[2] = 0.0

while(time <= final_time):

    time = time + Dt
    step = step + 1
    fluid_model_part.CloneTimeStep(time)

    print "STEP = ", step
    print "TIME = ", time

    if(step >= 3):
        fluid_solver.Solve()
        
        if(step < 4):
	    for k in range(0,ProjectParameters.divergence_cleareance_step):
		if(mpi.rank == 0):
		  print "DOING DIVERGENCE CLEAREANCE"
		buffer_size = fluid_model_part.GetBufferSize()
		for i in range(0,buffer_size):
		  for node in fluid_model_part.Nodes:
		    vel = node.GetSolutionStepValue(VELOCITY)
		    node.SetSolutionStepValue(VELOCITY,i,vel)
		    node.SetSolutionStepValue(PRESSURE,i,0.0)		    
		  if(SolverSettings.solver_type == "trilinos_vms_monolithic_solver"):
		    for node in fluid_model_part.Nodes:
		      node.SetSolutionStepValue(ACCELERATION,i,zero_vector)
		  if(ProjectParameters.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras"):
		    for node in fluid_model_part.Nodes:
		      visc = node.GetSolutionStepValue(VISCOSITY)
		      node.SetSolutionStepValue(VISCOSITY,i,visc)
		      
		fluid_solver.Solve()

        
        graph_printer.PrintGraphs(time)

    if(output_time <= out):
	if(ProjectParameters.VolumeOutput == True):
	    
	    local_nodes = fluid_model_part.GetCommunicator().LocalMesh().Nodes
	    PrintResults(local_nodes)
	    PrintDrag(drag_list,drag_file_output_list,fluid_model_part,time)
	    gid_io.Flush()
	    out = 0
	else:
	    cut_model_part.CloneTimeStep(time)
	    Cut_App.UpdateCutData(cut_model_part,fluid_model_part)
	    	    
	    PrintResults(cut_model_part.Nodes)
	    PrintDrag(drag_list,drag_file_output_list,fluid_model_part,time)
	    
	    out = 0	    

    out = out + Dt

gid_io.FinalizeResults()
    
          
        
