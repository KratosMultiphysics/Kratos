#################################################################
##################################################################
#import the configuration data as read from the GiD
import ProjectParameters
import define_output

def PrintResults(model_part):
    print "Writing results. Please run Gid for viewing results of analysis."
    for variable_name in ProjectParameters.nodal_results:
        gid_io.WriteNodalResults(variables_dictionary[variable_name],model_part.Nodes,time_t,0)
    for variable_name in ProjectParameters.gauss_points_results:
        gid_io.PrintOnGaussPoints(variables_dictionary[variable_name],model_part,time_t)

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

##################################################################
##################################################################
# functions needed for Mesh movement

def SelectNodes(moving_nodes,model_part):
    for node in model_part.Nodes:
        if(node.Y > 0.49 and node.Y < 0.51):
            if(node.X > 0.1 and node.X < 0.9 ):
                moving_nodes.append(node)
                node.Fix(DISPLACEMENT_X);
                node.Fix(DISPLACEMENT_Y);
                node.Fix(DISPLACEMENT_Z);
                
def ApplyDisplacementConditions(model_part):
    for node in model_part.Nodes:
        if(node.Y > 0.9999 or node.Y < 0.0001 or node.X < 0.0001 or node.X > 0.999):
            node.Fix(DISPLACEMENT_X);
            node.Fix(DISPLACEMENT_Y);
            node.Fix(DISPLACEMENT_Z);
                
def MoveNodes(moving_nodes,time_t):
    disp = Vector(3)
    disp[0] = 0.0; disp[2] = 0.0
    disp[1] = 0.05*math.sin(time_t*1.0)
    print disp
    for node in moving_nodes:
        node.SetSolutionStepValue(DISPLACEMENT,0,disp);

##################################################################
##################################################################
import sys
import math
import time
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

## defining variables to be used

variables_dictionary = {"PRESSURE" : PRESSURE,
                        "VELOCITY" : VELOCITY,
                        "REACTION" : REACTION,
                        "DISTANCE" : DISTANCE}

#defining a model part for the fluid 
fluid_model_part = ModelPart("FluidPart");

# Adding variables
fluid_model_part.AddNodalSolutionStepVariable(PRESSURE)
fluid_model_part.AddNodalSolutionStepVariable(VELOCITY)
fluid_model_part.AddNodalSolutionStepVariable(REACTION)
fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)


#############################################
##importing the solvers needed
import vms_fractional_step_solver as fluid_solver
import trilinos_mesh_solver as trilinos_mesh_solver

# Add Varialbes to solver
fluid_solver.AddVariables(fluid_model_part)
trilinos_mesh_solver.AddVariables(fluid_model_part)

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

if ProjectParameters.GiDMultiFileFlag == "Single":
    multifile = MultiFileFlag.SingleFile
elif ProjectParameters.GiDMultiFileFlag == "Multiples":
    multifile = MultiFileFlag.MultipleFiles
else:
    print "Unknown GiD multiple file mode, assuming Single"
    multifile = MultiFileFlag.SingleFile
    
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

##adding dofs
fluid_solver.AddDofs(fluid_model_part)
trilinos_mesh_solver.AddDofs(fluid_model_part)
       
# If Lalplacian form = 2, free all pressure Dofs
laplacian_form = ProjectParameters.laplacian_form 
if(laplacian_form >= 2):
    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

# decoupled solver schemes: need solver for pressure and solver fol velocity

if ProjectParameters.SolverType in ["FractionalStep"]:
    # Velocity preconditioner
    try:
        if ProjectParameters.Velocity_Preconditioner_type == 'Diagonal':
            vel_precond = DiagonalPreconditioner()
        elif ProjectParameters.Velocity_Preconditioner_type == 'ILU0':
            vel_precond = ILU0Preconditioner()
    except AttributeError:
        # If we are using a direct solver Velocity_Preconditioner_type will not
        # exist. This is not an error.
        vel_precond = None

    # Velocity solver
    if ProjectParameters.Velocity_Linear_Solver == "Conjugate gradient":
        if ProjectParameters.Velocity_Preconditioner_type == 'None':
            velocity_linear_solver =  CGSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration)
        else:
            velocity_linear_solver =  CGSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration,vel_precond)
    elif ProjectParameters.Velocity_Linear_Solver == "BiConjugate gradient stabilized":
        if ProjectParameters.Velocity_Preconditioner_type == 'None':
            velocity_linear_solver =  BICGSTABSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration)
        else:
            velocity_linear_solver =  BICGSTABSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration,vel_precond)
    elif ProjectParameters.Velocity_Linear_Solver == "GMRES":
        if ProjectParameters.Velocity_Preconditioner_type == 'None':
            velocity_linear_solver =  GMRESSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration)
        else:
            velocity_linear_solver =  GMRESSolver(ProjectParameters.Velocity_Iterative_Tolerance, ProjectParameters.Velocity_Solver_Max_Iteration,vel_precond)
    elif ProjectParameters.Velocity_Linear_Solver == "Skyline LU factorization":
        velocity_linear_solver = SkylineLUFactorizationSolver()
    elif ProjectParameters.Velocity_Linear_Solver == "Super LU":
        velocity_linear_solver = SuperLUSolver()
    elif ProjectParameters.Velocity_Linear_Solver == "Parallel MKL Pardiso":
        from KratosMultiphysics.MKLSolversApplication import MKLPardisoSolver
        velocity_linear_solver = MKLPardisoSolver()

    # Pressure preconditioner
    try:
        if ProjectParameters.Pressure_Preconditioner_type == 'Diagonal':
            press_precond = DiagonalPreconditioner()
        elif ProjectParameters.Pressure_Preconditioner_type == 'ILU0':
            press_precond = ILU0Preconditioner()
    except AttributeError:
        # If we are using a direct solver Pressure_Preconditioner_type will not
        # exist. This is not an error.
        press_precond = None

    # Pressure solver
    if ProjectParameters.Pressure_Linear_Solver == "Conjugate gradient":
        if ProjectParameters.Pressure_Preconditioner_type == 'None':
            pressure_linear_solver =  CGSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration)
        else:
            pressure_linear_solver =  CGSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration,press_precond)
    elif ProjectParameters.Pressure_Linear_Solver == "BiConjugate gradient stabilized":
        if ProjectParameters.Pressure_Preconditioner_type == 'None':
            pressure_linear_solver =  BICGSTABSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration)
        else:
            pressure_linear_solver =  BICGSTABSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration,press_precond)
    elif ProjectParameters.Pressure_Linear_Solver == "GMRES":
        if ProjectParameters.Pressure_Preconditioner_type == 'None':
            pressure_linear_solver =  GMRESSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration)
        else:
            pressure_linear_solver =  GMRESSolver(ProjectParameters.Pressure_Iterative_Tolerance, ProjectParameters.Pressure_Solver_Max_Iteration,press_precond)
    elif ProjectParameters.Pressure_Linear_Solver == "Skyline LU factorization":
        pressure_linear_solver = SkylineLUFactorizationSolver()
    elif ProjectParameters.Pressure_Linear_Solver == "Super LU":
        pressure_linear_solver = SuperLUSolver()
    elif ProjectParameters.Pressure_Linear_Solver == "Parallel MKL Pardiso":
        from KratosMultiphysics.MKLSolversApplication import MKLPardisoSolver
        pressure_linear_solver = MKLPardisoSolver()
        
dynamic_tau = ProjectParameters.use_dt_in_stabilization
oss_switch = ProjectParameters.use_orthogonal_subscales
#creating the solvers
#fluid solver
fluidSolver = fluid_solver.IncompressibleFluidSolver(fluid_model_part,domain_size)
fluidSolver.max_val_its = ProjectParameters.max_vel_its
fluidSolver.max_press_its = ProjectParameters.max_press_its
fluidSolver.predictor_corrector = ProjectParameters.predictor_corrector            
fluidSolver.vel_toll = ProjectParameters.velocity_relative_tolerance
fluidSolver.press_toll = ProjectParameters.pressure_relative_tolerance
fluidSolver.dynamic_tau = float(dynamic_tau)
fluidSolver.compute_reactions = ProjectParameters.Calculate_reactions     

fluidSolver.Initialize()     
print "fluid solver created"

#creating a mesh solver
reform_dofs_at_each_step = False
mesh_sol = trilinos_mesh_solver.TrilinosMeshSolver(fluid_model_part,domain_size,reform_dofs_at_each_step)
mesh_sol.Initialize()
print "mesh solver created"

#selecting nodes to be moved
moving_nodes = []
SelectNodes(moving_nodes,fluid_model_part);
ApplyDisplacementConditions(fluid_model_part);
print len(moving_nodes);

# Initialize mesh for results
mesh_name = 0.0 #if we want the mesh to change at each time step then ****mesh_name = time****
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh( fluid_model_part.GetMesh() )
gid_io.FinalizeMesh()      
print fluid_model_part

#######################################33
#preparing output of point graphs
import point_graph_printer

output_nodes_list = define_output.DefineOutputPoints()
graph_printer = point_graph_printer.PrintGraphPrinter(output_nodes_list, fluid_model_part, variables_dictionary, domain_size)


# load settings from project parameter file
Re = ProjectParameters.Re
for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,1.0/Re)
time_t = ProjectParameters.Start_time
final_time = ProjectParameters.max_time
full_Dt= ProjectParameters.Dt
initial_Dt = 0.001 * full_Dt #0.05 #0.01
step = 0

# Solution
max_solution_time = 0
while(time_t <= final_time):

    if(step < 3):
        Dt = initial_Dt
    else:
        Dt = full_Dt
        
    time_t = time_t + Dt
    step = step + 1
    fluid_model_part.CloneTimeStep(time_t)

    print "STEP = ", step
    print "TIME = ", time_t

    if(step >= 3):

        # move nodes is corresponding to a solution of a structural problem
        MoveNodes(moving_nodes,time_t)

	# regularize fluid mesh
	start_time = time.time()
        mesh_sol.Solve();   
	stopped_time = time.time() - start_time
	if(stopped_time > max_solution_time):
		max_solution_time = stopped_time		

	# solve fluid problem
        fluidSolver.Solve()
	
      	#results writing options
       	gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time_t,0)
       	gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time_t,0)
       	gid_io.WriteNodalResults(MESH_VELOCITY,fluid_model_part.Nodes,time_t,0)
       	gid_io.WriteNodalResults(DISPLACEMENT,fluid_model_part.Nodes,time_t,0)

print "##############################################################"
print "Maximal time needed for solution of  mesh movement: ",max_solution_time, "seconds"
print "##############################################################"

gid_io.FinalizeResults()
