from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# import the configuration data as read from the GiD
import ProjectParameters
import problem_settings


def PrintResults(model_part):
    print("Writing results. Please run Gid for viewing results of analysis.")
    for variable_name in ProjectParameters.nodal_results:
        gid_io.WriteNodalResults(variables_dictionary[variable_name],model_part.Nodes,time,0)
    for variable_name in ProjectParameters.gauss_points_results:
        gid_io.PrintOnGaussPoints(variables_dictionary[variable_name],model_part,time)
        
def OutputResults(model_part,initialization):
  
    for node in model_part.Nodes:
	if ( node.GetSolutionStepValue(TRIPLE_POINT) != 0.0 and node.X > 0.0 ):
	    advancing_angle = node.GetSolutionStepValue(CONTACT_ANGLE,0)
	    node.SetSolutionStepValue (CONTACT_ANGLE, 0, advancing_angle)
	    advancing_vel =  node.GetSolutionStepValue(VELOCITY_X,0)
	    node.SetSolutionStepValue (VELOCITY_X, 0, advancing_vel)
	if ( node.GetSolutionStepValue(TRIPLE_POINT) != 0.0 and node.X < 0.0 ):
	    receding_angle = node.GetSolutionStepValue(CONTACT_ANGLE)
	    node.SetSolutionStepValue (CONTACT_ANGLE, 0, receding_angle)
	    receding_vel =  node.GetSolutionStepValue(VELOCITY_X,0)
	    node.SetSolutionStepValue (VELOCITY_X, 0, receding_vel)

    if (initialization == True):
	file = open('contact_angle_values_trouble_shooting','w+')
	file.write('\n \n')
	file.write('Channel Geometry (Length X width X hight):1mmX1mmX60mm')
	file.write('\n \n')
	file.write('inlet velocity = 8.3 m/s')
	file.write('\n \n')
	file.write('initial contact angle values (in degrees):')
	file.write('\n')
	file.write('average_contact_angle = '+str(contact_angle))
	file.write('\n')
	file.write('advancing_angle (maximum advancing angle in the x-direction) = '+str(advancing_angle))
	file.write('\n')
	file.write('receding_angle ( receding angle in the x-direction) =  '+str(receding_angle))
	file.write('advancing_vel (maximum advancing vel in the x-direction) = '+str(advancing_vel))
	file.write('\n')
	file.write('receding_vel ( receding vel in the x-direction) =  '+str(receding_vel))
	file.write('\n============================================================\n')
	file.write('\n')
	file.close()

	# here below to create a file that can be transfered into XL_format
	file = open('contact_angle_values_XL_format','w+')
	file.write('\n \n')
	file.write('Channel Geometry (Length X width X hight):1mmX1mmX60mm')
	file.write('\n \n')
	file.write('initial contact angle values (in degrees):')
	file.write('\n')
	file.write('average_contact_angle = '+str(contact_angle))
	file.write('\n')
	file.write('advancing_angle (maximum advancing angle in the x-direction) = '+str(advancing_angle))
	file.write('\n')
	file.write('receding_angle (receding angle in the x-direction) = '+str(receding_angle))
	file.write('\n')
	file.write('advancing_vel (maximum advancing angle in the x-direction) = '+str(advancing_vel))
	file.write('\n')
	file.write('receding_vel (receding angle in the x-direction) = '+str(receding_vel))
	file.write('\n')
	file.close()
    else:
	file = open('contact_angle_values_trouble_shooting','a+')
	file.write('TIME '+str(time))
	file.write('\n \n')
	file.write(' \n')
	file.write('max_x_angle_advanced_x = '+str(advancing_angle))
	file.write('\n \n')
	file.write(' \n')
	file.write('min_x_angle_receding_x = '+str(receding_angle))
	file.write('\n \n')
	file.write(' \n')
	file.write('max_x_vel_advanced_x = '+str(advancing_vel))
	file.write('\n \n')
	file.write(' \n')
	file.write('receding_vel = '+str(receding_vel))
	file.write('\n===========================================\n')
	file.close()
	
	file = open('contact_angle_values_XL_format','a+')
	file.write(' \n')
	file.write(str(time)+' '+str(receding_angle)+' '+str(advancing_angle)+' '+str(receding_vel)+' '+str(advancing_vel)+' ')
	file.close()

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

##################################################################
##################################################################
import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.PFEMApplication import *
from KratosMultiphysics.ALEApplication import *

## defining variables to be used

variables_dictionary = {"PRESSURE" : PRESSURE,
                        "VELOCITY" : VELOCITY,
                        "REACTION" : REACTION,
                        "DISTANCE" : DISTANCE,
			 "AUX_VEL" : AUX_VEL,                        
                        "DISPLACEMENT" : DISPLACEMENT,
                        "IS_INTERFACE" : IS_INTERFACE,
                        "IS_STRUCTURE" : IS_STRUCTURE,
                        "VISCOUS_STRESSX": VISCOUS_STRESSX,
                        "VISCOUS_STRESSY": VISCOUS_STRESSY,
                        "IS_WATER": IS_WATER,
                        "DENSITY": DENSITY,
                        "VISCOSITY": VISCOSITY}

#defining a model part for the fluid 
eulerian_model_part = ModelPart("FluidPart");
reduced_model_part = ModelPart("RedPart");
lagrangian_model_part = ModelPart("LagrangianPart");
aux_model_part = ModelPart("AuxPart");
interface_part = ModelPart("InterfacePart");

SolverType=problem_settings.SolverType
if (SolverType=="Incompressible_Modified_FracStep"):
    fluid_only_model_part = ModelPart("FluidOnlyPart");


eulerian_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
eulerian_model_part.AddNodalSolutionStepVariable(DISABLE)
eulerian_model_part.AddNodalSolutionStepVariable(DISTANCE)
eulerian_model_part.AddNodalSolutionStepVariable(ACCELERATION)
eulerian_model_part.AddNodalSolutionStepVariable(AUX_VEL)
eulerian_model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
eulerian_model_part.AddNodalSolutionStepVariable(VISCOUS_STRESSX)
eulerian_model_part.AddNodalSolutionStepVariable(VISCOUS_STRESSY)


lagrangian_model_part.AddNodalSolutionStepVariable(DISTANCE)
lagrangian_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
lagrangian_model_part.AddNodalSolutionStepVariable(VELOCITY)
lagrangian_model_part.AddNodalSolutionStepVariable(VISCOSITY)
lagrangian_model_part.AddNodalSolutionStepVariable(DENSITY)
lagrangian_model_part.AddNodalSolutionStepVariable(BODY_FORCE)
lagrangian_model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
lagrangian_model_part.AddNodalSolutionStepVariable(VISCOUS_STRESSX)
lagrangian_model_part.AddNodalSolutionStepVariable(VISCOUS_STRESSY)
lagrangian_model_part.AddNodalSolutionStepVariable(AUX_VEL)
lagrangian_model_part.AddNodalSolutionStepVariable(CONTACT_ANGLE)
lagrangian_model_part.AddNodalSolutionStepVariable(IS_WATER)

#############################################
##importing the solvers needed
SolverType = ProjectParameters.SolverType
if(SolverType == "FractionalStep"):
    import vms_fractional_step_solver as solver
    solver.AddVariables(eulerian_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    import vms_monolithic_solver as solver
    solver.AddVariables(eulerian_model_part)
 
solver.AddVariables(eulerian_model_part)

## Choosing element type for lagrangian_model_part
element_type = problem_settings.lagrangian_element
#import monolithic_embedded_solver as solver_lagr
import vms_monolithic_solver_disp_ptfe_2 as solver_lagr
#solver_lagr.AddVariables(lagrangian_model_part)
SolverSettings = ProjectParameters.FluidSolverConfiguration
solver_lagr = import_solver(SolverSettings)
solver_lagr.AddVariables(lagrangian_model_part, SolverSettings)
import mesh_solver
mesh_solver.AddVariables(lagrangian_model_part)
    
#introducing input file name
input_file_name = ProjectParameters.problem_name

compute_reactions=False
#reading the fluid part

# initialize GiD  I/O
if ProjectParameters.GiDPostMode == "Binary":
    gid_mode = GiDPostMode.GiD_PostBinary
elif ProjectParameters.GiDPostMode == "Ascii":
    gid_mode = GiDPostMode.GiD_PostAscii
elif ProjectParameters.GiDPostMode == "AsciiZipped":
    gid_mode = GiDPostMode.GiD_PostAsciiZipped
else:
    print("Unknown GiD post mode, assuming Binary")
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
    print("Unknown GiD multiple file mode, assuming Single")
    multifile = MultiFileFlag.SingleFile

#READ EULERIAN PART
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(eulerian_model_part)

#READ LAGRANGIAN PART
#model_part_io_structure = ModelPartIO("drop_1mm2D")
model_part_io_structure = ModelPartIO("ptfe-2d1X1-25-275e-5")
model_part_io_structure.ReadModelPart(lagrangian_model_part)

compute_reactions=0

if(element_type == "ulf"):
    ulf_fluid.AddDofs(lagrangian_model_part, compute_reactions)
elif(element_type == "VMS"):
    solver_lagr.AddDofs(lagrangian_model_part, SolverSettings)
    mesh_solver.AddDofs(lagrangian_model_part)

#setting the limits of the bounding box
box_corner1 = Vector(3); 
box_corner1[0]=problem_settings.bounding_box_corner1_x; box_corner1[1]=problem_settings.bounding_box_corner1_y; box_corner1[2]=problem_settings.bounding_box_corner1_z;
box_corner2 = Vector(3); 
box_corner2[0]=problem_settings.bounding_box_corner2_x; box_corner2[1]=problem_settings.bounding_box_corner2_y; box_corner2[2]=problem_settings.bounding_box_corner2_z;
#here we write the convergence data..,
outstring2 = "convergence_info.txt"
outputfile1 = open(outstring2, 'w')

add_nodes=problem_settings.adaptive_refinement
bulk_modulus=problem_settings.bulk_modulus
density=problem_settings.density
FSI=problem_settings.FSI

zeta_dissapative_JM = 1.0
zeta_dissapative_BM = 0.0
zeta_dissapative_SM = 0.0
gamma_sl = 0
gamma_sv = 0


eul_model_part = 0
gamma = 0.072 		#surface tension coefficient [N m-1]
contact_angle = 109.0 	#contact angle [deg]
lag_solver = solver_lagr.CreateSolver(lagrangian_model_part, SolverSettings, eul_model_part, gamma, contact_angle, zeta_dissapative_JM, zeta_dissapative_BM, zeta_dissapative_SM, gamma_sl, gamma_sv)
# Mesh solver:
reform_dofs_at_each_step = False
mesh_solver = mesh_solver.MeshSolver(lagrangian_model_part, 2, reform_dofs_at_each_step)
pDiagPrecond = DiagonalPreconditioner()
mesh_solver.linear_solver = CGSolver(1e-3, 300, pDiagPrecond)
mesh_solver.time_order = 2
mesh_solver.Initialize()
mesh_solver.solver.SetEchoLevel(0)
print("mesh solver created")

lag_solver.alpha_shape = problem_settings.alpha_shape;
lag_solver.echo_level = 2;
lag_solver.Initialize()
print("lagrangian solver created")

solver.echo_level = 2;

inlet_vel=Vector(3)
inlet_vel[0]=0.0
inlet_vel[1]=0.0
inlet_vel[2]=0.0

dummy=LagrangianInletProcess(lagrangian_model_part, 0.0, inlet_vel)

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
eulerian_model_part.SetBufferSize(3)
reduced_model_part.SetBufferSize(3)
lagrangian_model_part.SetBufferSize(3)

##adding dofs

solver.AddDofs(eulerian_model_part)

########################### SIMULATION PARAMETERS ####################################
density_water = 1000.0
viscosity_water = 0.00015

density_air = 1.2
viscosity_air = 0.0001

channel_length = 0.06
left_edge = -0.03
right_edge = 0.03
upper_edge = 0.001
lower_edge = 0.0
######################################################################################

for node in lagrangian_model_part.Nodes:
  node.SetSolutionStepValue(DENSITY,0, density_water) 
  node.SetSolutionStepValue(VISCOSITY,0, viscosity_water)
  node.SetSolutionStepValue(BODY_FORCE_X, 0, 0.0)
  node.SetSolutionStepValue(BODY_FORCE_Y, 0, -9.81)
  node.SetSolutionStepValue(PRESSURE,0, 0.0)
  node.SetSolutionStepValue(IS_WATER,0, 1.0) 
  if (node.GetSolutionStepValue(IS_BOUNDARY) != 0.0):
      node.SetSolutionStepValue(IS_INTERFACE,0, 1.0)
      node.SetSolutionStepValue(IS_FREE_SURFACE,0, 1.0)
      node.SetSolutionStepValue(FLAG_VARIABLE,0, 1.0)
  if (node.Y <=0.0):
      node.SetSolutionStepValue(IS_STRUCTURE,0, 1.0)
      node.SetSolutionStepValue(VELOCITY_X,0, 0.0)
      node.SetSolutionStepValue(VELOCITY_Y,0, 0.0)
      node.Free(VELOCITY_X)
      node.Fix(VELOCITY_Y)


for node in eulerian_model_part.Nodes:
    node.SetSolutionStepValue(DENSITY, 0, density_air)
    node.SetSolutionStepValue(VISCOSITY, 0, viscosity_air)
    if (node.X <=left_edge):# and node.Y > 0.00000001 and node.Y<0.029999): #### elaf changed node.X < -0.059 to (node.X <= -0.06)
      node.SetSolutionStepValue(VELOCITY_X,0, 1.0);
      node.SetSolutionStepValue(VELOCITY_Y,0, 0.0)
      node.Fix(VELOCITY_X)
      node.Fix(VELOCITY_Y)
        
#PRESCRIBE OUTLET PRESSURE!!!
for node in eulerian_model_part.Nodes:
    if (node.X>=right_edge):
      node.SetSolutionStepValue(PRESSURE, 0, 0.0)
      node.Fix(PRESSURE)
        
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
        
elif "monolithic_solver_eulerian": # single coupled solver
    # preconditioner
    try:
        if ProjectParameters.Monolithic_Preconditioner_type == 'Diagonal':
            precond = DiagonalPreconditioner()
        elif ProjectParameters.Monolithic_Preconditioner_type == 'ILU0':
            precond = ILU0Preconditioner()
    except AttributeError:
        # If we are using a direct solver Pressure_Preconditioner_type will not
        # exist. This is not an error.
        precond = None

    # solver
    if ProjectParameters.Monolithic_Linear_Solver == "Conjugate gradient":
        if ProjectParameters.Monolithic_Preconditioner_type == 'None':
            monolithic_linear_solver =  CGSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration)
        else:
            monolithic_linear_solver =  CGSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration,precond)
    elif ProjectParameters.Monolithic_Linear_Solver == "BiConjugate gradient stabilized":
        if ProjectParameters.Monolithic_Preconditioner_type == 'None':
            monolithic_linear_solver =  BICGSTABSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration)
        else:
            monolithic_linear_solver =  BICGSTABSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration,precond)
    elif ProjectParameters.Monolithic_Linear_Solver == "GMRES":
        if ProjectParameters.Monolithic_Preconditioner_type == 'None':
            monolithic_linear_solver =  GMRESSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration)
        else:
            monolithic_linear_solver =  GMRESSolver(ProjectParameters.Monolithic_Iterative_Tolerance, ProjectParameters.Monolithic_Solver_Max_Iteration,precond)
    elif ProjectParameters.Monolithic_Linear_Solver == "Skyline LU factorization":
        monolithic_linear_solver = SkylineLUFactorizationSolver()
    elif ProjectParameters.Monolithic_Linear_Solver == "Super LU":
        monolithic_linear_solver = SuperLUSolver()
    elif ProjectParameters.Monolithic_Linear_Solver == "MixedUP":
        velocity_linear_solver = BICGSTABSolver(1e-5, 500) #SuperLUIterativeSolver()
        pressure_linear_solver =  BICGSTABSolver(1e-3, 500) 
        m = 5
        max_it = m
        monolithic_linear_solver = MixedUPLinearSolver(velocity_linear_solver,pressure_linear_solver,1e-5,max_it,m)
    elif ProjectParameters.Monolithic_Linear_Solver == "SuperLUIterativeSolver":    
        monolithic_linear_solver = ScalingSolver( SuperLUIterativeSolver() , True )
    elif ProjectParameters.Monolithic_Linear_Solver == "Parallel MKL Pardiso":
        from KratosMultiphysics.MKLSolversApplication import MKLPardisoSolver
        monolithic_linear_solver = MKLPardisoSolver()
        
#copy Y_WALL
for node in eulerian_model_part.Nodes:
    y = node.GetSolutionStepValue(Y_WALL,0)
    node.SetValue(Y_WALL,y)

dynamic_tau = ProjectParameters.use_dt_in_stabilization
oss_switch = ProjectParameters.use_orthogonal_subscales
#creating the solvers
#fluid solver
if(SolverType == "FractionalStep"):
    fluid_solver = solver.IncompressibleFluidSolver(eulerian_model_part,domain_size)
    fluid_solver.max_val_its = ProjectParameters.Velocity_Solver_Max_Iteration
    fluid_solver.max_press_its = ProjectParameters.Pressure_Solver_Max_Iteration
    fluid_solver.predictor_corrector = ProjectParameters.predictor_corrector            
    fluid_solver.vel_toll = ProjectParameters.velocity_relative_tolerance
    fluid_solver.press_toll = ProjectParameters.pressure_relative_tolerance
    fluid_solver.dynamic_tau = float(dynamic_tau)
    fluid_solver.compute_reactions = ProjectParameters.Calculate_reactions
    
    fluid_solver_red = solver.IncompressibleFluidSolver(reduced_model_part,domain_size)
    fluid_solver_red.max_vel_its = ProjectParameters.Velocity_Solver_Max_Iteration
    fluid_solver_red.max_press_its = ProjectParameters.Pressure_Solver_Max_Iteration
    fluid_solver_red.predictor_corrector = ProjectParameters.predictor_corrector            
    fluid_solver_red.vel_toll = ProjectParameters.velocity_relative_tolerance
    fluid_solver_red.press_toll = ProjectParameters.pressure_relative_tolerance
    fluid_solver_red.dynamic_tau = float(dynamic_tau)
    fluid_solver_red.compute_reactions = ProjectParameters.Calculate_reactions
    
elif (SolverType == "monolithic_solver_eulerian"):
    eul_model_part = 1
    fluid_solver = solver.MonolithicSolver(eulerian_model_part,domain_size,eul_model_part, gamma, contact_angle)  
    fluid_solver.oss_switch = int(oss_switch)
    fluid_solver.dynamic_tau = float(dynamic_tau)
    fluid_solver.rel_vel_tol = ProjectParameters.velocity_relative_tolerance
    fluid_solver.abs_vel_tol = ProjectParameters.velocity_absolute_tolerance
    fluid_solver.rel_pres_tol = ProjectParameters.pressure_relative_tolerance
    fluid_solver.abs_pres_tol = ProjectParameters.pressure_absolute_tolerance
    fluid_solver.max_iter = ProjectParameters.max_iterations
    fluid_solver.compute_reactions = ProjectParameters.Calculate_reactions
    fluid_solver.linear_solver = monolithic_linear_solver

    fluid_solver_red = solver.MonolithicSolver(reduced_model_part,domain_size,eul_model_part, gamma, contact_angle)
    fluid_solver_red.oss_switch = int(oss_switch)
    fluid_solver_red.dynamic_tau = float(dynamic_tau)
    fluid_solver_red.rel_vel_tol = ProjectParameters.velocity_relative_tolerance
    fluid_solver_red.abs_vel_tol = ProjectParameters.velocity_absolute_tolerance
    fluid_solver_red.rel_pres_tol = ProjectParameters.pressure_relative_tolerance
    fluid_solver_red.abs_pres_tol = ProjectParameters.pressure_absolute_tolerance
    fluid_solver_red.max_iter = ProjectParameters.max_iterations
    fluid_solver_red.compute_reactions = ProjectParameters.Calculate_reactions
    fluid_solver_red.linear_solver = monolithic_linear_solver

fluid_solver.Initialize()     
fluid_solver_red.Initialize()   

print("fluid solver created")
##########################################################################################################
time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
lin_solver =  SkylineLUFactorizationSolver()

utilities = VariableUtils()
AAA = MeshTransfer2D()
CoupledEulerianUlfUtils=CoupledEulerianUlfUtils()


ReformDofSet=True
ProjDirichletLinStrat=ResidualBasedLinearStrategy(aux_model_part, time_scheme, lin_solver, False, ReformDofSet, False, False)
ProjDirichletLinStrat.SetEchoLevel(2)

proj_dirichlet_process=ApplyProjDirichletProcess()
#for saving the pseudo lagrangian part (representation of lagrangian part within Eulerian mesh)
pseudo_lag_process=PseudoLagPartProcess
##########################################################################################################

cut_model_part = ModelPart("CutPart");
Multifile = True
# Initialize .post.list file (GiD postprocess list)
f = open(ProjectParameters.problem_name+'.post.lst','w')
f.write('Multiple\n')

#######################################33
# Stepping and time settings
Dt = ProjectParameters.Dt 
Nsteps  = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0

neigh_finder = FindNodalNeighboursProcess(aux_model_part,9,18)
elem_neigh_f=FindElementalNeighboursProcess(aux_model_part, 2, 10)
cond_neigh_f=FindConditionsNeighboursProcess(aux_model_part, 2, 10)

Multifile=True

for node in eulerian_model_part.Nodes:
  if(node.Y<= lower_edge): #non-slip boundary condition for lower wall
    node.SetSolutionStepValue(VELOCITY_X,0,0.0)
    node.SetSolutionStepValue(VELOCITY_Y,0,0.0)
    node.Fix(VELOCITY_X)
    node.Fix(VELOCITY_Y)
  if(node.Y>= upper_edge): #non-slip boundary condition for upper wall
    node.SetSolutionStepValue(VELOCITY_X,0,0.0)
    node.SetSolutionStepValue(VELOCITY_Y,0,0.0)
    node.Fix(VELOCITY_X)
    node.Fix(VELOCITY_Y)    

OutputResults(lagrangian_model_part,True)
    
while(time <= final_time):
  
    time = time + Dt
    step = step + 1	
    eulerian_model_part.CloneTimeStep(time)
    lagrangian_model_part.CloneTimeStep(time)
    reduced_model_part.ProcessInfo=eulerian_model_part.ProcessInfo
    aux_model_part.ProcessInfo=eulerian_model_part.ProcessInfo
    interface_part.ProcessInfo=eulerian_model_part.ProcessInfo	

    print("LINEAR SOLVER IS ------------------------------------------------------>", ProjectParameters.Monolithic_Linear_Solver)
    print("STEP = ", step)
    print("TIME = ", time)

    if(step >= 3 and step <200):
	lag_solver.Solve()
	mesh_solver.Solve()
	
	OutputResults(lagrangian_model_part,False)
      
    if(step >= 200):
	const = 0.1
	vel_inlet = (1.0/const)*time + 5.0
	#u_const = - 4.0 * vel_inlet * 2 / (H * H)
	#velocity_inlet = u_const * node.Y * (node.Y - H)
	if (vel_inlet > 25.0):
	  vel_inlet = 25.0
	#u_const = - 4.0 * vel_inlet * 2 / (H * H)
	#velocity_inlet = u_const * node.Y * (node.Y - H)
	for node in eulerian_model_part.Nodes:
	    if (node.X <=left_edge and node.Y >=lower_edge):    #### #### elaf changed
		node.SetSolutionStepValue(VELOCITY_X,0, vel_inlet);    
		node.Fix(VELOCITY_X)
	
	lag_solver.Solve()
	mesh_solver.Solve()
	
	OutputResults(lagrangian_model_part,False)
	
       ####################################################################################################
        for node in eulerian_model_part.Nodes:
	    node.SetSolutionStepValue(AUX_VEL_X,0,0.0)
	    node.SetSolutionStepValue(AUX_VEL_Y,0,0.0)
	    node.Free(AUX_VEL_X)
	    node.Free(AUX_VEL_Y)
	    if (node.X > left_edge and node.X< right_edge and node.Y> lower_edge and node.Y< upper_edge):
		node.Free(VELOCITY_X)
		node.Free(VELOCITY_Y)
		node.Free(AUX_VEL_X)
		node.Free(AUX_VEL_Y)
		node.Free(FRACT_VEL_X)
		node.Free(FRACT_VEL_Y)
	    if (node.Y >= upper_edge):
		node.Free(VELOCITY_X)
        for node in reduced_model_part.Nodes:
	    node.SetSolutionStepValue(AUX_VEL_X,0,0.0)
	    node.SetSolutionStepValue(AUX_VEL_Y,0,0.0)
	    node.Free(AUX_VEL_X)
	    node.Free(AUX_VEL_Y)
	    if (node.X >left_edge and node.X< right_edge and node.Y>lower_edge and node.Y<upper_edge):
		node.Free(VELOCITY_X)
		node.Free(VELOCITY_Y)
		node.Free(AUX_VEL_X)
		node.Free(AUX_VEL_Y)
		node.Free(FRACT_VEL_X)
		node.Free(FRACT_VEL_Y)
	    if (node.Y >= upper_edge):
		node.Free(VELOCITY_X)	      
	######################################################################################	
	# We set the velocity of the lagrangian boundary to 0 and project it to eulerian part
	# (we neglect the velocity effect of lagrangian domain onto eulerian domain
	vel=Vector(3)
	zerovel=Vector(3)
	zerovel[0]=0.0;
	zerovel[1]=0.0;
	zerovel[2]=0.0;
	for node in lagrangian_model_part.Nodes:
	    if (node.GetSolutionStepValue(IS_BOUNDARY)!=0):
		vel=node.GetSolutionStepValue(VELOCITY);
		node.SetSolutionStepValue(FRACT_VEL,0, vel)
		node.SetSolutionStepValue(VELOCITY,0, zerovel)	
        CoupledEulerianUlfUtils.FindIntersectionOfEdges(lagrangian_model_part, eulerian_model_part, aux_model_part, interface_part)
        	      
        CoupledEulerianUlfUtils.DisableSubdomain(eulerian_model_part, reduced_model_part)
	neigh_finder.Execute()
        elem_neigh_f.Execute()
        cond_neigh_f.Execute()

        if (aux_model_part.Conditions.Size()>0):
            print("SOLVING THE INTERFACE B.C. PROBLEM")            
            ProjDirichletLinStrat.Solve()       
            
	proj_dirichlet_process.ApplyProjDirichlet(reduced_model_part)
	for node in lagrangian_model_part.Nodes:
	    node.SetSolutionStepValue(EXTERNAL_PRESSURE, 0, 0)	  
	AAA.DirectScalarVarInterpolation(reduced_model_part, lagrangian_model_part, PRESSURE, EXTERNAL_PRESSURE);
	
	#mapping of the viscous stress:
	for node in lagrangian_model_part.Nodes:
	    node.SetSolutionStepValue(VISCOUS_STRESSX_X,0,0.0)
	    node.SetSolutionStepValue(VISCOUS_STRESSX_Y,0,0.0)
	  
	AAA.DirectVectorialVarInterpolation(reduced_model_part, lagrangian_model_part, VISCOUS_STRESSX, VISCOUS_STRESSX);

	for node in reduced_model_part.Nodes:
	    node.SetSolutionStepValue(VISCOUS_STRESSX_X,0,0.0)
	    node.SetSolutionStepValue(VISCOUS_STRESSX_Y,0,0.0)
 
	if (step > 200):
	    fluid_solver_red.Solve()			
	
	# Now we recover the lagrangian boundary velocity
	for node in lagrangian_model_part.Nodes:
	    if (node.GetSolutionStepValue(IS_BOUNDARY)!=0):
		vel=node.GetSolutionStepValue(FRACT_VEL);
		node.SetSolutionStepValue(VELOCITY,0, vel)

##################################################
##################################################

    if(output_time <= out):
      if (step>=3):
        res_name2 = str("AIR")
        gid_io.ChangeOutputName(res_name2)
        gid_io.InitializeMesh( time )
        gid_io.WriteMesh( reduced_model_part.GetMesh() )
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(time,(reduced_model_part).GetMesh())
        PrintResults(reduced_model_part)
        out = 0
        
        gid_io.FinalizeResults()
        f.write(ProjectParameters.problem_name+'_'+str(time)+'.post.bin\n')
  
        res_name1 = str("water")
        gid_io.ChangeOutputName(res_name1)
        gid_io.InitializeMesh( time );
        gid_io.WriteNodeMesh((lagrangian_model_part).GetMesh());
        gid_io.WriteMesh((lagrangian_model_part).GetMesh());
        gid_io.FinalizeMesh();
        gid_io.InitializeResults(time, (lagrangian_model_part).GetMesh());
    
        gid_io.WriteNodalResults(CONTACT_ANGLE,lagrangian_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(CURVATURE,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISPLACEMENT,lagrangian_model_part.Nodes,time,0)        
        #gid_io.WriteNodalResults(FORCE,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(EXTERNAL_PRESSURE,lagrangian_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(IS_FREE_SURFACE,lagrangian_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(IS_INTERFACE,lagrangian_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(IS_STRUCTURE,lagrangian_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(NODAL_H,lagrangian_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(NODAL_LENGTH,lagrangian_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(NORMAL,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,lagrangian_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(TRIPLE_POINT,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,lagrangian_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(VISCOUS_STRESSX,lagrangian_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(VISCOUS_STRESSY,lagrangian_model_part.Nodes,time,0)
        
        gid_io.Flush()
        gid_io.FinalizeResults();

    out = out + Dt

if Multifile:
    f.close()
else:
    gid_io.FinalizeResults()
    