from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# import the configuration data as read from the GiD
import ProjectParameters
#import define_output
import problem_settings


def PrintResults(model_part):
    print("Writing results. Please run Gid for viewing results of analysis.")
    for variable_name in ProjectParameters.nodal_results:
        gid_io.WriteNodalResults(variables_dictionary[variable_name],model_part.Nodes,time,0)
    for variable_name in ProjectParameters.gauss_points_results:
        gid_io.PrintOnGaussPoints(variables_dictionary[variable_name],model_part,time)

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

##################################################################
##################################################################
import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
#from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ULFApplication import *
#from KratosMultiphysics.StructuralApplication import *
#from KratosMultiphysics.PFEMApplication import *
#from KratosMultiphysics.ALEApplication import *

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
lagrangian_model_part = ModelPart("LagrangianPart");


SolverType=problem_settings.SolverType
if (SolverType=="Incompressible_Modified_FracStep"):
    fluid_only_model_part = ModelPart("FluidOnlyPart");



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
#lagrangian_model_part.AddNodalSolutionStepVariable(INITIAL_AVERAGE_NODAL_AREA)

#############################################
SolverType = ProjectParameters.SolverType

## Choosing element type for lagrangian_model_part
element_type = problem_settings.lagrangian_element
#import monolithic_embedded_solver as solver_lagr
import SurfaceTension_monolithic_injection_solver as solver_lagr
#solver_lagr.AddVariables(lagrangian_model_part)
SolverSettings = ProjectParameters.FluidSolverConfiguration
solver_lagr = import_solver(SolverSettings)
solver_lagr.AddVariables(lagrangian_model_part, SolverSettings)
#import mesh_solver
#mesh_solver.AddVariables(lagrangian_model_part)
    
#introducing input file name

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
    
input_file_name = "90deg_1mm_2"
    
gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)

#READ LAGRANGIAN PART
model_part_io_structure = ModelPartIO(input_file_name)
model_part_io_structure.ReadModelPart(lagrangian_model_part)

compute_reactions=0

if(element_type == "ulf"):
    ulf_fluid.AddDofs(lagrangian_model_part, compute_reactions)
elif(element_type == "SurfaceTension"):
    solver_lagr.AddDofs(lagrangian_model_part, SolverSettings)
    #mesh_solver.AddDofs(lagrangian_model_part)

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


zeta_dissapative_JM = 0.0
zeta_dissapative_BM = 0.0
zeta_dissapative_SM = 0.0
gamma_sl = 0.0
gamma_sv = 0.0


eul_model_part = 0
#gamma = 0.072 		#surface tension coefficient [N m-1]


surface_temp = 291.2

#gamma = 0.072 		#surface tension coefficient [N m-1]
gamma = 0.07275*(1.0-0.00212*(surface_temp-293.15))
contact_angle = 360.0 	#contact angle [deg]


contact_angle = 90.0 	#contact angle [deg]
#lag_solver = solver_lagr.CreateSolverEmbedded(lagrangian_model_part, SolverSettings, eul_model_part,box_corner1, box_corner2, add_nodes, gamma, contact_angle)

#lag_solver = solver_lagr.CreateSolver(lagrangian_model_part, SolverSettings, eul_model_part, box_corner1, box_corner2, add_nodes, gamma, contact_angle, zeta_dissapative_JM, zeta_dissapative_BM, zeta_dissapative_SM)

lag_solver = solver_lagr.CreateSolver(lagrangian_model_part, SolverSettings, eul_model_part, box_corner1, box_corner2, add_nodes, gamma, contact_angle, zeta_dissapative_JM, zeta_dissapative_BM, zeta_dissapative_SM, surface_temp)

#lag_solver = solver_lagr.CreateSolver(lagrangian_model_part, SolverSettings, eul_model_part)
# Mesh solver:
reform_dofs_at_each_step = False
pDiagPrecond = DiagonalPreconditioner()
##mesh_solver = mesh_solver.MeshSolver(lagrangian_model_part, 2, reform_dofs_at_each_step)
#pDiagPrecond = DiagonalPreconditioner()
#mesh_solver.linear_solver = CGSolver(1e-3, 300, pDiagPrecond)
#mesh_solver.time_order = 2
#mesh_solver.Initialize()
#mesh_solver.solver.SetEchoLevel(0)
#print("mesh solver created")

lag_solver.alpha_shape = problem_settings.alpha_shape;
lag_solver.echo_level = 1;
lag_solver.Initialize()
print("lagrangian solver created")
lagrangian_model_part.SetBufferSize(3)

#solver.echo_level = 2;

##############################################################################
### List of parameters:
v_inlet = 0.04


### Water inlet velocity coefficients (for parabolic profile):
v_max = 0.342
D_pore = 600	#in um

C_coef = v_max
A_coef = -C_coef/(((D_pore/2.0)*1.0E-6)**2)


##############################################################################

inlet_vel=Vector(3)
inlet_vel[0]=0.0
inlet_vel[1]=v_inlet
inlet_vel[2]=0.0

dummy=LagrangianInletProcess(lagrangian_model_part, 0.0, inlet_vel)

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
lagrangian_model_part.SetBufferSize(3)

#solver.AddDofs(eulerian_model_part)



for node in lagrangian_model_part.Nodes:
  node.SetSolutionStepValue(DENSITY,0, 1000.0) 
  node.SetSolutionStepValue(VISCOSITY,0, 0.00015)
  node.SetSolutionStepValue(BODY_FORCE_X, 0, 0.0)
  node.SetSolutionStepValue(BODY_FORCE_Y, 0, -9.81)
  node.SetSolutionStepValue(PRESSURE,0, 0.0)
  node.SetSolutionStepValue(IS_WATER,0, 1.0) 
  #if (node.Y < 0.0000001 and (node.X > -50e-6 and node.X < 50e-6)): #100um diameter pore
  #if (node.Y < 0.0000001 and (node.X > -15e-5 and node.X < 15e-5)): #0.6mm drop
  if (node.Y < 0.0000001 and (node.X >= -65e-5 and node.X <= 65e-5)): #1.2mm drop
      node.SetSolutionStepValue(IS_LAGRANGIAN_INLET,0, 1.0)
      node.SetSolutionStepValue(IS_FREE_SURFACE,0, 0.0)
      node.SetSolutionStepValue(IS_INTERFACE,0, 0.0)
      node.SetSolutionStepValue(FLAG_VARIABLE,0, 0.0)
      v_parab = (0.342-3800000*node.X*node.X)*0.01;
      #node.SetSolutionStepValue(VELOCITY_Y,0, v_parab)
      v_parab = 0.1
      node.Fix(VELOCITY_Y)
      node.SetSolutionStepValue(DISPLACEMENT_X,0, 0.0)
      node.SetSolutionStepValue(DISPLACEMENT_Y,0, 0.0)
      node.Fix(DISPLACEMENT_X)
      node.Fix(DISPLACEMENT_Y)
  else:
      node.SetSolutionStepValue(IS_LAGRANGIAN_INLET,0, 0.0)
  if (node.GetSolutionStepValue(IS_LAGRANGIAN_INLET) == 1.0):
      node.SetSolutionStepValue(IS_STRUCTURE,0, 0.0)
      node.SetSolutionStepValue(VELOCITY_X,0, 0.0)
      node.SetSolutionStepValue(VELOCITY_Y,0, 0.0)
      node.Fix(VELOCITY_X)
      node.Fix(VELOCITY_Y)




##########################################################################################################
time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
lin_solver =  SkylineLUFactorizationSolver()

utilities = VariableUtils()
AAA = MeshTransfer2D()



ReformDofSet=True


###############################################################

cut_model_part = ModelPart("CutPart");
Multifile = True
# Initialize .post.list file (GiD postprocess list)
f = open(ProjectParameters.problem_name+'.post.lst','w')
f.write('Multiple\n')

#######################################


#ProjDirichletLinStrat=ResidualBasedLinearStrategy(aux_model_part, time_scheme, lin_solver, False, ReformDofSet, False, False)
#ProjDirichletLinStrat.SetEchoLevel(2)

Multifile = True


#######################################

# Stepping and time settings
Dt_lag = ProjectParameters.Dt

# time adaptivity, the wrong way:
const_dt = 1
Dt = const_dt*Dt_lag
# time adaptivity, the good way:
#Dt = ProjectParameters.Dt 
Nsteps  = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0
    
step_lag = 0	   
time_lag = 0
#solve_air = 1 #CHECK radius in coupled_eulerian_ulf_utilities!!!!!!!!!!!!!!
max_step_lag = 10
min_step = 3


num_nodes = 0
total_area = 0.0


#for node in lagrangian_model_part.Nodes:
        #if (node.GetSolutionStepValue(IS_BOUNDARY)!=1.0 and node.GetSolutionStepValue(INITIAL_AVERAGE_NODAL_AREA) <= 0.0035):
            #node.Set(TO_ERASE,True)
        
        #initial_relative_element_size = node.GetSolutionStepValue(INITIAL_AVERAGE_NODAL_AREA)
        #node.SetSolutionStepValue(INITIAL_AVERAGE_NODAL_AREA,0, initial_relative_element_size)
        #node.Fix(INITIAL_AVERAGE_NODAL_AREA)

while(time <= final_time):

    time = time + Dt
    step = step + 1


    print("LINEAR SOLVER IS ------------------------------------------------------>", ProjectParameters.Monolithic_Linear_Solver)
    print("STEP = ", step)
    print("TIME = ", time)

    if(step >= min_step):

        while(step_lag < const_dt):
            time_lag = time_lag + Dt_lag
            step_lag = step_lag + 1
            lagrangian_model_part.CloneTimeStep(time_lag)	
            lag_solver.Solve(dummy)
	    
        step_lag = 0
        
        #for node in lagrangian_model_part.Nodes:
            #if (node.GetSolutionStepValue(IS_BOUNDARY)!=1.0 and node.GetSolutionStepValue(INITIAL_AVERAGE_NODAL_AREA) <= 0.0035):
                #node.Set(TO_ERASE,True)
        
        #num_nodes = 0
        #total_area = 0.0
        #for node in lagrangian_model_part.Nodes:
            #num_nodes = num_nodes + 1
            #total_area = node.GetSolutionStepValue(NODAL_AREA) + total_area
        #average_area = total_area / num_nodes
        #print("Total area is: ",total_area)
        #print("Relative element size is: ",average_area/total_area)
        #node.SetSolutionStepValue(NODAL_AREA, 0, 0.005)
        #node.Fix(NODAL_AREA)
    

   
    if(output_time <= out):

                out = 0

                gid_io.FinalizeResults()
                #f.write(ProjectParameters.problem_name+'_'+str(time)+'.post.bin\n')

                res_name1 = str("water")
                gid_io.ChangeOutputName(res_name1)
                gid_io.InitializeMesh( time );
                gid_io.WriteNodeMesh((lagrangian_model_part).GetMesh());
                gid_io.WriteMesh((lagrangian_model_part).GetMesh());
                gid_io.FinalizeMesh();
                gid_io.InitializeResults(time, (lagrangian_model_part).GetMesh());

                gid_io.WriteNodalResults(CONTACT_ANGLE,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(MEAN_CURVATURE_2D,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(DISPLACEMENT,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(IS_BOUNDARY,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(IS_FREE_SURFACE,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(IS_INTERFACE,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(IS_STRUCTURE,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(FLAG_VARIABLE,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(FORCE,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(NORMAL,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(PRESSURE,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(TRIPLE_POINT,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(VELOCITY,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(IS_LAGRANGIAN_INLET,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(NODAL_H,lagrangian_model_part.Nodes,time,0)
                gid_io.WriteNodalResults(NODAL_AREA,lagrangian_model_part.Nodes,time,0)
                #gid_io.WriteNodalResults(INITIAL_AVERAGE_NODAL_AREA,lagrangian_model_part.Nodes,time,0)
	 
                
                gid_io.Flush()
                gid_io.FinalizeResults();
                
                out = 0

    out = out + Dt

#if Multifile:
    #f.close()
else:
    gid_io.FinalizeResults()
