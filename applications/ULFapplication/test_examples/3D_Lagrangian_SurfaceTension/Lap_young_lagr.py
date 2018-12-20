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

##############################################################
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

import math

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

#lagrangian_model_part.AddNodalSolutionStepVariable(ADHESION_FORCE)
#lagrangian_model_part.AddNodalSolutionStepVariable(AUX_VEL)
#lagrangian_model_part.AddNodalSolutionStepVariable(BODY_FORCE)
#lagrangian_model_part.AddNodalSolutionStepVariable(CONTACT_ANGLE)
#lagrangian_model_part.AddNodalSolutionStepVariable(CURVATURE)
#lagrangian_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
#lagrangian_model_part.AddNodalSolutionStepVariable(DISTANCE)
#lagrangian_model_part.AddNodalSolutionStepVariable(DENSITY)
#lagrangian_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
#lagrangian_model_part.AddNodalSolutionStepVariable(GAUSSIAN_CURVATURE)
#lagrangian_model_part.AddNodalSolutionStepVariable(IS_FLUID)
#lagrangian_model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
#lagrangian_model_part.AddNodalSolutionStepVariable(IS_WATER)
#lagrangian_model_part.AddNodalSolutionStepVariable(MEAN_CURVATURE)
#lagrangian_model_part.AddNodalSolutionStepVariable(NORMAL_EQ)
#lagrangian_model_part.AddNodalSolutionStepVariable(NORMAL_GEOMETRIC)
#lagrangian_model_part.AddNodalSolutionStepVariable(NORMAL_TP)
#lagrangian_model_part.AddNodalSolutionStepVariable(PRINCIPAL_CURVATURE_1)
#lagrangian_model_part.AddNodalSolutionStepVariable(PRINCIPAL_CURVATURE_2)
#lagrangian_model_part.AddNodalSolutionStepVariable(SOLID_FRACTION_GRADIENT)
#lagrangian_model_part.AddNodalSolutionStepVariable(VELOCITY)
#lagrangian_model_part.AddNodalSolutionStepVariable(VISCOSITY)
#lagrangian_model_part.AddNodalSolutionStepVariable(VISCOUS_STRESSX)
#lagrangian_model_part.AddNodalSolutionStepVariable(VISCOUS_STRESSY)

#############################################
#importing the solvers needed

## Choosing element type for lagrangian_model_part
element_type = problem_settings.lagrangian_element
#import monolithic_embedded_solver as solver_lagr
import SurfaceTension_Monolithic_Solver_3D as solver_lagr
#solver_lagr.AddVariables(lagrangian_model_part)
SolverSettings = ProjectParameters.FluidSolverConfiguration
solver_lagr = import_solver(SolverSettings)
solver_lagr.AddVariables(lagrangian_model_part, SolverSettings)
#import mesh_solver
#mesh_solver.AddVariables(lagrangian_model_part)
    
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

#READ LAGRANGIAN PART
lag_model_part = "sphere" # "sphere" | "step3D" | "cube3D" | "cube3D_coarse"
gid_io = GidIO(lag_model_part,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io_structure = ModelPartIO(lag_model_part)
model_part_io_structure.ReadModelPart(lagrangian_model_part)

compute_reactions=0

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

eul_model_part = 0
gamma = 1.0		#surface tension coefficient [N m-1]
contact_angle = 95.0 #contact angle [deg]


#dissipative force variables: for example, the first one is the JM model in the x direction, and so on... here we are using the value of zero or one to enable which model and in which direction we are going to use our model.
zeta_dissapative_JM_x = 0.0
zeta_dissapative_JM_y = 0.0
zeta_dissapative_JM_z = 0.0

zeta_dissapative_BM_x = 0.0
zeta_dissapative_BM_y = 0.0
zeta_dissapative_BM_z = 0.0

zeta_dissapative_SM_x = 0.0
zeta_dissapative_SM_y = 0.0
zeta_dissapative_SM_z = 0.0

#these variable are for testing purposes: For example, (testfacta) is either zero or one, this is mainly to see where I am adding the dissipative forces in the surface_tension.h file, and which one will be working fine. 
testfacta = 0.0
testfactb = 0.0
testfactc = 0.0
testfactd = 0.0
testfacte = 0.0

################################################################

### the dissipative force is added utilizing the power law model, utilizing the capillary number when the velocity is the tangential component at the contact line,

##you can choose one of the three models as:

# option 1: assign: zeta_dissapative_JM = 1.0, for Jiang's Model : tanh(4.96 Ca^(0.702))

## or

# option 2: assign: zeta_dissapative_BM = 1.0, for Bracke's model : 2.24 ca ^(0.54)

## or 

# option 3 :assign: zeta_dissapative_SM = 1.0, for Seeberg's model: 2.24 ca ^(0.54) for Ca > 10^(-3), otherwise, 4.47 Ca^(0.42)

##or

# option 4: no dissipative force is added, by assigning zeta_dissapative_JM, zeta_dissapative_JB, and zeta_dissapative_SM to zero
###references:

#Manservisi S, Scardovelli R. A variational approach to the contact angle dynamics of spreading droplets. Computers & Fluids. 2009 Feb 1;38(2):406-24.
#Buscaglia GC, Ausas RF. Variational formulations for surface tension, capillarity and wetting. Computer Methods in Applied Mechanics and Engineering. 2011 Oct 15;200(45-46):3011-25.


################################################################


#lag_solver = solver_lagr.CreateSolver(lagrangian_model_part, SolverSettings, eul_model_part, gamma, contact_angle, zeta_dissapative_JM, zeta_dissapative_BM, zeta_dissapative_SM)

lag_solver = solver_lagr.CreateSolver(lagrangian_model_part, SolverSettings, eul_model_part, gamma, contact_angle, zeta_dissapative_JM_x, zeta_dissapative_BM_x, zeta_dissapative_SM_x, zeta_dissapative_JM_y, zeta_dissapative_BM_y, zeta_dissapative_SM_y, zeta_dissapative_JM_z, zeta_dissapative_BM_z, zeta_dissapative_SM_z, testfacta, testfactb, testfactc, testfactd, testfacte)
# Mesh solver:
reform_dofs_at_each_step = False
#mesh_solver = mesh_solver.MeshSolver(lagrangian_model_part, 2, reform_dofs_at_each_step)
#pDiagPrecond = DiagonalPreconditioner()
#mesh_solver.linear_solver = CGSolver(1e-3, 300, pDiagPrecond)
#mesh_solver.time_order = 2
#mesh_solver.Initialize()
#mesh_solver.solver.SetEchoLevel(0)
#print("mesh solver created")

lag_solver.alpha_shape = problem_settings.alpha_shape;
lag_solver.echo_level = 2;
lag_solver.Initialize()
print("lagrangian solver created")

#lag_solver.echo_level = 2;

#setting up the buffer size: SHOULD BE DONE AFTER READING!!!
lagrangian_model_part.SetBufferSize(3)


for node in lagrangian_model_part.Nodes:
  node.SetSolutionStepValue(DENSITY,0, 1000.0) 
  node.SetSolutionStepValue(VISCOSITY,0, 0.00015)
  node.SetSolutionStepValue(BODY_FORCE_X, 0, 0.0)
  node.SetSolutionStepValue(BODY_FORCE_Y, 0, 0.0)
  node.SetSolutionStepValue(BODY_FORCE_Z, 0, 0.0)
  node.SetSolutionStepValue(PRESSURE,0, 0.0)
  node.SetSolutionStepValue(IS_WATER,0, 1.0) 
  node.SetSolutionStepValue(IS_FLUID,0, 1.0) 
  if (node.GetSolutionStepValue(IS_BOUNDARY) != 0.0):
      node.SetSolutionStepValue(IS_INTERFACE,0, 1.0)
      node.SetSolutionStepValue(IS_FREE_SURFACE,0, 1.0)
      node.SetSolutionStepValue(FLAG_VARIABLE,0, 1.0)


###########################################################################################################
# Stepping and time settings
Dt = ProjectParameters.Dt 
Nsteps  = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0

Multifile=True   
	    
while(time <= final_time):
    
    time = time + Dt
    step = step + 1
    lagrangian_model_part.CloneTimeStep(time)
    
    print("LINEAR SOLVER IS ------------------------------------------------------>", ProjectParameters.Monolithic_Linear_Solver)
    print("STEP = ", step)
    print("TIME = ", time)
    
    #loop to reset some variables
    for node in lagrangian_model_part.Nodes:
        node.SetSolutionStepValue(FORCE_X,0,0.0)
        node.SetSolutionStepValue(FORCE_Y,0,0.0)
        node.SetSolutionStepValue(FORCE_Z,0,0.0)    

    if(step >= 3):
        lag_solver.Solve()
	#mesh_solver.Solve()

    if(output_time <= out):
        out = 0
        
        res_name1 = str("droplet")
        gid_io.ChangeOutputName(res_name1)
        gid_io.InitializeMesh( time );
        gid_io.WriteNodeMesh((lagrangian_model_part).GetMesh());
        gid_io.WriteMesh((lagrangian_model_part).GetMesh());
        gid_io.FinalizeMesh();
        gid_io.InitializeResults(time, (lagrangian_model_part).GetMesh());

        gid_io.WriteNodalResults(DISPLACEMENT,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(FORCE,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_FLUID,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(MEAN_CURVATURE_3D,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,lagrangian_model_part.Nodes,time,0)  
        
        gid_io.Flush()
        gid_io.FinalizeResults();

    out = out + Dt 
