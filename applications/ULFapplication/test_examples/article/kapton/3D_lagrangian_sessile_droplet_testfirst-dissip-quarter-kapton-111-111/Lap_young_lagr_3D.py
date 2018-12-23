from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# import the configuration data as read from the GiD
import math
import ProjectParameters
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
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ULFApplication import *
#from KratosMultiphysics.StructuralApplication import *
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


#############################################
##importing the solvers needed
SolverType = ProjectParameters.SolverType

## Choosing element type for lagrangian_model_part
element_type = problem_settings.lagrangian_element
import SurfaceTension_Monolithic_Solver_3D as solver_lagr
SolverSettings = ProjectParameters.FluidSolverConfiguration
solver_lagr = import_solver(SolverSettings)
solver_lagr.AddVariables(lagrangian_model_part, SolverSettings)

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


#input_file_name = "sessile-cube_5"

# 7 below is having 3e-5 overall and (all edgs points and lines of 2.5e-5)
#input_file_name = "sessile-cube-quarter-7" 

#below is having only below points and lines of 3.2e-5 and (overall) is 3.7 e-5 , hard and zero hard
#below is having only below points 3.4e-5 and lines of 3.6e-5 and (overall) is 3.9 e-5 , hard and zero hard
#(only lines 3)e-5
input_file_name = "sessile-cube-quarter-10"


gid_io = GidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)

#READ LAGRANGIAN PART
model_part_io_structure = ModelPartIO(input_file_name)
model_part_io_structure.ReadModelPart(lagrangian_model_part)

compute_reactions=0

if(element_type == "ulf"):
    ulf_fluid.AddDofs(lagrangian_model_part, compute_reactions)
elif(element_type == "SurfaceTension"):
    solver_lagr.AddDofs(lagrangian_model_part, SolverSettings)

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

eul_model_part = 0
gamma = 0.072 		#surface tension coefficient [N m-1]
contact_angle = 75.0 #contact angle [deg]


#dissipative force variables: for example, the first one is the JM model in the x direction, and so on... here we are using the value of zero or one to enable which model and in which direction we are going to use our model.

zeta_dissapative_JM_x = 0.0
zeta_dissapative_JM_y = 0.0
zeta_dissapative_JM_z = 0.0

#zeta_dissapative_JM_x = 1.0*(math.cos(contact_angle*0.0174533)+1.0)
#zeta_dissapative_JM_y = 1.0*(math.cos(contact_angle*0.0174533)+1.0)
#zeta_dissapative_JM_z = 1.0*(math.cos(contact_angle*0.0174533)+1.0)

#########################################################################

#zeta_dissapative_BM_x = 2.0*(math.cos(contact_angle*0.0174533)+1.0)
#zeta_dissapative_BM_y = 2.0*(math.cos(contact_angle*0.0174533)+1.0)
#zeta_dissapative_BM_z = 2.0*(math.cos(contact_angle*0.0174533)+1.0)

zeta_dissapative_BM_x = 0.0
zeta_dissapative_BM_y = 0.0
zeta_dissapative_BM_z = 0.0

##################################################

#zeta_dissapative_SM_x = 0.0
#zeta_dissapative_SM_y = 0.0
#zeta_dissapative_SM_z = 0.0

zeta_dissapative_SM_x = (math.cos(contact_angle*0.0174533)+1.0)/2.0
zeta_dissapative_SM_y = (math.cos(contact_angle*0.0174533)+1.0)/2.0
zeta_dissapative_SM_z = (math.cos(contact_angle*0.0174533)+1.0)/2.0

#these variable are for testing purposes: For example, (testfacta) is either zero or one, this is mainly to see where I am adding the dissipative forces in the surface_tension.h file, and which one will be working fine. 
testfacta = 1.0
testfactb = 1.0
testfactc = 1.0
testfactd = 1.0
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


#lag_solver = solver_lagr.CreateSolver(lagrangian_model_part, SolverSettings, eul_model_part, gamma, contact_angle)
# Mesh solver:
reform_dofs_at_each_step = True

pDiagPrecond = DiagonalPreconditioner()

print("mesh solver created")

lag_solver.alpha_shape = problem_settings.alpha_shape;
lag_solver.echo_level = 1;
lag_solver.Initialize()
print("lagrangian solver created")

lagrangian_model_part.SetBufferSize(3)
    
for node in lagrangian_model_part.Nodes:
  node.SetSolutionStepValue(DENSITY,0, 1000.0) 
  node.SetSolutionStepValue(VISCOSITY, 0, 8.90 * 1e-4)
  node.SetSolutionStepValue(BODY_FORCE_X, 0, 0.0)
  node.SetSolutionStepValue(BODY_FORCE_Y, 0, 0.0)
  node.SetSolutionStepValue(BODY_FORCE_Z, 0, -9.81)
  node.SetSolutionStepValue(PRESSURE,0, 0.0)
  node.SetSolutionStepValue(IS_WATER,0, 1.0) 
  node.SetSolutionStepValue(IS_FLUID,0, 1.0) 
  #z_min = min(node.Z for node in lagrangian_model_part.Nodes)
  #if (node.Z <= z_min ):
      #node.SetSolutionStepValue(IS_STRUCTURE,0,1.0)
  #if (node.Z <= 1e-5 ):
      #node.SetSolutionStepValue(IS_STRUCTURE,0,1.0)
  if (node.GetSolutionStepValue(IS_BOUNDARY) != 0.0 and node.GetSolutionStepValue(IS_STRUCTURE) != 1.0):
      node.SetSolutionStepValue(IS_INTERFACE,0, 1.0)
      node.SetSolutionStepValue(IS_FREE_SURFACE,0, 1.0)
      node.SetSolutionStepValue(FLAG_VARIABLE,0, 1.0)
  if (node.GetSolutionStepValue(IS_BOUNDARY) != 1.0 and node.GetSolutionStepValue(TRIPLE_POINT)!= 0.0):
        node.Set(TO_ERASE,True)
  if (node.GetSolutionStepValue(IS_BOUNDARY) != 1.0 and node.GetSolutionStepValue(IS_STRUCTURE)!= 0.0):
        node.Set(TO_ERASE,True)
  #if (node.GetSolutionStepValue(IS_STRUCTURE) != 0.0):
    #node.SetSolutionStepValue(VELOCITY_X,0, 0.0)
    #node.SetSolutionStepValue(VELOCITY_Y,0, 0.0)
    #node.SetSolutionStepValue(VELOCITY_Z,0, 0.0)
    #node.Fix(VELOCITY_Z)
    #node.Fix(VELOCITY_X)
    #node.Fix(VELOCITY_Y)
  #if (node.GetSolutionStepValue(IS_STRUCTURE) != 0.0 and node.GetSolutionStepValue(TRIPLE_POINT) != 0.0):
        #if ((node.GetSolutionStepValue(CONTACT_ANGLE) < 78.5) and (node.GetSolutionStepValue(CONTACT_ANGLE) > 77.5)):
            #ca=node.GetSolutionStepValue(CONTACT_ANGLE,0)
            #ca = ca
            #node.SetSolutionStepValue(CONTACT_ANGLE,0,ca)
            #node.Fix(CONTACT_ANGLE);
            #node.Fix(VELOCITY_X);
            #node.Fix(VELOCITY_Y);
            #node.Fix(VELOCITY_Z);
            #node.Fix(DISPLACEMENT_X);
            #node.Fix(DISPLACEMENT_Y);
            #node.Fix(DISPLACEMENT_Z);
            #node.Fix(CONTACT_ANGLE);

add_nodes=problem_settings.adaptive_refinement
bulk_modulus=problem_settings.bulk_modulus
density=problem_settings.density
FSI=problem_settings.FSI
##########################################################################################################
Multifile = True
# Initialize .post.list file (GiD postprocess list)
f = open(ProjectParameters.problem_name+'.post.lst','w')
f.write('Multiple\n')   

################################################################

# Stepping and time settings
Dt = ProjectParameters.Dt 
Nsteps  = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0

## Obtain the initial mesh size for the droplet base. This size will be used in tetgen_pfem_refine in order to mantain the element size in the base
##CalculateNodalLength().CalculateNodalLength3D(lagrangian_model_part)
##for node in lagrangian_model_part.Nodes:
  ##if (node.GetSolutionStepValue(IS_STRUCTURE) != 0.0):
    ##mesh_size = node.GetSolutionStepValue(NODAL_LENGTH)
    ##node.SetSolutionStepValue(INITIAL_MESH_SIZE,0,mesh_size)
    ##node.Fix(INITIAL_MESH_SIZE)

while(time <= final_time):
    
    time = time + Dt
    step = step + 1
    lagrangian_model_part.CloneTimeStep(time)

    print("LINEAR SOLVER IS ------------------------------------------------------>", ProjectParameters.Monolithic_Linear_Solver)
    print("STEP = ", step)
    print("TIME = ", time)

    if(step >= 3):
      
      for node in lagrangian_model_part.Nodes:
        #if (node.Z <= 1e-5 ):
            #node.SetSolutionStepValue(IS_STRUCTURE,0,1.0)
        if (node.GetSolutionStepValue(IS_BOUNDARY) != 1.0 and node.GetSolutionStepValue(TRIPLE_POINT)!= 0.0):
            node.Set(TO_ERASE,True)
        #x_min = min(node.X for node in lagrangian_model_part.Nodes)
        #x_max = max(node.X for node in lagrangian_model_part.Nodes)
        #if (node.X >= x_min and node.X <= x_max):
            ##for node in lagrangian_model_part.Nodes:
                #y_min = min(node.Y for node in lagrangian_model_part.Nodes)
                #y_max = max(node.Y for node in lagrangian_model_part.Nodes)
                #if (node.Y >= y_min and node.Y <= y_max):
                    ##for node in lagrangian_model_part.Nodes:
                        #z_min = min(node.Z for node in lagrangian_model_part.Nodes)
                        #if (node.Z <= (z_min + 1e-6) ):
                            #z_min = node.Z
                            #if (node.Z <= (z_min + 1e-6) ):
                                #node.SetSolutionStepValue(IS_STRUCTURE,0,1.0)
                                #node.SetSolutionStepValue(IS_STRUCTURE,1,1.0)
                                #print(node.Id, node.X, node.Y, node.Z, z_min)

        if (node.GetSolutionStepValue(IS_BOUNDARY) != 0.0 and node.GetSolutionStepValue(IS_STRUCTURE) != 1.0):
            node.SetSolutionStepValue(IS_INTERFACE,0, 1.0)
            node.SetSolutionStepValue(IS_FREE_SURFACE,0, 1.0)
            node.SetSolutionStepValue(FLAG_VARIABLE,0, 1.0)
        if (node.GetSolutionStepValue(IS_BOUNDARY) != 1.0 and node.GetSolutionStepValue(TRIPLE_POINT)!= 0.0):
            node.Set(TO_ERASE,True)
        if (node.GetSolutionStepValue(IS_BOUNDARY) != 1.0 and node.GetSolutionStepValue(IS_STRUCTURE)!= 0.0):
            node.Set(TO_ERASE,True)
            
        #if (node.GetSolutionStepValue(IS_STRUCTURE) != 0.0 and node.GetSolutionStepValue(TRIPLE_POINT) != 0.0):
          #if ((node.GetSolutionStepValue(CONTACT_ANGLE) < 78.5) and (node.GetSolutionStepValue(CONTACT_ANGLE) > 77.5)):
            #ca=node.GetSolutionStepValue(CONTACT_ANGLE,0)
            #ca = ca
            #node.SetSolutionStepValue(CONTACT_ANGLE,0,ca)
            #node.Fix(CONTACT_ANGLE);
            #node.Fix(VELOCITY_X);
            #node.Fix(VELOCITY_Y);
            #node.Fix(VELOCITY_Z);
            #node.Fix(DISPLACEMENT_X);
            #node.Fix(DISPLACEMENT_Y);
            #node.Fix(DISPLACEMENT_Z);
            #node.Fix(CONTACT_ANGLE);
      
      lag_solver.Solve();
  

##################################################
##################################################

    if(output_time <= out):
      
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
        gid_io.WriteNodalResults(DISPLACEMENT,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(FLAG_VARIABLE,lagrangian_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(FORCE,lagrangian_model_part.Nodes,time,0)
        #gid_io.WriteNodalResults(INITIAL_MESH_SIZE,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_BOUNDARY,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_FREE_SURFACE,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_INTERFACE,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_STRUCTURE,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(MEAN_CURVATURE_3D,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(NODAL_H,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(NODAL_LENGTH,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(NORMAL,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(TRIPLE_POINT,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,lagrangian_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(IS_FLUID,lagrangian_model_part.Nodes,time,0)
        
        gid_io.Flush()
        gid_io.FinalizeResults();

    out = out + Dt

if Multifile:
    f.close()
else:
    gid_io.FinalizeResults()
