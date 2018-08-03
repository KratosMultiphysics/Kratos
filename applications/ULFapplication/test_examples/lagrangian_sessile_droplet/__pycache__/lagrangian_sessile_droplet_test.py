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


#############################################
##importing the solvers needed
SolverType = ProjectParameters.SolverType

## Choosing element type for lagrangian_model_part
element_type = problem_settings.lagrangian_element
import SurfaceTension_monolithic_solver as solver_lagr
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

eul_model_part = 0
gamma = 0.072 		#surface tension coefficient [N m-1]
contact_angle = 45.0	#contact angle [deg]

zeta_dissapative_JM = 1.0
zeta_dissapative_BM = 0.0
zeta_dissapative_SM = 0.0



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


lag_solver = solver_lagr.CreateSolver(lagrangian_model_part, SolverSettings, eul_model_part, gamma, contact_angle, zeta_dissapative_JM, zeta_dissapative_BM, zeta_dissapative_SM)



reform_dofs_at_each_step = False
pDiagPrecond = DiagonalPreconditioner()

lag_solver.alpha_shape = problem_settings.alpha_shape;
lag_solver.echo_level = 1;
lag_solver.Initialize()
print("lagrangian solver created")

lagrangian_model_part.SetBufferSize(3)

for node in lagrangian_model_part.Nodes:
  node.SetSolutionStepValue(DENSITY,0, 1000.0)
  node.SetSolutionStepValue(VISCOSITY,0, 8.90 * 1e-4)
  node.SetSolutionStepValue(BODY_FORCE_X, 0, 0.0)
  node.SetSolutionStepValue(BODY_FORCE_Y, 0, -9.81)
  node.SetSolutionStepValue(PRESSURE,0, 0.0)
  node.SetSolutionStepValue(IS_FLUID,0, 1.0)
  if (node.GetSolutionStepValue(IS_BOUNDARY) != 0.0 and node.GetSolutionStepValue(IS_STRUCTURE) != 1.0):
    node.SetSolutionStepValue(IS_FREE_SURFACE,0, 1.0)
    node.SetSolutionStepValue(IS_INTERFACE,0, 1.0)
    node.SetSolutionStepValue(FLAG_VARIABLE,0, 1.0)
  if (node.GetSolutionStepValue(IS_STRUCTURE) != 0.0):
    node.SetSolutionStepValue(VELOCITY_Y,0, 0.0)
    node.Fix(VELOCITY_Y)
    node.Free(VELOCITY_X)



##########################################################################################################
Multifile = True

################################################################

# Stepping and time settings
Dt = ProjectParameters.Dt 
Nsteps  = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0

while(time <= final_time):
    
    time = time + Dt
    step = step + 1
    lagrangian_model_part.CloneTimeStep(time)

    print("LINEAR SOLVER IS ------------------------------------------------------>", ProjectParameters.Monolithic_Linear_Solver)
    print("STEP = ", step)
    print("TIME = ", time)

    if(step >= 3):
      
      for node in lagrangian_model_part.Nodes:
        if (node.GetSolutionStepValue(IS_BOUNDARY) != 0.0 and node.GetSolutionStepValue(IS_STRUCTURE) != 1.0):
            node.SetSolutionStepValue(IS_FREE_SURFACE,0, 1.0)
            node.SetSolutionStepValue(IS_INTERFACE,0, 1.0)
            node.SetSolutionStepValue(FLAG_VARIABLE,0, 1.0)
        if (node.GetSolutionStepValue(IS_STRUCTURE) != 0.0):
            node.SetSolutionStepValue(VELOCITY_Y,0, 0.0)
            node.Fix(VELOCITY_Y)
            if (node.GetSolutionStepValue(TRIPLE_POINT) != 1.0):
                node.SetSolutionStepValue(VELOCITY_X,0, 0.0)
                node.Fix(VELOCITY_X)
            else:
                node.Free(VELOCITY_X)

      lag_solver.Solve()	

##################################################
##################################################

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
        
        gid_io.Flush()
        gid_io.FinalizeResults();

    out = out + Dt

#if Multifile:
    #f.close()
else:
    gid_io.FinalizeResults()
