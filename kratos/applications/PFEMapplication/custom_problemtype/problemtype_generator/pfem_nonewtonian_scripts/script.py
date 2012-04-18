#import the configuration data as read from the GiD
import problem_settings

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = problem_settings.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
import sys
sys.path.append(problem_settings.kratos_path)

#importing Kratos main library
from KratosMultiphysics import *

from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.PFEMApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.MKLSolversApplication import*

# from now on the order is not anymore crucial
##################################################################
##################################################################

#defining a model part
model_part = ModelPart("FluidPart");  

#importing the solver files and adding the variables
import monolithic_solver_lagrangian_nonnewtonian
monolithic_solver_lagrangian_nonnewtonian.AddVariables(model_part)

#reading a model
name = problem_settings.problem_name

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io = ModelPartIO(name)
model_part_io.ReadModelPart(model_part)

##assigning density viscosity and gravity acceleration
for node in model_part.Nodes:
    node.SetSolutionStepValue(DENSITY,0,problem_settings.Density)
    node.SetSolutionStepValue(VISCOSITY,0,problem_settings.Viscosity) 
    node.SetSolutionStepValue(BODY_FORCE_X,0,problem_settings.Gravity_X) 
    node.SetSolutionStepValue(BODY_FORCE_Y,0,problem_settings.Gravity_Y) 
    node.SetSolutionStepValue(BODY_FORCE_Z,0,problem_settings.Gravity_Z)
    node.SetSolutionStepValue(INTERNAL_FRICTION_ANGLE,0,problem_settings.Friction_Angle_Tangent)
    node.SetSolutionStepValue(YIELD_STRESS,0,problem_settings.Yield_Stress)
    
#initializing the mesh
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()

print model_part
print model_part.Properties

#setting the limits of the bounding box
box_corner1 = Vector(3);
box_corner1[0]=problem_settings.min_x;
box_corner1[1]=problem_settings.min_y;
box_corner1[2]=problem_settings.min_z;
box_corner2 = Vector(3);
box_corner2[0]=problem_settings.max_x;
box_corner2[1]=problem_settings.max_y;
box_corner2[2]=problem_settings.max_z;

#time setting
output_Dt = problem_settings.output_Dt
max_dt = problem_settings.max_dt
min_dt = problem_settings.min_dt
safety_factor = problem_settings.safety_factor
final_time = problem_settings.max_time

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)


#adding dofs
monolithic_solver_lagrangian_nonnewtonian.AddDofs(model_part)
#defining the solver
solver = monolithic_solver_lagrangian_nonnewtonian.MonolithicSolver(model_part,domain_size,box_corner1,box_corner2)
#setting solver parameters
solver.oss_switch = int(problem_settings.use_oss) #OSS = 1; ASGS = 0;
solver.dynamic_tau = problem_settings.dynamic_tau 
solver.regularization_coef = problem_settings.m_coef #m regularization coefficient in the exponential law of viscosity
solver.echo_level = 2

solver.max_iter = 8
print "max iterations non-newtonian solver = 8"
#Initializing the solver
solver.Initialize(output_Dt)

step = 0
time = 0.0
##for step in range(0,nsteps):
while(time < final_time):
##    print "line49"

    
    if(step < 10):
        time = time + min_dt
    else:
        new_Dt = solver.EstimateDeltaTime(min_dt,max_dt)
        time = time + new_Dt*safety_factor
##        print 'Delta Time  = ', new_Dt
##        print 'Safety factor  = ', safety_factor
           
##    print 'Time  = ', time
    
    model_part.CloneTimeStep(time)

    print time

    #solving the fluid problem
    if(step > 3):
        solver.Solve(time,gid_io)

    step = step + 1

          
        
