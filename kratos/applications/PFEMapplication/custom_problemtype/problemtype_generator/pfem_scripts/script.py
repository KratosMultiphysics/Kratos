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
from KratosMultiphysics.ExternalSolversApplication import *

## from now on the order is not anymore crucial
##################################################################
##################################################################

#defining a model part
model_part = ModelPart("FluidPart");  

#importing the solver files and adding the variables
SolverType = problem_settings.SolverType

if(SolverType == "pfem_solver_ale"):
    import pfem_solver_ale
    pfem_solver_ale.AddVariables(model_part)
elif(SolverType == "monolithic_solver_lagrangian"):
    import monolithic_solver_lagrangian
    monolithic_solver_lagrangian.AddVariables(model_part)
else:
    raise "solver type not supported: options are FractionalStep - Monolithic"


#reading a model
name = problem_settings.problem_name

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io = ModelPartIO(name)
model_part_io.ReadModelPart(model_part)

##check to ensure that no node has zero density or pressure
for node in model_part.Nodes:
    node.SetSolutionStepValue(BODY_FORCE_X,0,problem_settings.Gravity_X) 
    node.SetSolutionStepValue(BODY_FORCE_Y,0,problem_settings.Gravity_Y) 
    node.SetSolutionStepValue(BODY_FORCE_Z,0,problem_settings.Gravity_Z)
    node.SetSolutionStepValue(DENSITY,0,problem_settings.density)
    node.SetSolutionStepValue(VISCOSITY,0,problem_settings.viscosity)



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
##nsteps = int(problem_settings.nsteps)
final_time = problem_settings.max_time
#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

if(SolverType == "pfem_solver_ale"):
    #adding dofs
    pfem_solver_ale.AddDofs(model_part)
    #creating a fluid solver object
    solver = pfem_solver_ale.PFEMSolver(model_part,name,box_corner1,box_corner2,domain_size)
    solver.laplacian_form = int(problem_settings.laplacian_form)
    solver.echo_level = 0
    solver.prediction_order = 1
    solver.predictor_corrector = True
    solver.smooth = True
    solver.alpha_shape = 1.2
    solver.max_press_its = 3;  
    #initializing the solver
    initial_dt = 0.001*min_dt
    solver.Initialize(initial_dt,output_Dt)
elif(SolverType == "monolithic_solver_lagrangian"):
    #adding dofs
    monolithic_solver_lagrangian.AddDofs(model_part)
    solver = monolithic_solver_lagrangian.MonolithicSolver(model_part,domain_size,box_corner1,box_corner2)
    solver.oss_swith = int(problem_settings.use_oss)
    solver.dynamic_tau = int(problem_settings.dynamic_tau)
    solver.echo_level = 2

    solver.Initialize(output_Dt)




###############################################################
time = 0.0
step = 0
##for step in range(0,nsteps):
while(time < final_time):
    print "line49"

    new_Dt = solver.EstimateDeltaTime(min_dt,max_dt)

    time = time + new_Dt*safety_factor
    
    model_part.CloneTimeStep(time)

    print time

    #solving the fluid problem
    if(step > 3):
        solver.Solve(time,gid_io)

    step = step + 1 
          
        
