#import the configuration data as read from the GiD
import pfem_nonewtonian_var

##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = pfem_nonewtonian_var.domain_size

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path            = pfem_nonewtonian_var.kratos_path + '/libs' ##kratos_root/libs
kratos_applications_path    = pfem_nonewtonian_var.kratos_path + '/applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.Import_PFEMApplication = True
applications_interface.Import_MeshingApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)
from KratosStructuralApplication import *
from KratosIncompressibleFluidApplication import *
## from now on the order is not anymore crucial
##################################################################
##################################################################

#defining a model part
model_part = ModelPart("FluidPart");  

#importing the solver files and adding the variables
SolverType = pfem_nonewtonian_var.SolverType

if(SolverType == "pfem_solver_ale"):
    import pfem_solver_ale
    pfem_solver_ale.AddVariables(model_part)
elif(SolverType == "monolithic_solver_lagrangian"):
    import monolithic_solver_lagrangian_nonnewtonian
    monolithic_solver_lagrangian_nonnewtonian.AddVariables(model_part)
else:
    raise "solver type not supported: options are FractionalStep - Monolithic"


#reading a model
name = pfem_nonewtonian_var.problem_name

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io = ModelPartIO(name)
model_part_io.ReadModelPart(model_part)

##check to ensure that no node has zero density or pressure
for node in model_part.Nodes:
    node.SetSolutionStepValue(DENSITY,0,pfem_nonewtonian_var.Density)
    node.SetSolutionStepValue(VISCOSITY,0,pfem_nonewtonian_var.Viscosity) 
    node.SetSolutionStepValue(BODY_FORCE_X,0,pfem_nonewtonian_var.Gravity_X) 
    node.SetSolutionStepValue(BODY_FORCE_Y,0,pfem_nonewtonian_var.Gravity_Y) 
    node.SetSolutionStepValue(BODY_FORCE_Z,0,pfem_nonewtonian_var.Gravity_Z)
    node.SetSolutionStepValue(YIELD_STRESS,0,pfem_nonewtonian_var.YieldStress)
##    if(node.GetSolutionStepValue(DENSITY) == 0.0):
##        print "node ",node.Id," has zero density!"
##        raise 'node with zero density found'
##    if(node.GetSolutionStepValue(VISCOSITY) == 0.0):
##        print "node ",node.Id," has zero viscosity!"
##        raise 'node with zero density found'


mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()

print model_part
print model_part.Properties

#setting the limits of the bounding box
box_corner1 = Vector(3);
box_corner1[0]=pfem_nonewtonian_var.min_x;
box_corner1[1]=pfem_nonewtonian_var.min_y;
box_corner1[2]=pfem_nonewtonian_var.min_z;
box_corner2 = Vector(3);
box_corner2[0]=pfem_nonewtonian_var.max_x;
box_corner2[1]=pfem_nonewtonian_var.max_y;
box_corner2[2]=pfem_nonewtonian_var.max_z;

#time setting
output_Dt = pfem_nonewtonian_var.output_Dt
max_dt = pfem_nonewtonian_var.max_dt
min_dt = pfem_nonewtonian_var.min_dt
safety_factor = pfem_nonewtonian_var.safety_factor
nsteps = pfem_nonewtonian_var.nsteps

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

if(SolverType == "pfem_solver_ale"):
    #adding dofs
    pfem_solver_ale.AddDofs(model_part)
    #creating a fluid solver object
    name = str("dam2d")
    solver = pfem_solver_ale.PFEMSolver(model_part,name,box_corner1,box_corner2,domain_size)
    solver.laplacian_form = pfem_nonewtonian_var.laplacian_form
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
    monolithic_solver_lagrangian_nonnewtonian.AddDofs(model_part)
    solver = monolithic_solver_lagrangian_nonnewtonian.MonolithicSolver(model_part,domain_size,box_corner1,box_corner2)
    oss_swith = pfem_nonewtonian_var.use_oss
    dynamic_tau = pfem_nonewtonian_var.dynamic_tau
    solver.echo_level = 2
    model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_swith);				
    model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    solver.Initialize(output_Dt)


time = 0.0
for step in range(0,nsteps):
    print "line49"

    new_Dt = solver.EstimateDeltaTime(min_dt,max_dt)

    time = time + new_Dt*safety_factor
    
    model_part.CloneTimeStep(time)

    print time

    #solving the fluid problem
    if(step > 3):
        solver.Solve(time,gid_io)


          
        
