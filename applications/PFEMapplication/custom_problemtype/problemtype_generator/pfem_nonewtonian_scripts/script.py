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
##from KratosStructuralApplication import *
from KratosIncompressibleFluidApplication import *
from KratosPFEMApplication import *
from KratosMeshingApplication import *
# from now on the order is not anymore crucial
##################################################################
##################################################################

#defining a model part
model_part = ModelPart("FluidPart");  

#importing the solver files and adding the variables
import monolithic_solver_lagrangian_nonnewtonian
monolithic_solver_lagrangian_nonnewtonian.AddVariables(model_part)

#reading a model
name = pfem_nonewtonian_var.problem_name

gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO(name,gid_mode,multifile,deformed_mesh_flag, write_conditions)
model_part_io = ModelPartIO(name)
model_part_io.ReadModelPart(model_part)

##assigning density viscosity and gravity acceleration
for node in model_part.Nodes:
    node.SetSolutionStepValue(DENSITY,0,pfem_nonewtonian_var.Density)
    node.SetSolutionStepValue(VISCOSITY,0,pfem_nonewtonian_var.Viscosity) 
    node.SetSolutionStepValue(BODY_FORCE_X,0,pfem_nonewtonian_var.Gravity_X) 
    node.SetSolutionStepValue(BODY_FORCE_Y,0,pfem_nonewtonian_var.Gravity_Y) 
    node.SetSolutionStepValue(BODY_FORCE_Z,0,pfem_nonewtonian_var.Gravity_Z)
    node.SetSolutionStepValue(INTERNAL_FRICTION_ANGLE,0,pfem_nonewtonian_var.Friction_Angle_Tangent)
    node.SetSolutionStepValue(YIELD_STRESS,0,pfem_nonewtonian_var.Yield_Stress)
    
#initializing the mesh
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
final_time = pfem_nonewtonian_var.max_time

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)


#adding dofs
monolithic_solver_lagrangian_nonnewtonian.AddDofs(model_part)
#defining the solver
solver = monolithic_solver_lagrangian_nonnewtonian.MonolithicSolver(model_part,domain_size,box_corner1,box_corner2)
#setting solver parameters
solver.oss_switch = int(pfem_nonewtonian_var.use_oss) #OSS = 1; ASGS = 0;
solver.dynamic_tau = pfem_nonewtonian_var.dynamic_tau 
solver.regularization_coef = pfem_nonewtonian_var.m_coef #m regularization coefficient in the exponential law of viscosity
solver.echo_level = 2
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

          
        
