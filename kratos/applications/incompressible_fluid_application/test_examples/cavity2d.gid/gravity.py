##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs/' ##kratos_root/libs
kratos_applications_path = '../../../../applications/' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_IncompressibleFluidApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosIncompressibleFluidApplication import *


#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here when the "very fix container will be ready"
import incompressible_fluid_solver
incompressible_fluid_solver.AddVariables(model_part)

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("cavity2d",gid_mode,use_multifile,deformed_print_flag,write_conditions)
write_conditions = WriteConditionsFlag.WriteElementsOnly
#gid_io.ReadMesh(model_part.GetMesh())
gid_io.ReadModelPart(model_part)
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name)
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh();
print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

##add Degrees of Freedom to all of the nodes
incompressible_fluid_solver.AddDofs(model_part)

gravity = Vector(3);
gravity[0] = 0.0; gravity[1] = -9.81; gravity[2] = 0.0
zero = Vector(3);
zero[0] = 0.0; zero[1] = 0.0; zero[2] = 0.0
for node in model_part.Nodes:
    if(node.Y > 0.99 and node.X > 0.001 and node.X < 0.999):
        node.Free(VELOCITY_X);
        node.Free(VELOCITY_Y);
        node.Free(VELOCITY_Z);
    node.Free(PRESSURE)
    node.SetSolutionStepValue(BODY_FORCE,0,gravity);
    node.SetSolutionStepValue(VELOCITY,0,zero);
        
#creating a fluid solver object
fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(model_part,domain_size)
fluid_solver.laplacian_form = 3;
fluid_solver.vel_toll = 1e-6
fluid_solver.time_order = 1
fluid_solver.max_press_its = 10;
fluid_solver.predictor_corrector = False

##pILUPrecond = ILU0Preconditioner() 
##fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)

fluid_solver.Initialize()

#settings to be changed
Re = 100.0
nsteps = 100
output_step = 1

for node in model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,0.000001)

if(fluid_solver.laplacian_form == 1):
    for node in model_part.Nodes:
        node.Free(PRESSURE)
        if(node.Y > 0.99):
            node.Fix(PRESSURE)
    

Dt = 0.001
out = 0

gid_io.InitializeResults( 0.0, model_part.GetMesh() )

for step in range(0,nsteps):
    print "line49"

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 3):
        fluid_solver.Solve()


    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESS_PROJ,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(CONV_PROJ,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        out = 0
    out = out + 1

          
gid_io.FinalizeResults();
        

