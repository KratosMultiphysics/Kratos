##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2

##################################################################
##################################################################
## ATTENTION: here the order is important

#including kratos path
kratos_libs_path = '../../../../libs' ##kratos_root/libs
kratos_applications_path = '../../../../applications' ##kratos_root/applications
import sys
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)

#importing Kratos main library
from Kratos import *
kernel = Kernel()   #defining kernel

#importing applications
import applications_interface
applications_interface.Import_ConvectionDiffusionApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosConvectionDiffusionApplication import *

#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here when the "very fix container will be ready"

#reading a model
#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("square",gid_mode,multifile,deformed_mesh_flag, write_conditions)
gid_io.ReadModelPart(model_part)

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()
print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#importing the solver files
import convection_diffusion_solver
convection_diffusion_solver.AddVariables(model_part)

    
#creating a fluid solver object
solver = convection_diffusion_solver.ConvectionDiffusionSolver(model_part,domain_size)
solver.time_order = 2
#pDiagPrecond = DiagonalPreconditioner()
#solver.linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)
#solver.linear_solver = SkylineLUFactorizationSolver();

solver.Initialize()

#assigning the fluid properties
conductivity = 0.1;
density = 1000.0;
specific_heat = 1.0;
for node in model_part.Nodes:
    node.SetSolutionStepValue(CONDUCTIVITY,0,conductivity);
    node.SetSolutionStepValue(DENSITY,0,density);
    node.SetSolutionStepValue(SPECIFIC_HEAT,0,specific_heat);

#assigning a zero velocity field
##vel = Vector(3);
##vel[0] = 0.0; vel[1]= 0.0; vel[2]=0.0
##for node in model_part.Nodes:
##    node.SetSolutionStepValue(VELOCITY,0,vel);
##    node.Free(TEMPERATURE);
##
vel = Vector(3);
for node in model_part.Nodes:
##    vel[0] = -node.Y
##    vel[1] = node.X
##    vel[2] = 0.00
##    node.SetSolutionStepValue(VELOCITY,0,vel);
    if(node.X**2 + node.Y**2 < 0.249999):
        vel[0] = -node.Y
        vel[1] = node.X
        vel[2] = 0.00
        node.SetSolutionStepValue(VELOCITY,0,vel);
    else:
        vel[0] = 0.0
        vel[1] = 0.0
        vel[2] = 0.0
        node.SetSolutionStepValue(VELOCITY,0,vel);

for node in model_part.Nodes:
    node.Free(TEMPERATURE);
    
#applying a temperature of 100
for node in model_part.Nodes:
    if(node.X > 0.499):
        node.Fix(TEMPERATURE)
        node.SetSolutionStepValue(TEMPERATURE,0,100.0);
    if(node.X < -0.499):
        node.Fix(TEMPERATURE)
        node.SetSolutionStepValue(TEMPERATURE,0,0.0);

 


#settings to be changed
nsteps = 300
output_step = 10

Dt = 0.1

out = 0


for step in range(0,nsteps):
    print "line49"

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 3):
        solver.Solve()

    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(TEMP_CONV_PROJ,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(MESH_VELOCITY,model_part.Nodes,time,0)
        out = 0
    out = out + 1

          
        

