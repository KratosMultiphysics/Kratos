##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 3

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

#reading a model
gid_io = GidIO("dam3d",GiDPostMode.GiD_PostBinary)
gid_io.ReadMesh(model_part.GetMesh())
gid_io.WriteMesh((model_part).GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);
print model_part

g = Array3()
g[2] = -9.81;
for node in model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,0.000001)
    node.SetSolutionStepValue(DENSITY,0,1.000)
    node.SetSolutionStepValue(BODY_FORCE,0,g)

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#importing the solver files
import incompressible_fluid_solver
incompressible_fluid_solver.AddVariables(model_part)


#creating a fluid solver object
fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(model_part,domain_size)
fluid_solver.predictor_corrector = True
fluid_solver.laplacian_form =2;
fluid_solver.vel_toll = 1e-3
fluid_solver.time_order = 0

pILUPrecond = ILU0Preconditioner() 
fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pILUPrecond)


fluid_solver.Initialize()

#settings to be changed
Re = 100.0
nsteps = 100
output_step = 1


Dt = 0.001
out = 0


for step in range(0,nsteps):
    print "line49"

    time = Dt*step
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 4):
        fluid_solver.Solve()


    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        out = 0
    out = out + 1

          
        

