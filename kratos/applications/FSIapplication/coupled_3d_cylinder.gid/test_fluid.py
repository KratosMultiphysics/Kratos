#path for applications s
import sys
sys.path.append('/home/rrossi/kratosR1/applications/incompressible_fluid_application/python_scripts')
sys.path.append('/home/rrossi/kratosR1/applications/incompressible_fluid_application/Linux')

#importing the Kratos Library
from Kratos import *
from KratosIncompressibleFluidApplication import *


kernel = Kernel()
print Hello()

#############################################
#defining the domain size
domain_size = 3

#importing the applications library and adding them to the kernel
incompressible_fluid_application = KratosIncompressibleFluidApplication()
kernel.AddApplication(incompressible_fluid_application)

#dynamic renumbering of variables to ensure consistency 
kernel.Initialize()
kernel.InitializeApplication(incompressible_fluid_application);

#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here when the "very fix container will be ready"

#reading a model
gid_io = GidIO("coupled_3d_cylinder_fluid",GiDPostMode.GiD_PostBinary)
gid_io.ReadMesh(model_part.GetMesh())
gid_io.WriteMesh((model_part).GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);
print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(3)

#importing the solver files
import incompressible_fluid_solver
incompressible_fluid_solver.AddVariables(model_part)


#creating a fluid solver object
fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(model_part,domain_size)
##fluid_solver.laplacian_form =2;
##fluid_solver.vel_toll = 1e-3
##fluid_solver.time_order = 2
##
##pILUPrecond = ILU0Preconditioner() 
##fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)


fluid_solver.Initialize()

#settings to be changed
Re = 100.0
nsteps = 100
output_step = 20

Dt = 1000.0/Re
if(Dt > 0.1):
    Dt = 0.1
out = 0

viscosity = 0.00001;
density = 1.21;
for node in model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,viscosity);
    node.SetSolutionStepValue(DENSITY,0,density);

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

          
        

