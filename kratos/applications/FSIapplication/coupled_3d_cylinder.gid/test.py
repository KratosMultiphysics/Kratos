#path for applications s
import sys
applications_path = '/home/rrossi/kratosR1/applications/'
sys.path.append(applications_path);

#importing the Kratos Library
from Kratos import *

sys.path.append('/home/rrossi/kratosR1/applications/incompressible_fluid_application/python_scripts')
sys.path.append('/home/rrossi/kratosR1/applications/incompressible_fluid_application/Linux')
from KratosIncompressibleFluidApplication import *

sys.path.append(applications_path + 'ALEapplication/python_scripts' )
sys.path.append(applications_path + 'ALEapplication/Linux')
from KratosALEApplication import *


kernel = Kernel()
print Hello()

#############################################
#defining the domain size
domain_size = 3

#importing the applications library and adding them to the kernel
incompressible_fluid_application = KratosIncompressibleFluidApplication()
kernel.AddApplication(incompressible_fluid_application)

ale_app = KratosALEApplication();
kernel.AddApplication(ale_app)

#dynamic renumbering of variables to ensure consistency 
kernel.Initialize()
kernel.InitializeApplication(incompressible_fluid_application);
kernel.InitializeApplication(ale_app);

#defining a model part
fluid_model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here when the "very fix container will be ready"

#reading a model
gid_io = GidIO("coupled_3d_cylinder_fluid",GiDPostMode.GiD_PostBinary)
gid_io.ReadMesh(fluid_model_part.GetMesh())
gid_io.WriteMesh((fluid_model_part).GetMesh(),domain_size,GiDPostMode.GiD_PostBinary);
print fluid_model_part

#the buffer size should be set up here after the mesh is read for the first time
fluid_model_part.SetBufferSize(3)

#importing the solver files
import incompressible_fluid_solver
incompressible_fluid_solver.AddVariables(fluid_model_part)


#creating a fluid solver object
fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(fluid_model_part,domain_size)
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
for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,viscosity);
    node.SetSolutionStepValue(DENSITY,0,density);

for step in range(0,nsteps):
    print "line49"

    time = Dt*step
    fluid_model_part.CloneTimeStep(time)

    print time
    #print fluid_model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 4):
        fluid_solver.Solve()


    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(PRESSURE,fluid_model_part.Nodes,time,0)
        gid_io.WriteNodalResults(VELOCITY,fluid_model_part.Nodes,time,0)
        out = 0
    out = out + 1

          
        

