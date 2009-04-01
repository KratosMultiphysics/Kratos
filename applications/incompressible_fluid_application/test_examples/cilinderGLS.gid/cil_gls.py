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
#applications_interface.Import_ULFApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################

from KratosIncompressibleFluidApplication import *


#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here when the "very fix container will be ready"
import runge_kutta_frac_step_solver
runge_kutta_frac_step_solver.AddVariables(model_part)

model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT)

#reading a model
gid_mode = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("cilinderGLS",gid_mode,use_multifile,deformed_print_flag,write_conditions)
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
runge_kutta_frac_step_solver.AddDofs(model_part)

gravity = Vector(3);
gravity[0] = 0.0; gravity[1] = -1.0; gravity[2] = 0.0
zero = Vector(3);
zero[0] = 0.0; zero[1] = 0.0; zero[2] = 0.0
##for node in model_part.Nodes:    
##    if(node.Y > 0.99 and node.X > 0.001 and node.X < 0.999):
##        node.Free(VELOCITY_X);
##        node.Free(VELOCITY_Y);
##        #node.Free(VELOCITY_Z);
##    #node.Free(PRESSURE)
##    node.SetSolutionStepValue(BODY_FORCE,0,gravity);
##    node.SetSolutionStepValue(VELOCITY,0,zero);
        
#creating a fluid solver object
fluid_solver = runge_kutta_frac_step_solver.RungeKuttaFracStepSolver(model_part,domain_size)

##pILUPrecond = ILU0Preconditioner() 
##fluid_solver.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)

fluid_solver.Initialize()

#settings to be changed
Re = 100.0
nsteps = 5000
output_step = 50

for node in model_part.Nodes :
    node.SetSolutionStepValue(VISCOSITY,0,1.0/Re);
    #node.SetSolutionStepValue(VISCOSITY,0,0.000017);
    node.SetSolutionStepValue(DENSITY,0,1.0);
    #node.SetSolutionStepValue(BODY_FORCE_Y,0,-10.00)

    #node.SetSolutionStepValue(PRESSURE,0,(1.0-node.Y)*10.0)
    #node.SetSolutionStepValue(VELOCITY_X,0,-1.0)
    #node.SetSolutionStepValue(VELOCITY_Y,0,-1.0)
    #node.SetSolutionStepValue(BODY_FORCE_Y,0,-50.0)
    #node.SetSolutionStepValue(BODY_FORCE_X,0,0.0)

#now we compute the delta time using CFL law
CFL_time_estimate_process=CFLProcess(model_part)
CFL=0.8;

#value to initialize only
Dt = 0.0;
dt_max=0.001;
out = 0
time=0.0
gid_io.InitializeResults( 0.0, model_part.GetMesh() )

for step in range(0,nsteps):
    print "line49"
    Dt=(CFL_time_estimate_process).EstimateTime(CFL, dt_max)
    print "CFL gave this time step", Dt
           
    time = time + Dt
    model_part.CloneTimeStep(time)

    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem
    if(step > 3):
        fluid_solver.Solve()

    #print the results
    if(out == output_step):
        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(ACCELERATION,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(BODY_FORCE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISPLACEMENT,model_part.Nodes,time,0)

        out = 0
    out = out + 1

          
gid_io.FinalizeResults();
        

