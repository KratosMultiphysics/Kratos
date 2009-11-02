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
applications_interface.Import_PFEMApplication = True
applications_interface.Import_MeshingApplication = True
applications_interface.ImportApplications(kernel, kratos_applications_path)

## from now on the order is not anymore crucial
##################################################################
##################################################################

import math

#defining a model part
model_part = ModelPart("FluidPart");  

#adding of Variables to Model Part should be here when the "very fix container will be ready"
import pfem_solver_ale
pfem_solver_ale.AddVariables(model_part)

#reading a model
gid_mode_flag = GiDPostMode.GiD_PostBinary
use_multifile = MultiFileFlag.MultipleFiles
deformed_print_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO("dam3d",gid_mode_flag,use_multifile,deformed_print_flag,write_conditions)
gid_io.ReadModelPart(model_part)

print model_part

#the buffer size should be set up here after the mesh is read for the first time
model_part.SetBufferSize(2)

#importing the solver files
pfem_solver_ale.AddDofs(model_part)

#for node in model_part.Nodes:
#    node.GetSolutionStepValue(TEMPERATURE);

#setting the limits of the bounding box
box_corner1 = Vector(3); box_corner1[0]=-0.01; box_corner1[1]=-0.01; box_corner1[2]=-0.01;
box_corner2 = Vector(3); box_corner2[0]=1.01; box_corner2[1]=1.01;  box_corner2[2]=1.01;

#creating a fluid solver object
name = str("dam3d")

solver = pfem_solver_ale.PFEMSolver(model_part,name,box_corner1,box_corner2,domain_size)
solver.predictor_corrector = True
solver.prediction_order = 1;
solver.max_vel_its = 5;
solver.max_press_its = 3;
solver.set_echo_level = 0;
#pDiagPrecond = DiagonalPreconditioner()
#solver.pressure_linear_solver =  CGSolver(1e-4, 5000,pDiagPrecond)
#solver.pressure_linear_solver = SkylineLUFactorizationSolver()

g = Array3()
g[2] = -9.81;
for node in model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY,0,0.000001)
    node.SetSolutionStepValue(DENSITY,0,1000.000)
    node.SetSolutionStepValue(BODY_FORCE,0,g)
    
Dt = 0.005
nsteps = 6
output_Dt = 0.1
min_dt = Dt
max_dt = 10*Dt
tmax = 20.00

#initializing the solver
solver.Initialize(Dt,output_Dt)

scalarout = open("pressure_history.out", 'w')
scalarout.write( "time (0.0,0.093) (0.0,0.299) \n")


time = 0.0
step = 0
new_Dt = 0.00
#for step in range(0,nsteps):
while (time < tmax):
    step = step + 1

    new_Dt = solver.EstimateDeltaTime(min_dt,max_dt)
    if(step <= 3):
        new_Dt = 0.000001
 #   step = step + 1
    time = time + new_Dt
    print "time = ", time, " new_Dt= ",new_Dt," step = ", step

    

    #time = Dt*step
    model_part.CloneTimeStep(time)

   
    print time
    #print model_part.ProcessInfo()[TIME]

    #solving the fluid problem

##    solver.Remesh()
##    solver.OutputStep(time,gid_io);

    print " step = ",step
    if(step > 3):
        solver.Solve(time,gid_io)
##    else:
##        print "remeshing"
##        print "printing"
##        solver.OutputStep(time,gid_io)


##
##    #chapuza to identify the effect of the mesh motion
##    solver.close_result_file = False
##    for node in model_part.Nodes:
##        vel = node.GetSolutionStepValue(VELOCITY);
##        mesh_vel = node.GetSolutionStepValue(MESH_VELOCITY);
##        temp = (vel[0]-mesh_vel[0])**2 + (vel[1]-mesh_vel[1])**2 + (vel[2]-mesh_vel[2])**2
##        denom = (vel[0])**2 + (vel[1])**2 + (vel[2])**2
##
##        if (denom < 0.000001):
##            node.SetSolutionStepValue(TEMPERATURE,0,0.00);
##        else:
##            node.SetSolutionStepValue(TEMPERATURE,0,temp/denom);
##
##    gid_io.WriteNodalResults(TEMPERATURE, model_part.Nodes, time, 0);
##    gid_io.Flush()
##    gid_io.CloseResultFile();



          
        

