##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2  # 2D problem

#including kratos path
from KratosMultiphysics import *    #we import the KRATOS
from KratosMultiphysics.PFEM2Application import *        #and now our application. note that we can import as many as we need to solve our specific problem
from KratosMultiphysics.ExternalSolversApplication import *        #and now our application. note that we can import as many as we need to solve our specific problem
#from KratosMultiphysics.OpenCLApplication import *        #in case you want to use the gpu to solve the system

#defining a model part
model = Model()
model_part = model.CreateModelPart("ExampleModelPart");  #we create a model part

import pfem_2_solver_monolithic_fluid as pfem_2_solver           #we import the python file that includes the commands that we need

pfem_2_solver.AddVariables(model_part)

from math import sqrt
from math import sin
from math import cos
 #now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file
gid_mode = GiDPostMode.GiD_PostBinary  #we import the python file that includes the commands that we need
multifile = MultiFileFlag.SingleFile #MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("results_monolithic",gid_mode,multifile,deformed_mesh_flag,write_conditions)

model_part_io = ModelPartIO("obstacle_monolithic")             # we set the name of the .mdpa file
model_part_io.ReadModelPart(model_part)         # we load the info from the .mdpa


model_part.ProcessInfo.SetValue(DENSITY_WATER, 1000.0);
model_part.ProcessInfo.SetValue(DENSITY_AIR, 1.0);
model_part.ProcessInfo.SetValue(VISCOSITY, 0.0001);
model_part.ProcessInfo.SetValue(VISCOSITY_AIR, 0.000001);
model_part.ProcessInfo.SetValue(VISCOSITY_WATER, 0.001);
model_part.ProcessInfo.SetValue(GRAVITY_X, 0.0);
model_part.ProcessInfo.SetValue(GRAVITY_Y, -9.8);

pi=3.14159

for node in model_part.Nodes:
    node.Free(VELOCITY_X)
    node.Free(VELOCITY_Y)
    node.SetSolutionStepValue(DISTANCE,0,1.0)
    if ((node.Y)<0.292 and node.X<0.146):
        node.SetSolutionStepValue(DISTANCE,0,-1.0) #piscina
    if node.X>0.583999 or node.X<0.0001:
        node.Fix(VELOCITY_X)
        # node.Fix(DISTANCE)
        # node.SetSolutionStepValue(DISTANCE,-0.1)
    if node.Y<0.001 or node.Y>0.583999:
        node.Fix(VELOCITY_Y)

    if node.Y>0.499:
        node.Fix(VELOCITY_Y)
        node.Fix(DISTANCE)
        node.SetSolutionStepValue(DISTANCE,0,1.0)

    if node.Y>0.57:
        node.Fix(DISTANCE)

    if node.X>0.2919 and node.X<0.2921 and node.Y<0.0481:
        node.Fix(VELOCITY_X)
    if node.X>0.3159 and node.X<0.3161 and node.Y<0.0481:
        node.Fix(VELOCITY_X)
    if node.X>0.2919 and node.X<0.3161 and node.Y<0.0481 and node.Y>0.0479:
        node.Fix(VELOCITY_Y)
    node.SetSolutionStepValue(BODY_FORCE_Y,0,-9.8)
    node.SetSolutionStepValue(PRESS_PROJ_Y,0,-9.8)






#the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)
model_part.SetBufferSize(2)

 # we add the DoFs
pfem_2_solver.AddDofs(model_part)


#creating a solver object
maximum_nonlin_iterations=1 #nonlinear problem
solver = pfem_2_solver.PFEM2Solver(model_part,domain_size,maximum_nonlin_iterations)
solver.time_order = 1
solver.echo_level = 3
solver.Initialize()

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());

gid_io.FinalizeMesh()


gid_io.InitializeResults(mesh_name,(model_part).GetMesh())


nsteps=200
Dt=0.005
out=0
out_step=2


gid_io.WriteNodalResults(DISTANCE,model_part.Nodes,0,0)
gid_io.WriteNodalResults(YP,model_part.Nodes,0,0)
gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,0,0)
gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,0,0)

import time as timer
t1 = timer.time()
for step in range(1,nsteps):
    out=out+1
    print("---------------")
    print("   new step")
    time = Dt*step

    model_part.CloneTimeStep(time)

    if step>1:
            solver.Solve()

    #if step==500:
        #model_part.ProcessInfo.SetValue(GRAVITY_Y, 0.0);
        #model_part.ProcessInfo.SetValue(VISCOSITY, 0.0);

    if step==0:
        for node in model_part.Nodes:
            node.SetSolutionStepValue(VELOCITY_X,0,0.0)
            node.SetSolutionStepValue(VELOCITY_Y,0,0.0)
    if out==out_step:
        out=0
        print("printing a step")

        gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DISTANCE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PRESS_PROJ,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(PROJECTED_VELOCITY,model_part.Nodes,time,0)
        gid_io.Flush()
    print("Time in simulation: ",time," s")


t2=timer.time()
total_time=t2-t1

print ("total_time", total_time)

gid_io.FinalizeResults()
