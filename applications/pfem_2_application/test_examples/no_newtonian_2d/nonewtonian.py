##################################################################
##################################################################
#setting the domain size for the problem to be solved
domain_size = 2  # 2D problem

#including kratos path
from KratosMultiphysics import *    #we import the KRATOS
from KratosMultiphysics.PFEM2Application import *        #and now our application. note that we can import as many as we need to solve our specific problem
from KratosMultiphysics.ExternalSolversApplication import *
#defining a model part
model = Model()
model_part = model.CreateModelPart("ExampleModelPart");  #we create a model part

import pfem_2_solver_monolithic_fluid as pfem2_solver           #we import the python file that includes the commands that we need
#static_pressure_solver.AddVariables(model_part,linea_model_part)  #from the static_poisson_solver.py we call the function Addvariables so that the model part we have just created has the needed variables
pfem2_solver.AddVariables(model_part)

from math import sqrt
from math import sin
from math import cos
from math import fabs
 # (note that our model part does not have nodes or elements yet)

 #now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file
gid_mode = GiDPostMode.GiD_PostBinary  #we import the python file that includes the commands that we need
multifile = MultiFileFlag.SingleFile #MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("nonewtonian_test",gid_mode,multifile,deformed_mesh_flag,write_conditions)


model_part_io = ModelPartIO("nonewtonian")             # we set the name of the .mdpa file
model_part_io.ReadModelPart(model_part)         # we load the info from the .mdpa


gravity = -9.8
model_part.ProcessInfo.SetValue(GRAVITY_X, 0.0 );
model_part.ProcessInfo.SetValue(GRAVITY_Y, -9.8)
model_part.ProcessInfo.SetValue(GRAVITY_Z, 0.0 );

model_part.ProcessInfo.SetValue(VISCOSITY_WATER, 0.001 );
model_part.ProcessInfo.SetValue(DENSITY_WATER, 1600.0);
model_part.ProcessInfo.SetValue(VISCOSITY_AIR, 0.001 );
model_part.ProcessInfo.SetValue(DENSITY_AIR, 1.0);


#the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)
model_part.SetBufferSize(3)
list_nodes = [] #here we create an empty list

pfem2_solver.AddDofs(model_part)

node_id=0


for node in model_part.Nodes:
        dist=(-1.0)
        if node.X>node.Y+0.1:
            dist=node.X-0.2+0.0
        else:
            dist= node.Y-0.1+0.0
        node.SetSolutionStepValue(DISTANCE,0,dist)

        if node.Y<0.001 or node.Y>0.14999:
              node.Fix(VELOCITY_Y)
              node.Fix(VELOCITY_X)
        if node.X>0.49999 or node.X<0.0001:
             node.Fix(VELOCITY_X)
        if node.Y>0.1499:
             node.Fix(PRESSURE)

        node.SetSolutionStepValue(BODY_FORCE_Y,0,-9.8)


#creating a solver object
maximum_nonlin_iterations=10 #linear problem
solver = pfem2_solver.PFEM2Solver(model_part,domain_size,maximum_nonlin_iterations)
solver.time_order = 1
solver.echo_level = 1
solver.Initialize()

mesh_name = 0.0
gid_io.InitializeMesh( mesh_name );
gid_io.WriteMesh((model_part).GetMesh());
gid_io.FinalizeMesh()




gid_io.InitializeResults(mesh_name,(model_part).GetMesh())


nsteps=202
Dt=0.005
out=0
out_step=1

gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,0,0)
gid_io.WriteNodalResults(DISTANCE,model_part.Nodes,0,0)
gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,0,0)

time=0.0
for step in range(1,nsteps):
    out=out+1
    print("---NEW STEP---")
    #if time>0.3:
    #	Dt=0.001
    time = time+Dt
    #time = Dt*step
    model_part.CloneTimeStep(time)

    despl=0.0;
    if step>1:
        solver.Solve()

    if out==out_step:
       out=0
       print("printing a step")
       gid_io.WriteNodalResults(VELOCITY,model_part.Nodes,time,0)
       gid_io.WriteNodalResults(PRESSURE,model_part.Nodes,time,0)
       gid_io.WriteNodalResults(DISTANCE,model_part.Nodes,time,0)



gid_io.FinalizeResults()
