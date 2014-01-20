from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
# setting the domain size for the problem to be solved
domain_size = 2  # 2D problem

# including kratos path
from KratosMultiphysics import *  # we import the KRATOS
from KratosMultiphysics.PFEM2Application import *  # and now our application. note that we can import as many as we need to solve our specific problem
from KratosMultiphysics.ConvectionDiffusionApplication import *  # and now our application. note that we can import as many as we need to solve our specific problem

# defining a model part
model_part = ModelPart("ExampleModelPart")
# we create a model part
linea_model_part = ModelPart("linea")

import pfem_2_solver  # we import the python file that includes the commands that we need
# static_pressure_solver.AddVariables(model_part,linea_model_part)  #from the static_poisson_solver.py we call the function Addvariables so that the model part we have just created has the needed variables
pfem_2_solver.AddVariables(model_part, linea_model_part)

from math import sqrt
from math import sin
 # (note that our model part does not have nodes or elements yet)

 # now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file
gid_mode = GiDPostMode.GiD_PostBinary  # we import the python file that includes the commands that we need
multifile = MultiFileFlag.SingleFile  # MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("output", gid_mode, multifile, deformed_mesh_flag, write_conditions)

model_part_io = ModelPartIO("viscous_cylinder_in_air")             # we set the name of the .mdpa file
model_part_io.ReadModelPart(model_part)         # we load the info from the .mdpa
8

model_part.ProcessInfo.SetValue(DENSITY_WATER, 1000.0);
model_part.ProcessInfo.SetValue(DENSITY_AIR, 1.0);
model_part.ProcessInfo.SetValue(VISCOSITY_WATER, 100.0);
model_part.ProcessInfo.SetValue(VISCOSITY_AIR, 0.0001);
model_part.ProcessInfo.SetValue(GRAVITY_X, 0.0);
model_part.ProcessInfo.SetValue(GRAVITY_Y, -9.8);


for element in model_part.Elements:
    element.SetValue(IS_INACTIVE, False)


# the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)
model_part.SetBufferSize(2)

 # we add the DoFs
pfem_2_solver.AddDofs(model_part)

for node in model_part.Nodes:
    art = 1 + 1
    node.SetSolutionStepValue(DISTANCE, 0, 1.0)
    if ((node.Y) > 0.292 and node.X < 0.146 and node.X > 0.08):
    # if (node.Y>0.4):
        node.SetSolutionStepValue(DISTANCE, 0, -1.0)  # piscina
    if node.X > 0.583999 or node.X < 0.0001:
        node.Fix(VELOCITY_X)
        node.Fix(FRACT_VEL_X)
    if node.Y > 0.583999 or node.Y < 0.00001:
        node.Fix(VELOCITY_Y)
        node.Fix(FRACT_VEL_Y)
    # node.Fix(DISTANCE)
    #	node.SetSolutionStepValue(DISTANCE,-0.1)
    # if node.Y<0.001 or node.Y>0.583999:
    #	node.Fix(VELOCITY_X)
    if node.Y > 0.5839:
        node.Fix(VELOCITY_Y)
        node.Fix(VELOCITY_X)
        node.Fix(FRACT_VEL_Y)
        node.Fix(FRACT_VEL_X)
        # node.Fix(PRESSURE)
        # node.Fix(DISTANCE)
        # node.SetSolutionStepValue(DISTANCE,0,1.0)

    # node.SetSolutionStepValue(DISTANCE,0,node.Y-0.3) #piscina
    # if ((node.Y)<0.5 and  (node.Y)>0.4 and node.X>0.1 and node.X<0.5):
    # if (node.Y>0.5 or node.X>0.5):
    # node.SetSolutionStepValue(DISTANCE,0,-1.0) #piscina

# creating a solver object
solver = pfem_2_solver.PFEM2Solver(model_part, domain_size)
solver.time_order = 1
solver.echo_level = 4
solver.Initialize()


print_lagrangian_mesh = False;
# print_lagrangian_mesh=True;

mesh_name = 0.0
gid_io.InitializeMesh(mesh_name);
if (print_lagrangian_mesh):
    gid_io.WriteNodeMesh((linea_model_part).GetMesh());
else:
    gid_io.WriteMesh((model_part).GetMesh());

gid_io.FinalizeMesh()


if (print_lagrangian_mesh):
    gid_io.InitializeResults(mesh_name, (linea_model_part).GetMesh())
else:
    gid_io.InitializeResults(mesh_name, (model_part).GetMesh())


nsteps = 4000
Dt = 0.005
out = 0
out_step = 1

if (print_lagrangian_mesh):
    gid_io.WriteNodalResults(TEMPERATURE, linea_model_part.Nodes, 0, 0)
    gid_io.WriteNodalResults(ERASE_FLAG, linea_model_part.Nodes, 0, 0)
    gid_io.WriteNodalResults(DISPLACEMENT, linea_model_part.Nodes, 0, 0)
    gid_io.WriteNodalResults(VELOCITY, linea_model_part.Nodes, 0, 0)
else:
    print("dgsdsfhdshfdshgidhsfh")
    gid_io.WriteNodalResults(TEMPERATURE, model_part.Nodes, 0, 0)
    gid_io.WriteNodalResults(YP, model_part.Nodes, 0, 0)
    gid_io.WriteNodalResults(G_VALUE, model_part.Nodes, 0, 0)
    gid_io.WriteNodalResults(VELOCITY, model_part.Nodes, 0, 0)
    gid_io.WriteNodalResults(PRESSURE, model_part.Nodes, 0, 0)
    gid_io.WriteNodalResults(PRESS_PROJ, model_part.Nodes, 0, 0)


out_file1 = open("max_courant.dat", 'w')
out_file2 = open("mean_courant.dat", 'w')
out_file3 = open("water_fraction.dat", 'w')

calculate_water_fraction = CalculateWaterFraction(model_part)

for step in range(1, nsteps):
    out = out + 1
    print("new step")
    time = Dt * step

    model_part.CloneTimeStep(time)

    if step > 1:
        solver.Solve()

    # if step==500:
        # model_part.ProcessInfo.SetValue(GRAVITY_Y, 0.0);
        # model_part.ProcessInfo.SetValue(VISCOSITY, 0.0);

    if step == 0:
        for node in model_part.Nodes:
            node.SetSolutionStepValue(VELOCITY_X, 0, 0.0)
            node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)
    if out == out_step:
        out = 0
        print("printing a step")

        gid_io.WriteNodalResults(G_VALUE, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VELOCITY, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(PRESSURE, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(FRACT_VEL, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(ACCELERATION, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DISTANCE, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(YP, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(PRESS_PROJ, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(PRESS_PROJ_NO_RO, model_part.Nodes, time, 0)
        # gid_io.WriteNodalResults(TEMPERATURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(RHS, model_part.Nodes, time, 0)
        # gid_io.WriteNodalResults(FLAG_VARIABLE,model_part.Nodes,time,0)
        # gid_io.WriteNodalResults(NORMAL,model_part.Nodes,time,0)
        # gid_io.WriteNodalResults(DENSITY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(MESH_VELOCITY, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(PREVIOUS_ITERATION_PRESSURE, model_part.Nodes, time, 0)
        # gid_io.WriteNodalResults(FIRST_ITERATION_PRESSURE,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(DELTA_VELOCITY, model_part.Nodes, time, 0)

    # break;
    if step > 2:  # otherwise we have no data and it will die
        max_courant = calculate_water_fraction.CalculateMaxCourant()
        out_line1 = str(time - Dt) + '	' + str(max_courant) + '\n'
        out_file1.write(out_line1)
        out_file1.flush()

        water_fraction = calculate_water_fraction.Calculate()
        out_line2 = str(time - Dt) + '	' + str(water_fraction) + '\n'
        out_file2.write(out_line2)
        out_file2.flush()

        mean_courant = calculate_water_fraction.CalculateMeanCourant()
        out_line3 = str(time - Dt) + '	' + str(mean_courant) + '\n'
        out_file3.write(out_line3)
        out_file3.flush()
        # break;


gid_io.FinalizeResults()
