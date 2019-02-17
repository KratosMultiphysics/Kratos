from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import pfem_var

#
#
# setting the domain size for the problem to be solved
domain_size = pfem_var.domain_size

#
#
# ATTENTION: here the order is important

# including kratos path
# including kratos path
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.PFEMApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *


if(Kratos_Structural_Application_var.LinearSolver == "SuperLUSolver"):
    from KratosMultiphysics.ExternalSolversApplication import *

# if(Kratos_Structural_Application_var.SolverType == "ParallelSolver"):
#     from KratosMultiphysics.MKLSolversApplication import *


# defining a model part for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart")
# solid_model_part = ModelPart("StructurePart")

#
#


def NodeFinder(nodes_list, x, y, z):
    print("inside finfer")
    for node in nodes_list:
        a = (node.X - x) ** 2
        a = a + (node.Y - y) ** 2
        a = a + (node.Z - z) ** 2
        if(a < 0.000001):
            print(node)
            return node
#
#


def PrintPressure(time, filename, node):
    output = str(time) + " "
    output += str(node.GetSolutionStepValue(WATER_PRESSURE)) + " "
    output += str(node.GetSolutionStepValue(AIR_PRESSURE)) + "\n"
    filename.write(output)


#
#


#
# importing the solvers needed
SolverType = pfem_var.SolverType
if(SolverType == "FractionalStep"):
    import incompressible_fluid_solver
    incompressible_fluid_solver.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    import monolithic_solver_eulerian
    monolithic_solver_eulerian.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian_compressible"):
    import monolithic_solver_eulerian_compressible
    monolithic_solver_eulerian_compressible.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_lagrangian"):
    import monolithic_solver_lagrangian_contact
    monolithic_solver_lagrangian_contact.AddVariables(fluid_model_part)
else:
    raise "solver type not supported: options are FractionalStep - Monolithic"

# introducing input file name
input_file_name = pfem_var.problem_name

# reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = TwoFluidGidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
# gid_io = TwoFluidGidIO(input_file_name,gid_mode,multifile,deformed_mesh_flag, write_conditions)

model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

# adding dofs
if(SolverType == "FractionalStep"):
    incompressible_fluid_solver.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    monolithic_solver_eulerian.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian_compressible"):
    monolithic_solver_eulerian_compressible.AddDofs(fluid_model_part)
elif(SolverType == "monolithic_solver_lagrangian"):
    monolithic_solver_lagrangian_contact.AddDofs(fluid_model_part)


# check to ensure that no node has zero density or pressure
# for node in fluid_model_part.Nodes:
#    if((node.GetSolutionStepValue(DENSITY_WATER) == 0.0) and (node.GetSolutionStepValue(DENSITY_AIR) == 0.0)):
#        print "node ",node.Id," has zero density!"
#        raise 'node with zero density found'
#    if((node.GetSolutionStepValue(VISCOSITY_WATER) == 0.0) and (node.GetSolutionStepValue(VISCOSITY_AIR) == 0.0)):
#        print "node ",node.Id," has zero viscosity!"
#        raise 'node with zero VISCOSITY found'

mesh_name = 0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh((fluid_model_part).GetMesh())
gid_io.FinalizeMesh()

# Constitutive law
for prop in fluid_model_part.Properties:
    prop.SetValue(CONSTITUTIVE_LAW, Isotropic2D())
    # prop.SetValue(CONSTITUTIVE_LAW, PlaneStressJ2() )
    # prop.SetValue(DENSITY,2700.0000);
    # prop.SetValue(POISSON_RATIO, 0.3);
    # prop.SetValue(THICKNESS, 0.000899);
    # prop.SetValue(YIELD_STRESS, 2.76e8);
    # prop.SetValue(PLASTIC_MODULUS, 0.0); #6.34317e8);
    # prop.SetValue(YOUNG_MODULUS, 70000e6);


print(fluid_model_part)
print(fluid_model_part.Properties)

# setting the limits of the bounding box
box_corner1 = Vector(3)
box_corner1[0] = pfem_var.min_x
box_corner1[1] = pfem_var.min_y
box_corner1[2] = pfem_var.min_z
box_corner2 = Vector(3);
box_corner2[0] = pfem_var.max_x;
box_corner2[1] = pfem_var.max_y;
box_corner2[2] = pfem_var.max_z;

# time setting
output_Dt = pfem_var.output_Dt
max_dt = pfem_var.max_dt
min_dt = pfem_var.min_dt
safety_factor = pfem_var.safety_factor
nsteps = pfem_var.nsteps


# creating the solvers
# fluid solver
if(SolverType == "FractionalStep"):
    fluid_solver = incompressible_fluid_solver.IncompressibleFluidSolver(fluid_model_part, domain_size)
    fluid_solver.laplacian_form = laplacian_form;  # standard laplacian form
    fluid_solver.predictor_corrector = fluid_only_var.predictor_corrector
    fluid_solver.max_press_its = fluid_only_var.max_press_its
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_eulerian"):
    fluid_solver = monolithic_solver_eulerian.MonolithicSolver(fluid_model_part, domain_size)
    oss_swith = fluid_only_var.use_oss
    dynamic_tau = fluid_only_var.dynamic_tau
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_swith);
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.Initialize()
elif(SolverType == "monolithic_solver_lagrangian"):
    fluid_solver = monolithic_solver_lagrangian_contact.MonolithicSolver(fluid_model_part, domain_size, box_corner1, box_corner2)
    # fluid_solver.remeshing_flag = False
    oss_swith = pfem_var.use_oss
    dynamic_tau = pfem_var.dynamic_tau
    fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_swith);
    fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau);
    fluid_solver.Initialize(output_Dt)
    fluid_solver.remeshing_flag = False
    fluid_solver.contact_mesh_flag = False

#
# test_node_1 = NodeFinder(fluid_model_part.Nodes, -0.1016,0.1373, 0.0)
# print "**************************FIRST sensor*******************"
# print test_node_1.GetSolutionStepValue(FLAG_VARIABLE)
# pressureout_1 = open("pressure_history_1.out", 'w')

# test_node_2 = NodeFinder(fluid_model_part.Nodes, 0.1016,0.1373, 0.0)
# print "**************************second sensor*******************"
# print test_node_2.GetSolutionStepValue(FLAG_VARIABLE)
# pressureout_2 = open("pressure_history_2.out", 'w')

# test_node_3 = NodeFinder(fluid_model_part.Nodes, 0.0,0.1373, -0.1016)
# print "**************************third sensor*******************"
# print test_node_3.GetSolutionStepValue(FLAG_VARIABLE)
# pressureout_3 = open("pressure_history_3.out", 'w')

# test_node_4 = NodeFinder(fluid_model_part.Nodes, 0.0,0.1373, 0.1016)
# print "**************************forth sensor*******************"
# print test_node_4.GetSolutionStepValue(FLAG_VARIABLE)
# pressureout_4 = open("pressure_history_4.out", 'w')

#

for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(FORCE_X, 0, 0.0)
    node.SetSolutionStepValue(FORCE_Y, 0, 10000.0)
    node.SetSolutionStepValue(FORCE_Z, 0, 0.0)

fluid_solver.output_time_increment = 0.00000000005
time = 0.0
new_Dt = max_dt
for step in range(0, nsteps):
    print("line49")

   # new_Dt = fluid_solver.EstimateDeltaTime(min_dt,max_dt)
    # new_Dt = max_dt

    time = fluid_model_part.ProcessInfo[TIME]

    time = time + new_Dt * safety_factor
    print(safety_factor)
    fluid_model_part.CloneTimeStep(time)

    print("new_time= ", time)

    # solving the fluid problem
    if(step > 3):
        fluid_solver.Solve(time, gid_io)
        # PrintPressure(time,pressureout_1,test_node_1)
        # pressureout_1.flush()
        # PrintPressure(time,pressureout_2,test_node_2)
        # pressureout_2.flush()
        # PrintPressure(time,pressureout_3,test_node_3)
        # pressureout_3.flush()
        # PrintPressure(time,pressureout_4,test_node_4)
        # pressureout_4.flush()
    new_Dt = fluid_model_part.ProcessInfo[DELTA_TIME]
