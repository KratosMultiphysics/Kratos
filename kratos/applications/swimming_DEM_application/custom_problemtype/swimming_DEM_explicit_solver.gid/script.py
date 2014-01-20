from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
import math
import matplotlib
from numpy import *
from pylab import *

# Importing the problem parameters
import UD_var as Param

# Including kratos path
sys.path.append(Param.kratos_path)

from KratosMultiphysics import *
from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.PFEMApplication import *
from KratosMultiphysics.StructuralApplication import *

import DEM_procedures as DEMProc
import swimming_DEM_procedures as SwimProc
import swimming_sphere_strategy as SwimStrat

my_timer = Timer()

# Constructing a DEM_procedures object
DEM_proc = DEMProc.Procedures(Param)

# Constructing an IOTools object
IO_tools = SwimProc.IOTools(Param)

# Constructing fluid-DEM strategy (the constructor automatically imports a particle strategy and a fluid strategy):
solver_strategy = SwimStrat.ULFDEMStrategy(Param)

# Defining model parts for the fluid and one for the structure
fluid_model_part = ModelPart("FluidPart")
structure_model_part = ModelPart("StructurePart")
combined_model_part = ModelPart("CombinedPart")
mixed_model_part = ModelPart("MixedPart")

if (Param.SolverType == "Incompressible_Modified_FracStep" or Param.SolverType == "FracStep"):
    fluid_only_model_part = ModelPart("FluidOnlyPart")

# Defining a model part for the DEM part
balls_model_part = ModelPart("SolidPart")
solver_strategy.AddVariables(balls_model_part, fluid_model_part)

# Reading the solid part: binary or ascii, multifile or single

if (Param.OutputFileType == "Binary"):
    gid_mode = GiDPostMode.GiD_PostBinary

else:
    gid_mode = GiDPostMode.GiD_PostAscii

if (Param.Multifile == "multiple_files"):
    DEM_multifile = MultiFileFlag.MultipleFiles

else:
    DEM_multifile = MultiFileFlag.SingleFile

deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
DEM_gid_io = GidIO(Param.DEM_problem_name, gid_mode, DEM_multifile, deformed_mesh_flag, write_conditions)

# Introducing input file name
input_file_name = Param.problem_name

# Reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions = WriteConditionsFlag.WriteConditions
gid_io = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_origin = ModelPartIO(input_file_name)

# Reading fluid model part
model_part_io_origin.ReadModelPart(fluid_model_part)

print("ULF model read correctly")

# Setting the fluid buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

# Reading the DEM model part
model_part_io_solid = ModelPartIO(Param.DEM_problem_name)
model_part_io_solid.ReadModelPart(balls_model_part)

print("DEM model read correctly")

# Setting the DEM buffer size: SHOULD BE DONE AFTER READING!!!
balls_model_part.SetBufferSize(3)

# Adding dofs to the nodes of each model part
solver_strategy.AddDofs(balls_model_part, fluid_model_part)

if (Param.SolverType == "Quasi_Inc_Constant_Pressure" or Param.SolverType == "Quasi_Inc_Linear_Pressure"):

    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

# Setting the limits of the bounding box
box_corner_1 = Vector(3)
box_corner_1[0] = Param.bounding_box_corner1_x
box_corner_1[1] = Param.bounding_box_corner1_y
box_corner_1[2] = Param.bounding_box_corner1_z

box_corner_2 = Vector(3)
box_corner_2[0] = Param.bounding_box_corner2_x
box_corner_2[1] = Param.bounding_box_corner2_y
box_corner_2[2] = Param.bounding_box_corner2_z

# Here we write the convergence data...
outstring2 = "convergence_info.txt"
outputfile1 = open(outstring2, 'w')

# Creating solvers

# Fluid solver
fluid_solver = solver_strategy.FluidStrategy(outputfile1, fluid_only_model_part, fluid_model_part, structure_model_part, combined_model_part, box_corner_1, box_corner_2)
fluid_solver.Initialize()

print("Fluid solver created and initialized")

# DEM_solver
DEM_solver = solver_strategy.DEMStrategy(balls_model_part)  # here, solver variables initialize as default
DEM_proc.GiDSolverTransfer(balls_model_part, DEM_solver)
DEM_solver.Initialize()

print("DEM solver created and initialized")

# Checking to ensure that no node has zero density or pressure
is_fsi_interf = 0.0
[inverted_elements, domain_volume] = fluid_solver.CheckForInvertedElements()

# Calculating porosity
porosity_utils = DEMProc.GranulometryUtils(domain_volume, balls_model_part)
porosity_utils.PrintCurrentData()
n_particles_in_depth = int(math.sqrt(porosity_utils.number_of_balls / domain_volume))

# Constitutive laws

if (Param.FSI):

    if (Param.domain_size == 2):
        fluid_model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D())

    elif (Param.domain_size == 3):
        fluid_model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D())

    else:
        raise "Domain size error. It should be 2D or 3D"

for node in fluid_model_part.Nodes:

    if (node.GetSolutionStepValue(DENSITY) == 0.0):
        print("node ", node.Id, " has zero density!")
        raise 'node with zero density found'

    if (node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print("node ", node.Id, " has zero viscosity!")
        raise 'node with zero VISCOSITY found'

    if (Param.FSI):
        is_fsi_interf += node.GetSolutionStepValue(IS_INTERFACE)

if (Param.SolverType == "Incompressible_Modified_FracStep" and Param.FSI):

    if (is_fsi_interf == 0):
        raise 'For running FSI using the Modified Frac Step Solver you must prescribe IS_INTERFACE flag at the surface/outer contour of your structure'

# Adding dofs
solver_strategy.AddDofs(balls_model_part, fluid_model_part)

# Choosing names for directories to be created
dir_names = []
dir_names.append('post_files')
dir_names.append('post_lists')
dir_names.append('neigh_lists')
dir_names.append('data_and_results')
dir_names.append('graphs')
dir_names.append('MPI_results')
dir_names.append('fluid_results')
dir_names.append('mixed_results')
main_path = os.getcwd()
directories = IO_tools.CreateProblemDirectories(main_path, dir_names)

# Creating a variable for each directory name 'name' with the name 'name_path'

for i in range(len(dir_names)):
    vars()[dir_names[i] + '_path'] = directories[i]

os.chdir(data_and_results_path)
results = open('results.txt', 'w')
summary_results = open('summary_results.txt', 'w')

force_list = []
force_list_2 = []
time_list = []

os.chdir(post_lists_path)
DEM_multifile = open(Param.DEM_problem_name + '_all' + '.post.lst', 'w')
multifile_5 = open(Param.DEM_problem_name + '_5' + '.post.lst', 'w')
multifile_10 = open(Param.DEM_problem_name + '_10' + '.post.lst', 'w')
multifile_50 = open(Param.DEM_problem_name + '_50' + '.post.lst', 'w')

DEM_multifile.write('Multiple\n')
multifile_5.write('Multiple\n')
multifile_10.write('Multiple\n')
multifile_50.write('Multiple\n')

first_print = True
index_5 = 1
index_10 = 1
index_50 = 1
prev_time = 0.0
control = 0.0

os.chdir(post_files_path)

if (Param.Multifile == "single_file"):
    DEM_gid_io.InitializeMesh(0.0)
    DEM_gid_io.WriteSphereMesh(balls_model_part.GetMesh())
    DEM_gid_io.FinalizeMesh()
    DEM_gid_io.InitializeResults(0.0, balls_model_part.GetMesh())

os.chdir(main_path)

# For plotting the graph:
velocity_node_y = 0.0

# Initializations
initial_dt = 0.001 * Param.Dt
safety_factor = 0.5
time = 0.0
step = 0

DEM_dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME)
DEM_time = 0.0
DEM_step = 0
DEM_time_old_print = 0.0

initial_pr_time = timer.clock()
initial_real_time = timer.time()

print(('\n' + 'Calculation starts at instant: ' + str(initial_pr_time) + '\n'))

total_steps_expected = int(Param.FinalTime / DEM_dt)

print(('Total number of TIME STEPs expected in the calculation is: ' + str(total_steps_expected) + ' if time step is kept ' + '\n'))

inlet_vel = Vector(3)

if (Param.lagrangian_nodes_inlet == 1):

    for node in fluid_model_part.Nodes:

        if (node.GetSolutionStepValue(IS_LAGRANGIAN_INLET)):

            inlet_vel = node.GetSolutionStepValue(VELOCITY, 0)
            print("Lagrangian inlet(s) velocity  is ", inlet_vel)
            break
else:
    inlet_vel[0] = 0.0
    inlet_vel[1] = 0.0
    inlet_vel[2] = 0.0

inlet_process = LagrangianInletProcess(fluid_model_part, 0.0, inlet_vel)

# Creation of projection module and initial projection
h_min = 0.01
projection_module = SwimProc.ProjectionModule(fluid_model_part, balls_model_part, Param.domain_size, n_particles_in_depth)
projection_module.UpdateDatabase(h_min)
projection_module.ProjectFromFluid()

# Temporal loop

while (time < Param.max_time):
    DEM_dt = balls_model_part.ProcessInfo.GetValue(DELTA_TIME)
    DEM_time = DEM_time + DEM_dt

    if (DEM_time >= time):
        step += 1

        if (step <= 3):
            new_dt = 0.0000001
            time = time + new_dt * safety_factor

        # Solving the fluid problem

        if (step > 3):
            new_dt = fluid_solver.EstimateDeltaTime(Param.Dt, Param.domain_size)
            time = time + new_dt * safety_factor
            combined_model_part.CloneTimeStep(time)

            # Updating containers database
            projection_module.UpdateDatabase(h_min)

            # Coupling DEM to fluid
            projection_module.ProjectFromParticles()

            # Solving fluid
            fluid_solver.Solve(inlet_process)

            # Coupling fluid to DEM
            projection_module.ProjectFromFluid()

    # Solving particles
    DEM_solver.Solve()

    incremental_time = (timer.time() - initial_real_time) - prev_time
    balls_model_part.CloneTimeStep(DEM_time)
    balls_model_part.ProcessInfo[TIME_STEPS] = DEM_step
    total_force = 0
    force_node = 0

    # Writing lists to be printed
    os.chdir(data_and_results_path)
    force_list.append(total_force)
    time_list.append(DEM_time)

    os.chdir(main_path)
    total_force = 0
    force_node = 0

    # Echoes
    IO_tools.ControlEcho(step, incremental_time, total_steps_expected)

    if (first_print and incremental_time > 60.0):
        IO_tools.CalculationLengthEstimationEcho(step, incremental_time, total_steps_expected)
        first_print = False

    # Writting
    os.chdir(post_lists_path)
    DEM_multifile.write(Param.DEM_problem_name + '_' + str(DEM_time) + '.post.bin\n')
    os.chdir(main_path)

    DEM_step += 1

  # GiD IO        ################################################################################

    time_to_print = DEM_time - DEM_time_old_print

    if (time_to_print >= Param.OutputTimeStep):
        os.chdir(main_path)

        if (Param.PrintNeighbourLists == "ON"):  # Printing neighbours id's
            os.chdir(neigh_lists_path)
            neighbours_list = open('neigh_list_' + str(DEM_time), 'w')

            for elem in balls_model_part.Elements:
                ID = (elem.Id)
                Neigh_ID = elem.GetValue(NEIGHBOURS_IDS)

                for i in range(len(Neigh_ID)):
                    neighbours_list.write(str(ID) + ' ' + str(Neigh_ID[i]) + '\n')

            neighbours_list.close()

        if (Param.Multifile == "multiple_files"):
            os.chdir(mixed_results_path)
            ParticleUtils2D().VisualizationModelPart(mixed_model_part, fluid_model_part, balls_model_part)
            DEM_gid_io.InitializeMesh(time)
            DEM_gid_io.WriteSphereMesh(balls_model_part.GetMesh())
            DEM_gid_io.WriteMesh(mixed_model_part.GetMesh())
            DEM_gid_io.FinalizeMesh()
            DEM_gid_io.InitializeResults(10.0, mixed_model_part.GetMesh())
            DEM_gid_io.WriteNodalResults(VELOCITY, mixed_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(DISPLACEMENT, mixed_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(PRESSURE, fluid_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(PRESSURE_GRADIENT, fluid_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(DRAG_REACTION, fluid_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(RADIUS, balls_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(FLUID_VEL_PROJECTED, balls_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(DRAG_FORCE, balls_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(BUOYANCY, balls_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(PRESSURE_GRAD_PROJECTED, balls_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(TOTAL_FORCES, balls_model_part.Nodes, time, 0)

            if (Param.PostExportId == "1"):
                DEM_gid_io.WriteNodalResults(EXPORT_ID, balls_model_part.Nodes, time, 0)

            DEM_gid_io.Flush()
            DEM_gid_io.FinalizeResults()

        # PRINTING VARIABLES############

        os.chdir(data_and_results_path)

        if (index_5 == 5):
            multifile_5.write(Param.DEM_problem_name + '_' + str(DEM_time) + '.post.bin\n')
            index_5 = 0

        if (index_10 == 10):
            multifile_10.write(Param.DEM_problem_name + '_' + str(DEM_time) + '.post.bin\n')
            index_10 = 0

        if (index_50 == 50):
            multifile_50.write(Param.DEM_problem_name + '_' + str(DEM_time) + '.post.bin\n')
            index_50 = 0

        index_5 += 1
        index_10 += 1
        index_50 += 1

        if (Param.Multifile == "multiple_files"):
            DEM_gid_io.FinalizeResults()

        os.chdir(main_path)
        DEM_time_old_print = DEM_time

    # End of print loop
    os.chdir(main_path)

# End of temporal loop

if (Param.Multifile == "single_file"):
    DEM_gid_io.FinalizeResults()

os.chdir(data_and_results_path)
results.close()
summary_results.close()

os.chdir(post_lists_path)
DEM_multifile.close()
multifile_5.close()
multifile_10.close()
multifile_50.close()

os.chdir(main_path)
elapsed_pr_time = timer.clock() - initial_pr_time
elapsed_real_time = timer.time() - initial_real_time
print('Calculation ends at instant: ' + str(timer.time()))
print('Calculation ends at processing time instant: ' + str(timer.clock()))
print('Elapsed processing time: ' + str(elapsed_pr_time))
print('Elapsed real time: ' + str(elapsed_real_time))
print (my_timer)
print("ANALYSIS COMPLETED")
