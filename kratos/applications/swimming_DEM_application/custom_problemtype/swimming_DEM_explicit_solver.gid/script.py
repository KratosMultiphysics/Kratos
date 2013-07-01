import UD_var

import time as timer
import os
import sys
import math
import matplotlib
from numpy import *
from pylab import *

# Including kratos path
sys.path.append(UD_var.kratos_path)

from KratosMultiphysics import *
from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.PFEMApplication import *
from KratosMultiphysics.StructuralApplication import *

from DEM_explicit_solver_var import *
from DEM_procedures import *
from swimming_DEM_procedures import *

my_timer = Timer()

import swimming_sphere_strategy as solver_strategy

# Defining model parts for the fluid and one for the structure

fluid_model_part = ModelPart("FluidPart")
structure_model_part = ModelPart("StructurePart")
combined_model_part = ModelPart("CombinedPart")
mixed_model_part = ModelPart("MixedPart")

if (UD_var.SolverType == "Incompressible_Modified_FracStep" or UD_var.SolverType == "FracStep"):
    fluid_only_model_part = ModelPart("FluidOnlyPart")

# Defining a model part for the solid part

balls_model_part = ModelPart("SolidPart")
solver_strategy.AddVariables(balls_model_part, fluid_model_part)

# Reading the solid part: binary or ascii, multifile or single

if (UD_var.OutputFileType == "Binary"):
    gid_mode = GiDPostMode.GiD_PostBinary

else:
    gid_mode = GiDPostMode.GiD_PostAscii

if (UD_var.Multifile == "multiple_files"):
    DEM_multifile = MultiFileFlag.MultipleFiles

else:
    DEM_multifile = MultiFileFlag.SingleFile

deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
write_conditions   = WriteConditionsFlag.WriteConditions
DEM_gid_io         = GidIO(UD_var.DEM_problem_name, gid_mode, DEM_multifile, deformed_mesh_flag, write_conditions)

# Introducing input file name

input_file_name = UD_var.problem_name

# Reading the fluid part

gid_mode             = GiDPostMode.GiD_PostBinary
multifile            = MultiFileFlag.MultipleFiles
deformed_mesh_flag   = WriteDeformedMeshFlag.WriteDeformed
write_conditions     = WriteConditionsFlag.WriteConditions
gid_io               = GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)
model_part_io_origin = ModelPartIO(input_file_name)

# Reading fluid model part

model_part_io_origin.ReadModelPart(fluid_model_part)

print "ULF model read correctly"

# Setting the fluid buffer size: SHOULD BE DONE AFTER READING!!!

fluid_model_part.SetBufferSize(3)

# Reading the DEM model part

model_part_io_solid = ModelPartIO(UD_var.DEM_problem_name)
model_part_io_solid.ReadModelPart(balls_model_part)

print "DEM model read correctly"

# Setting the DEM buffer size: SHOULD BE DONE AFTER READING!!!

balls_model_part.SetBufferSize(3)

# Adding dofs to the nodes of each model part

solver_strategy.AddDofs(balls_model_part, fluid_model_part)

if (UD_var.SolverType == "Quasi_Inc_Constant_Pressure" or UD_var.SolverType == "Quasi_Inc_Linear_Pressure"):

    for node in fluid_model_part.Nodes:
        node.Free(PRESSURE)

# Setting the limits of the bounding box

box_corner1 = Vector(3)
box_corner1[0] = UD_var.bounding_box_corner1_x
box_corner1[1] = UD_var.bounding_box_corner1_y
box_corner1[2] = UD_var.bounding_box_corner1_z

box_corner2 = Vector(3)
box_corner2[0] = UD_var.bounding_box_corner2_x
box_corner2[1] = UD_var.bounding_box_corner2_y
box_corner2[2] = UD_var.bounding_box_corner2_z

# Here we write the convergence data...

outstring2   = "convergence_info.txt"
outputfile1  = open(outstring2, 'w')

# Creation of solvers

# Fluid solver

fluid_solver  = solver_strategy.FluidStrategy(outputfile1, fluid_only_model_part, fluid_model_part, structure_model_part, combined_model_part, box_corner1, box_corner2)
fluid_solver.Initialize()

print "Fluid solver created and initialized"

# DEM_solver

DEM_solver = solver_strategy.DEMStrategy(balls_model_part, UD_var.DEM_domain_size) #here, solver variables initialize as default
ProcGiDSolverTransfer(balls_model_part, DEM_solver)
DEM_solver.Initialize()

print "DEM solver created and initialized"

# Checking to ensure that no node has zero density or pressure

is_fsi_interf = 0.0

[inverted_elements, domain_volume] = fluid_solver.CheckForInvertedElements()

# Calculating porosity

porosity_utils = PorosityUtils(domain_volume, balls_model_part)
porosity_utils.PrintCurrentData()
n_particles_in_depth = int(math.sqrt(porosity_utils.number_of_balls / domain_volume))

# Constitutive laws

if (UD_var.FSI):

    if (UD_var.domain_size == 2):
        fluid_model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic2D())

    elif (UD_var.domain_size == 3):
        fluid_model_part.Properties[1].SetValue(CONSTITUTIVE_LAW, Isotropic3D())

    else:
        raise "Domain size error. It should be 2D or 3D"

for node in fluid_model_part.Nodes:

    if (node.GetSolutionStepValue(DENSITY) == 0.0):
        print "node ",node.Id," has zero density!"
        raise 'node with zero density found'

    if (node.GetSolutionStepValue(VISCOSITY) == 0.0):
        print "node ",node.Id," has zero viscosity!"
        raise 'node with zero VISCOSITY found'

    if (UD_var.FSI):
        is_fsi_interf += node.GetSolutionStepValue(IS_INTERFACE)

if (UD_var.SolverType == "Incompressible_Modified_FracStep" and UD_var.FSI):

    if (is_fsi_interf == 0):
        raise 'For running FSI using the Modified Frac Step Solver you must prescribe IS_INTERFACE flag at the surface/outer contour of your structure'

# Adding dofs

solver_strategy.AddDofs(balls_model_part, fluid_model_part)

# Paths:

main_path             = os.getcwd()
post_path             = str(main_path) + '/' + str(UD_var.DEM_problem_name) + '_Post_Files'
list_path             = str(main_path) + '/' + str(UD_var.DEM_problem_name) + '_Post_Lists'
neigh_list_path       = str(main_path) + '/' + str(UD_var.DEM_problem_name) + '_Neigh_Lists'
data_and_results_path = str(main_path) + '/' + str(UD_var.DEM_problem_name) + '_Results_and_Data'
graphs_path           = str(main_path) + '/' + str(UD_var.DEM_problem_name) + '_Graphs'
MPI_results_path      = str(main_path) + '/' + str(UD_var.DEM_problem_name) + '_MPI_results'
fluid_results_path    = str(main_path) + '/' + str(UD_var.DEM_problem_name) + '_Fluid_results'
mixed_results_path    = str(main_path) + '/' + str(UD_var.DEM_problem_name) + '_Mixed_results'

for directory in [post_path, list_path, neigh_list_path, data_and_results_path, graphs_path, MPI_results_path, fluid_results_path, mixed_results_path]:

    if (not os.path.isdir(directory)):
        os.makedirs(str(directory))

os.chdir(data_and_results_path)
results           = open('results.txt', 'w')
summary_results   = open('summary_results.txt', 'w')

force_list        = []
force_list_2      = []
time_list         = []

os.chdir(list_path)
DEM_multifile     = open(UD_var.DEM_problem_name + '_all' + '.post.lst', 'w')
multifile_5       = open(UD_var.DEM_problem_name + '_5'   + '.post.lst', 'w')
multifile_10      = open(UD_var.DEM_problem_name + '_10'  + '.post.lst', 'w')
multifile_50      = open(UD_var.DEM_problem_name + '_50'  + '.post.lst', 'w')

DEM_multifile.write('Multiple\n')
multifile_5.write('Multiple\n')
multifile_10.write('Multiple\n')
multifile_50.write('Multiple\n')

first_print = True
index_5     = 1
index_10    = 1
index_50    = 1
prev_time   = 0.0
control     = 0.0

os.chdir(post_path)

if (UD_var.Multifile == "single_file"):
    DEM_gid_io.InitializeMesh(0.0)
    DEM_gid_io.WriteSphereMesh(balls_model_part.GetMesh())
    DEM_gid_io.FinalizeMesh()
    DEM_gid_io.InitializeResults(0.0, balls_model_part.GetMesh())

os.chdir(main_path)

# For plotting the graph:

velocity_node_y = 0.0

# Initializations

initial_dt         = 0.001 * UD_var.Dt
safety_factor      = 0.5
time               = 0.0
step               = 0
inlet_vel          = Vector(3)

DEM_dt             = balls_model_part.ProcessInfo.GetValue(DELTA_TIME)
DEM_time           = 0.0
DEM_step           = 0
DEM_time_old_print = 0.0

initial_pr_time    = timer.clock()
initial_real_time  = timer.time()

print ('\n' + 'Calculation starts at instant: ' + str(initial_pr_time) + '\n')
total_steps_expected = int(UD_var.final_time / DEM_dt)
print ('Total number of TIME STEPs expected in the calculation is: ' + str(total_steps_expected) + ' if time step is kept ' + '\n' )

if (UD_var.lagrangian_nodes_inlet == 1):
    for node in fluid_model_part.Nodes:
        if (node.GetSolutionStepValue(IS_LAGRANGIAN_INLET)):
            inlet_vel = node.GetSolutionStepValue(VELOCITY, 0)
            print "Lagrangian Inlet(s) Velocity  is ", inlet_vel
            break
else:
    inlet_vel[0] = 0.0
    inlet_vel[1] = 0.0
    inlet_vel[2] = 0.0

inlet_process = LagrangianInletProcess(fluid_model_part, 0.0, inlet_vel)


# Creation of projection module and initial projection

h_min = 0.01
projection_module = ProjectionModule(fluid_model_part, balls_model_part, UD_var.domain_size, n_particles_in_depth)
projection_module.UpdateDatabase(h_min)
projection_module.ProjectFromFluid()

# Temporal loop

while (time < UD_var.max_time):
    DEM_dt   = balls_model_part.ProcessInfo.GetValue(DELTA_TIME)
    DEM_time = DEM_time + DEM_dt

    if (DEM_time >= time):
        step += 1

        if (step <= 3):
            new_dt = 0.0000001
            time   = time + new_dt * safety_factor

        # Solving the fluid problem

        if (step > 3):
            new_dt = fluid_solver.EstimateDeltaTime(UD_var.Dt, UD_var.domain_size)
            time   = time + new_dt * safety_factor
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

    if (incremental_time > UD_var.control_time):
        percentage = 100.0 * (float(DEM_step) / total_steps_expected)
        print 'Real time calculation: ' + str(timer.time() - initial_real_time)
        print 'Percentage Completed: '  + str(percentage) + ' %'
        print "TIME STEP = "            + str(DEM_step) + '\n'
        prev_time = (timer.time() - initial_real_time)

    if ((timer.time() - initial_real_time > 60.0) and first_print == True):
        first_print = False
        estimation_time = 60.0 * (total_steps_expected / DEM_step) #seconds
        print('The total calculation estimated time is ' + str(estimation_time) + 'seconds.' + '\n')
        print('In minutes :' + str(estimation_time / 60) + 'min.' + '\n')
        print('In hours :' + str(estimation_time / 3600) + 'hrs.' + '\n')
        print('In days :' + str(estimation_time / 86400) + 'days.' + '\n')

        if (estimation_time / 86400 > 2.0):
            print('WARNING!!!:       VERY LASTING CALCULATION'+'\n')

    os.chdir(list_path)
    DEM_multifile.write(UD_var.DEM_problem_name + '_' + str(DEM_time) + '.post.bin\n')
    os.chdir(main_path)

    DEM_step += 1

  ##############     GiD IO        ################################################################################

    time_to_print = DEM_time - DEM_time_old_print
    if (time_to_print >= UD_var.output_dt):
        os.chdir(main_path)

        if (UD_var.PrintNeighbourLists == "ON"): # Printing neighbours id's
            os.chdir(neigh_list_path)
            neighbours_list = open('neigh_list_' + str(DEM_time), 'w')

            for elem in balls_model_part.Elements:
                ID = (elem.Id)
                Neigh_ID = elem.GetValue(NEIGHBOURS_IDS)

                for i in range(len(Neigh_ID)):
                    neighbours_list.write(str(ID) + ' ' + str(Neigh_ID[i]) + '\n')

            neighbours_list.close()

        if (UD_var.Multifile == "multiple_files"):
            os.chdir(mixed_results_path)
            ParticleUtils2D().VisualizationModelPart(mixed_model_part, fluid_model_part, balls_model_part)
            DEM_gid_io.InitializeMesh(time)
            DEM_gid_io.WriteSphereMesh(balls_model_part.GetMesh())
            DEM_gid_io.WriteMesh(mixed_model_part.GetMesh())
            DEM_gid_io.FinalizeMesh()
            DEM_gid_io.InitializeResults(10.0, mixed_model_part.GetMesh())
            DEM_gid_io.WriteNodalResults(VELOCITY, mixed_model_part.Nodes,time, 0)
            DEM_gid_io.WriteNodalResults(DISPLACEMENT, mixed_model_part.Nodes,time, 0)
            DEM_gid_io.WriteNodalResults(PRESSURE, fluid_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(PRESSURE_GRADIENT, fluid_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(DRAG_REACTION, fluid_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(RADIUS, balls_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(FLUID_VEL_PROJECTED, balls_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(DRAG_FORCE, balls_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(BUOYANCY, balls_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(PRESSURE_GRAD_PROJECTED, balls_model_part.Nodes, time, 0)
            DEM_gid_io.WriteNodalResults(TOTAL_FORCES, balls_model_part.Nodes, time, 0)

            if (print_export_id == "1"):
                DEM_gid_io.WriteNodalResults(EXPORT_ID, balls_model_part.Nodes, time, 0)

            DEM_gid_io.Flush()
            DEM_gid_io.FinalizeResults()

        ##########PRINTING VARIABLES############

        os.chdir(data_and_results_path)

        if (index_5 == 5):
            multifile_5.write(UD_var.DEM_problem_name + '_' + str(DEM_time) + '.post.bin\n')
            index_5 = 0

        if (index_10 == 10):
            multifile_10.write(UD_var.DEM_problem_name + '_' + str(DEM_time) +  '.post.bin\n')
            index_10 = 0

        if (index_50 == 50):
            multifile_50.write(UD_var.DEM_problem_name + '_' + str(DEM_time) + '.post.bin\n')
            index_50 = 0

        index_5 += 1
        index_10 += 1
        index_50 += 1

        if (UD_var.Multifile == "multiple_files"):
            DEM_gid_io.FinalizeResults()

        os.chdir(main_path)
        DEM_time_old_print = DEM_time

    #End of print loop

    os.chdir(main_path)

# End of temporal loop

if (UD_var.Multifile == "single_file"):
    DEM_gid_io.FinalizeResults()

os.chdir(data_and_results_path)
results.close()
summary_results.close()

os.chdir(list_path)
DEM_multifile.close()
multifile_5.close()
multifile_10.close()
multifile_50.close()

os.chdir(main_path)
elapsed_pr_time = timer.clock() - initial_pr_time
elapsed_real_time = timer.time() - initial_real_time
print 'Calculation ends at instant: '                 + str(timer.time())
print 'Calculation ends at processing time instant: ' + str(timer.clock())
print 'Elapsed processing time: '                     + str(elapsed_pr_time)
print 'Elapsed real time: '                           + str(elapsed_real_time)
print (my_timer)
print "ANALYSIS COMPLETED"
