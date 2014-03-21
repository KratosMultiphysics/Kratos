from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import the configuration data as read from the GiD
import ProjectParameters
import define_output

# setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
import math
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
import DEM_explicit_solver_var as DEMParameters
import DEM_procedures
import swimming_DEM_procedures
import embedded
import mesh_motion

# PROJECT PARAMETERS (to be put in problem type)
ProjectParameters.projection_module_option         = 1
ProjectParameters.print_particles_results_option   = 0
ProjectParameters.project_at_every_substep_option  = 0
ProjectParameters.velocity_trap_option             = 0
ProjectParameters.inlet_option                     = 1
ProjectParameters.manually_imposed_drag_law_option = 0
ProjectParameters.stationary_problem_option        = 0 # inactive (0), stop calculating the fluid after it reaches the stationary state (1)
ProjectParameters.body_force_on_fluid_option       = 0
ProjectParameters.similarity_transformation_type   = 0 # no transformation (0), Tsuji (1)
ProjectParameters.dem_inlet_element_type           = "SphericSwimmingParticle3D"  # "SphericParticle3D", "SphericSwimmingParticle3D"
ProjectParameters.coupling_scheme_type             = "UpdatedFluid" # "UpdatedFluid", "UpdatedDEM"
ProjectParameters.coupling_weighing_type           = 2 # {fluid_to_DEM, DEM_to_fluid, Solid_fraction} = {lin, const, const} (0), {lin, lin, const} (1), {lin, lin, lin} (2)
ProjectParameters.buoyancy_force_type              = 1 # null buoyancy (0), standard (1) (but if drag_force_type is 2, buoyancy is always parallel to gravity) 
ProjectParameters.drag_force_type                  = 2 # null drag (0), standard (1), Weatherford (2), Ganser (3)
ProjectParameters.virtual_mass_force_type          = 0 # null virtual mass force (0)
ProjectParameters.lift_force_type                  = DEMParameters.consider_lift_force_option # Yes (1) (Weatherford recommendation), No (0)
ProjectParameters.drag_modifier_type               = DEMParameters.drag_modifier_type # Hayder (2), Chien (3)  # problemtype option
ProjectParameters.interaction_start_time           = 0.00
ProjectParameters.max_solid_fraction               = 0.6
ProjectParameters.initial_drag_force               = 0.0   # problemtype option
ProjectParameters.drag_law_slope                   = 0.0   # problemtype option
ProjectParameters.power_law_tol                    = 0.0
ProjectParameters.model_over_real_diameter_factor  = 1.0 # not active if similarity_transformation_type = 0
ProjectParameters.max_pressure_variation_rate_tol  = 1e-4 # for stationary problems, criterion to stop the fluid calculations
ProjectParameters.time_steps_per_stationarity_step = 15 # number of fluid time steps between consecutive assessment of stationarity steps

# variables to be printed
ProjectParameters.dem_nodal_results                = ["RADIUS", "FLUID_VEL_PROJECTED", "DRAG_FORCE", "LIFT_FORCE", "BUOYANCY", "PRESSURE_GRAD_PROJECTED"]
ProjectParameters.mixed_nodal_results              = ["VELOCITY", "DISPLACEMENT"]
ProjectParameters.variables_to_print_in_file       = ["DRAG_FORCE", "LIFT_FORCE", "BUOYANCY", "VELOCITY"]

# changes on PROJECT PARAMETERS for the sake of consistency
ProjectParameters.nodal_results.append("SOLID_FRACTION")
ProjectParameters.nodal_results.append("MESH_VELOCITY1")
ProjectParameters.nodal_results.append("BODY_FORCE")
ProjectParameters.nodal_results.append("DRAG_REACTION")
ProjectParameters.nodal_results.append("DISTANCE") # embedded

DEMParameters.project_from_particles_option *= ProjectParameters.projection_module_option
ProjectParameters.project_at_every_substep_option *= ProjectParameters.projection_module_option
ProjectParameters.time_steps_per_stationarity_step = max(1, int(ProjectParameters.time_steps_per_stationarity_step))  # it should never be smaller than 1!
ProjectParameters.stationary_problem_option *= not DEMParameters.project_from_particles_option

for var in ProjectParameters.mixed_nodal_results:

    if var in ProjectParameters.nodal_results:
        ProjectParameters.nodal_results.remove(var)

if (ProjectParameters.body_force_on_fluid_option==0):
    ProjectParameters.nodal_results.remove("PRESSURE")
# choosing the variables to be passed as a parameter to the constructor of the
# ProjectionModule class constructor that are to be filled with the other phase's
# information through the coupling process

fluid_vars_for_coupling = [DRAG_REACTION,
                           BODY_FORCE,
                           SOLID_FRACTION,
                           MESH_VELOCITY1]

balls_vars_for_coupling = [FLUID_VEL_PROJECTED,
                           FLUID_DENSITY_PROJECTED,
                           PRESSURE_GRAD_PROJECTED,
                           FLUID_VISCOSITY_PROJECTED,
                           SOLID_FRACTION_PROJECTED,
                           SHEAR_RATE_PROJECTED,
                           FLUID_VORTICITY_PROJECTED,
                           POWER_LAW_N,
                           POWER_LAW_K,
                           YIELD_STRESS,
                           GEL_STRENGTH,
                           DISTANCE] # <- REQUIRED BY EMBEDDED

# extra nodal variables to be added to the model parts (memory will be allocated for them)
fluid_variables_to_add = []
fluid_variables_to_add += fluid_vars_for_coupling
fluid_variables_to_add += [PRESSURE_GRADIENT,
                           SOLID_FRACTION_GRADIENT,
                           AUX_DOUBLE_VAR,
                           YIELD_STRESS,
                           BINGHAM_SMOOTHER,
                           POWER_LAW_N,
                           POWER_LAW_K,
                           GEL_STRENGTH,
                           DISTANCE]  # <- REQUIRED BY EMBEDDED

balls_variables_to_add = []
balls_variables_to_add += balls_vars_for_coupling
balls_variables_to_add += [SPHERICITY,
                           DRAG_FORCE,
                           LIFT_FORCE,
                           BUOYANCY,
                           REYNOLDS_NUMBER]

fem_dem_variables_to_add = [VELOCITY,
                            DISPLACEMENT,
                            FORCE,  # <- REQUIRED BY EMBEDDED
                            POSITIVE_FACE_PRESSURE,  # <- REQUIRED BY EMBEDDED
                            NEGATIVE_FACE_PRESSURE]  # <- REQUIRED BY EMBEDDED

DEM_inlet_variables_to_add = balls_variables_to_add

# constructing a DEM_procedures object
dem_procedures = DEM_procedures.Procedures(DEMParameters)
# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

# defining variables to be used
# GID IO IS NOT USING THIS NOW. TO BE REMOVED ONCE THE "PRINT IN POINTS"
# CODE IS NOT USING IT

variables_dictionary = {"PRESSURE": PRESSURE,
                        "VELOCITY": VELOCITY,
                        "MU": MU,  # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
                        "BUOYANCY": BUOYANCY,  # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
                        "DRAG_FORCE": DRAG_FORCE,  # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
                        "LIFT_FORCE": LIFT_FORCE}  # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#                        "REACTION": REACTION,
#                        "DISTANCE": DISTANCE, }

# defining a model part for the fluid part
fluid_model_part = ModelPart("FluidPart")

if "REACTION" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(REACTION)
if "DISTANCE" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

#
#
# importing the solvers needed
SolverSettings = ProjectParameters.FluidSolverConfiguration
solver_module = import_solver(SolverSettings)

#
# importing variables
print('Adding nodal variables to the fluid_model_part')  # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
sys.stdout.flush()

# caution with breaking up this block! # (memory allocation) {
solver_module.AddVariables(fluid_model_part, SolverSettings)
swimming_DEM_procedures.AddNodalVariables(fluid_model_part, fluid_variables_to_add)  # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
# }

# introducing input file name
input_file_name = ProjectParameters.problem_name

# reading the fluid part
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS


for node in fluid_model_part.Nodes:
    first_node_yield_stress = node.GetSolutionStepValue(YIELD_STRESS)
    break

if(first_node_yield_stress == 0.0):
    non_newtonian_option = 0

# defining model parts for the balls part and for the DEM-FEM interaction elements

balls_model_part = ModelPart("SolidPart")
fem_dem_model_part = ModelPart("RigidFace_Part");

import sphere_strategy as SolverStrategy
SolverStrategy.AddVariables(balls_model_part, DEMParameters)

print('Adding nodal variables to the balls_model_part, necessary for the DEM-Fluid interaction...')  # (memory allocation)
sys.stdout.flush()

swimming_DEM_procedures.AddNodalVariables(balls_model_part, balls_variables_to_add)

print('Adding nodal variables to the dem_fem_wall_model_part...')  # (memory allocation)
sys.stdout.flush()

swimming_DEM_procedures.AddNodalVariables(fem_dem_model_part, fem_dem_variables_to_add)

# defining a model part for the mixed part
mixed_model_part = ModelPart("MixedPart")

# reading the balls model part
model_part_io_solid = ModelPartIO(DEMParameters.problem_name + "DEM")
model_part_io_solid.ReadModelPart(balls_model_part)
rigidFace_mp_filename = DEMParameters.problem_name + "DEM_FEM_boundary"
model_part_io_solid = ModelPartIO(rigidFace_mp_filename)
model_part_io_solid.ReadModelPart(fem_dem_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
balls_model_part.SetBufferSize(1)

# adding nodal degrees of freedom
SolverStrategy.AddDofs(balls_model_part)

# adding extra process info variables
balls_model_part.ProcessInfo.SetValue(BUOYANCY_FORCE_TYPE, ProjectParameters.buoyancy_force_type)
balls_model_part.ProcessInfo.SetValue(DRAG_FORCE_TYPE, ProjectParameters.drag_force_type)
balls_model_part.ProcessInfo.SetValue(LIFT_FORCE_TYPE, ProjectParameters.lift_force_type)
balls_model_part.ProcessInfo.SetValue(VIRTUAL_MASS_FORCE_TYPE, ProjectParameters.virtual_mass_force_type)
# balls_model_part.ProcessInfo.SetValue(NON_NEWTONIAN_OPTION, non_newtonian_option)
balls_model_part.ProcessInfo.SetValue(MANUALLY_IMPOSED_DRAG_LAW_OPTION, ProjectParameters.manually_imposed_drag_law_option)
balls_model_part.ProcessInfo.SetValue(DRAG_MODIFIER_TYPE, ProjectParameters.drag_modifier_type)
balls_model_part.ProcessInfo.SetValue(INIT_DRAG_FORCE, ProjectParameters.initial_drag_force)
balls_model_part.ProcessInfo.SetValue(DRAG_LAW_SLOPE, ProjectParameters.drag_law_slope)
balls_model_part.ProcessInfo.SetValue(POWER_LAW_TOLERANCE, ProjectParameters.power_law_tol)

# creating a distance calculation process for the embedded technology
# (used to calculate elemental distances defining the structure embedded in the fluid mesh)
calculate_distance_process = CalculateSignedDistanceTo3DSkinProcess(fem_dem_model_part, fluid_model_part)
calculate_distance_process.Execute()

# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA


# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

# adding nodal degrees of freedom
solver_module.AddDofs(fluid_model_part, SolverSettings)

# If Lalplacian form = 2, free all pressure Dofs
# laplacian_form = ProjectParameters.laplacian_form
# if(laplacian_form >= 2):
    # for node in fluid_model_part.Nodes:
        # node.Free(PRESSURE)

# copy Y_WALL
for node in fluid_model_part.Nodes:
    y = node.GetSolutionStepValue(Y_WALL, 0)
    node.SetValue(Y_WALL, y)

#
#
# Creating the fluid solver
fluid_solver = solver_module.CreateSolver(
    fluid_model_part, SolverSettings)

# activate turbulence model
if(ProjectParameters.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras"):
    # apply the initial turbulent viscosity on all of the nodes
    turb_visc = SolverSettings.TurbulentViscosity
    for node in fluid_model_part.Nodes:
        node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, turb_visc)
        visc = node.GetSolutionStepValue(VISCOSITY)
        node.SetSolutionStepValue(MOLECULAR_VISCOSITY, 0, visc)
        if (node.IsFixed(VELOCITY_X) and node.GetSolutionStepValue(VELOCITY_X, 0) != 0.0) or \
           (node.IsFixed(VELOCITY_Y) and node.GetSolutionStepValue(VELOCITY_Y, 0) != 0.0) or \
           (node.IsFixed(VELOCITY_Z) and node.GetSolutionStepValue(VELOCITY_Z, 0) != 0.0):
            node.Fix(TURBULENT_VISCOSITY)

    # select nodes on the wall
    fluid_solver.wall_nodes = []
    for i in SolverSettings.SA_wall_group_ids:
        # get the nodes of the wall for SA.
        nodes = fluid_model_part.GetNodes(i)
        for node in nodes:
            fluid_solver.wall_nodes.append(node)
            node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, 0.0)
            node.Fix(TURBULENT_VISCOSITY)


fluid_solver.Initialize()
print("fluid solver created")
sys.stdout.flush()

# initialize GiD  I/O
from gid_output import GiDOutput
gid_io = GiDOutput(input_file_name,
                   ProjectParameters.VolumeOutput,
                   ProjectParameters.GiDPostMode,
                   ProjectParameters.GiDMultiFileFlag,
                   ProjectParameters.GiDWriteMeshFlag,
                   ProjectParameters.GiDWriteConditionsFlag)

if not ProjectParameters.VolumeOutput:
    cut_list = define_output.DefineCutPlanes()
    gid_io.define_cuts(fluid_model_part, cut_list)

# gid_io.initialize_results(fluid_model_part) # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
import swimming_DEM_gid_output
swimming_DEM_gid_io = swimming_DEM_gid_output.SwimmingDEMGiDOutput(
                                                                   input_file_name,
                                                                   ProjectParameters.VolumeOutput,
                                                                   ProjectParameters.GiDPostMode,
                                                                   ProjectParameters.GiDMultiFileFlag,
                                                                   ProjectParameters.GiDWriteMeshFlag,
                                                                   ProjectParameters.GiDWriteConditionsFlag)

swimming_DEM_gid_io.initialize_swimming_DEM_results(balls_model_part, fem_dem_model_part, mixed_model_part)
# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

# define the drag computation list
drag_list = define_output.DefineDragList()
drag_file_output_list = []
for it in drag_list:
    f = open(it[1], 'w')
    drag_file_output_list.append(f)
    tmp = "#Drag for group " + it[1] + "\n"
    f.write(tmp)
    tmp = "time RX RY RZ"
    f.write(tmp)
    f.flush()

print(drag_file_output_list)


def PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time):
    i = 0
    for it in drag_list:
        print(it[0])
        nodes = fluid_model_part.GetNodes(it[0])
        drag = Vector(3)
        drag[0] = 0.0
        drag[1] = 0.0
        drag[2] = 0.0
        for node in nodes:
            reaction = node.GetSolutionStepValue(REACTION, 0)
            drag[0] += reaction[0]
            drag[1] += reaction[1]
            drag[2] += reaction[2]

        output = str(time) + " " + str(drag[0]) + " " + str(
            drag[1]) + " " + str(drag[2]) + "\n"
        # print drag_file_output_list[i]
        # print output
        drag_file_output_list[i].write(output)
        drag_file_output_list[i].flush()
        i = i + 1


# preparing output of point graphs
import point_graph_printer

output_nodes_list = define_output.DefineOutputPoints()
graph_printer = point_graph_printer.PrintGraphPrinter(
    output_nodes_list,
    fluid_model_part,
    variables_dictionary,
    domain_size)

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
# applying changes to the physiscal properties of the model to adjust for
# the similarity transformation if required (fluid effects only).
swimming_DEM_procedures.ApplySimilarityTransformations(fluid_model_part, ProjectParameters.similarity_transformation_type, ProjectParameters.model_over_real_diameter_factor)

max_fluid_node_Id = swimming_DEM_procedures.FindMaxNodeIdInFLuid(fluid_model_part)

# creating a Post Utils object that executes several post-related tasks
post_utils = DEM_procedures.PostUtils(DEMParameters, balls_model_part)
post_utilities = post_utils.post_utilities

# creating an IOTools object to perform other printing tasks
io_tools = swimming_DEM_procedures.IOTools(ProjectParameters)

# creating a projection module for the fluid-DEM coupling
h_min = 0.01
n_balls = 1
fluid_volume = 10
n_particles_in_depth = int(math.sqrt(n_balls / fluid_volume))

if (ProjectParameters.projection_module_option):
    projection_module = swimming_DEM_procedures.ProjectionModule(fluid_model_part, balls_model_part, fem_dem_model_part, domain_size, ProjectParameters.max_solid_fraction, ProjectParameters.coupling_weighing_type, n_particles_in_depth, balls_vars_for_coupling, fluid_vars_for_coupling)
    #to do: the projection module should be created with a list of variables to be projected (one list for every direction)
    projection_module.UpdateDatabase(h_min)
    interaction_calculator = CustomFunctionsCalculator()

if (ProjectParameters.body_force_on_fluid_option):

    for node in fluid_model_part.Nodes:
        node.SetSolutionStepValue(BODY_FORCE_X, 0, DEMParameters.GravityX)
        node.SetSolutionStepValue(BODY_FORCE_Y, 0, DEMParameters.GravityY)
        node.SetSolutionStepValue(BODY_FORCE_Z, 0, DEMParameters.GravityZ)

# creating a CreatorDestructor object, encharged of any adding or removing of elements during the simulation
creator_destructor = ParticleCreatorDestructor()
creator_destructor.SetMaxNodeId(max_fluid_node_Id)

# creating a Solver object for the DEM part. It contains the sequence of function calls necessary for the evolution of the DEM system at every time step
dem_solver = SolverStrategy.ExplicitStrategy(balls_model_part, fem_dem_model_part, creator_destructor, DEMParameters)

# constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation

if (ProjectParameters.inlet_option):
    DEM_inlet_model_part = ModelPart("DEMInletPart")
    DEM_Inlet_filename = DEMParameters.problem_name + "DEM_Inlet"
    SolverStrategy.AddVariables(DEM_inlet_model_part, DEMParameters)
    swimming_DEM_procedures.AddNodalVariables(DEM_inlet_model_part, DEM_inlet_variables_to_add)
    model_part_io_demInlet = ModelPartIO(DEM_Inlet_filename)
    model_part_io_demInlet.ReadModelPart(DEM_inlet_model_part)

    # setting up the buffer size:
    DEM_inlet_model_part.SetBufferSize(1)

    # adding nodal degrees of freedom
    SolverStrategy.AddDofs(DEM_inlet_model_part)
    DEM_inlet_parameters = DEM_inlet_model_part.Properties

    DEM_inlet_model_part.ProcessInfo.SetValue(BUOYANCY_FORCE_TYPE, ProjectParameters.buoyancy_force_type)
    DEM_inlet_model_part.ProcessInfo.SetValue(DRAG_FORCE_TYPE, ProjectParameters.drag_force_type)
    DEM_inlet_model_part.ProcessInfo.SetValue(LIFT_FORCE_TYPE, ProjectParameters.lift_force_type)
    DEM_inlet_model_part.ProcessInfo.SetValue(VIRTUAL_MASS_FORCE_TYPE, ProjectParameters.virtual_mass_force_type)
    DEM_inlet_model_part.ProcessInfo.SetValue(MANUALLY_IMPOSED_DRAG_LAW_OPTION, ProjectParameters.manually_imposed_drag_law_option)
    DEM_inlet_model_part.ProcessInfo.SetValue(DRAG_MODIFIER_TYPE, ProjectParameters.drag_modifier_type)
    DEM_inlet_model_part.ProcessInfo.SetValue(INIT_DRAG_FORCE, ProjectParameters.initial_drag_force)
    DEM_inlet_model_part.ProcessInfo.SetValue(DRAG_LAW_SLOPE, ProjectParameters.drag_law_slope)
    DEM_inlet_model_part.ProcessInfo.SetValue(POWER_LAW_TOLERANCE, ProjectParameters.power_law_tol)

    # constructiong the inlet and intializing it
    DEM_inlet = DEM_Inlet(DEM_inlet_model_part)
    DEM_inlet.InitializeDEM_Inlet(balls_model_part, creator_destructor)

# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

# renumerating IDs if required

# swimming_DEM_procedures.RenumberNodesIdsToAvoidRepeating(fluid_model_part, balls_model_part, fem_dem_model_part)

# Stepping and time settings

Dt = ProjectParameters.Dt
Nsteps = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = Dt
step = 0

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
time_dem = 0.0
DEM_step = 0  # this variable is necessary to get a good random insertion of particles
stat_steps = 0
Dt_DEM = DEMParameters.MaxTimeStep
stationarity = False
dem_solver.Initialize()
# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

while(time <= final_time):

    time = time + Dt
    step += 1
    stat_steps += 1
    fluid_model_part.CloneTimeStep(time)
    print("\n", "TIME = ", time)
    sys.stdout.flush()

    if (ProjectParameters.coupling_scheme_type == "UpdatedDEM"):
        time_final_DEM_substepping = time + Dt
        time_dem = time + Dt_DEM

    else:
        time_final_DEM_substepping = time
        time_dem = time - Dt + Dt_DEM

    # walls movement    
    mesh_motion.MoveAllMeshes(fem_dem_model_part, time)
    
    # Calculate elemental distances defining the structure embedded in the fluid mesh
    calculate_distance_process.Execute()

    if(step >= 3):  # MA: because I think DISTANCE,1 (from previous time step) is not calculated correctly for step=1
        embedded.ApplyEmbeddedBCsToFluid(fluid_model_part)
        embedded.ApplyEmbeddedBCsToBalls(balls_model_part,DEMParameters)

    # solving the fluid part

    if (step >= 3 and not stationarity):
        print("Solving Fluid... (", fluid_model_part.NumberOfElements(0), " elements )")
        sys.stdout.flush()
        fluid_solver.Solve()               

    # assessing stationarity

        if (stat_steps >= ProjectParameters.time_steps_per_stationarity_step and ProjectParameters.stationary_problem_option):
            print("Assessing Stationarity...")
            sys.stdout.flush()
            stat_steps = 0
            stationarity = interaction_calculator.AssessStationarity(fluid_model_part, ProjectParameters.max_pressure_variation_rate_tol)  # in the first time step the 'old' pressure vector is created and filled

            if (stationarity):
                print("**************************************************************************************************")
                print()
                print("The model has reached a stationary state. The fluid calculation is suspended.")
                print()
                print("**************************************************************************************************")
                sys.stdout.flush()

    # printing if required

    if (ProjectParameters.print_particles_results_option):
        io_tools.PrintParticlesResults(ProjectParameters.variables_to_print_in_file, time, balls_model_part)
        graph_printer.PrintGraphs(time)
        PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

    if (output_time <= out and ProjectParameters.coupling_scheme_type == "UpdatedDEM"):

        if (ProjectParameters.projection_module_option):
            projection_module.ComputePostProcessResults(balls_model_part.ProcessInfo)

        if (ProjectParameters.GiDMultiFileFlag == "Multiples"):
            mixed_model_part.Elements.clear()
            mixed_model_part.Nodes.clear()
            post_utilities.AddModelPartToModelPart(mixed_model_part, fluid_model_part)
            post_utilities.AddModelPartToModelPart(mixed_model_part, balls_model_part)
            post_utilities.AddModelPartToModelPart(mixed_model_part, fem_dem_model_part)

        swimming_DEM_gid_io.write_swimming_DEM_results(time, fluid_model_part, balls_model_part, fem_dem_model_part, mixed_model_part, ProjectParameters.nodal_results, ProjectParameters.dem_nodal_results, ProjectParameters.mixed_nodal_results, ProjectParameters.gauss_points_results)
        out = 0

    # solving the DEM part

    if (time >= ProjectParameters.interaction_start_time and ProjectParameters.projection_module_option):
        interaction_calculator.CalculatePressureGradient(fluid_model_part)

    print("Solving DEM... (", balls_model_part.NumberOfElements(0), " elements)")
    sys.stdout.flush()
    first_dem_iter = True
    Dt_DEM_tolerance = 1e-9

    while time_dem < (time_final_DEM_substepping + Dt_DEM):

        if time_dem > time_final_DEM_substepping - Dt_DEM_tolerance:
            time_dem = time_final_DEM_substepping

        DEM_step += 1   # this variable is necessary to get a good random insertion of particles

        balls_model_part.ProcessInfo[TIME_STEPS] = DEM_step

        # applying fluid-to-DEM coupling if required

        if (time >= ProjectParameters.interaction_start_time and ProjectParameters.projection_module_option and (ProjectParameters.project_at_every_substep_option or first_dem_iter)):

            if (ProjectParameters.coupling_scheme_type == "UpdatedDEM"):
                projection_module.ProjectFromNewestFluid()

            else:
                projection_module.ProjectFromFluid((time_final_DEM_substepping - time_dem) / Dt)

        # performing the time integration of the DEM part

        balls_model_part.CloneTimeStep(time_dem)

        # actual_Dt_DEM = balls_model_part.ProcessInfo[DELTA_TIME] #uncomment if you want to use the actual Dt for the DEM (obtained by CloneTimeStep)
        # print "actual_Dt_DEM = ", actual_Dt_DEM

        dem_solver.Solve()

        # adding DEM elements by the inlet

        if (ProjectParameters.inlet_option):
            DEM_inlet.CreateElementsFromInletMesh(balls_model_part, DEM_inlet_model_part, creator_destructor, ProjectParameters.dem_inlet_element_type)  # After solving, to make sure that neighbours are already set.

        first_dem_iter = False

        time_dem += Dt_DEM

    # measuring mean velocities in a certain control volume (the 'velocity trap')

    if (ProjectParameters.velocity_trap_option):
        post_utils.ComputeMeanVelocitiesinTrap("Average_Velocity", time_dem)

    # applying DEM-to-fluid coupling

    if (time >= ProjectParameters.interaction_start_time and DEMParameters.project_from_particles_option):
        print("Projecting from particles to the fluid...")
        sys.stdout.flush()
        projection_module.ProjectFromParticles()

    # printing if required

    if (ProjectParameters.print_particles_results_option):
        io_tools.PrintParticlesResults(ProjectParameters.variables_to_print_in_file, time, balls_model_part)
        graph_printer.PrintGraphs(time)
        PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

    if (output_time <= out and ProjectParameters.coupling_scheme_type == "UpdatedFluid"):

        if (ProjectParameters.projection_module_option):
            projection_module.ComputePostProcessResults(balls_model_part.ProcessInfo)

        print("")
        print("*******************  PRINTING RESULTS FOR GID  ***************************")
        sys.stdout.flush()

        if (ProjectParameters.GiDMultiFileFlag == "Multiples"):
            mixed_model_part.Elements.clear()
            mixed_model_part.Nodes.clear()
            post_utilities.AddModelPartToModelPart(mixed_model_part, fem_dem_model_part) # order can be important here!! --> adding the fluid after the fem_dem ...
            post_utilities.AddModelPartToModelPart(mixed_model_part, fluid_model_part) #  ...overwrites the results on the repeated nodes !!! This is useful for IS_SLIP conditions
            post_utilities.AddModelPartToModelPart(mixed_model_part, balls_model_part)
            

        swimming_DEM_gid_io.write_swimming_DEM_results(time, fluid_model_part, balls_model_part, fem_dem_model_part, mixed_model_part, ProjectParameters.nodal_results, ProjectParameters.dem_nodal_results, ProjectParameters.mixed_nodal_results, ProjectParameters.gauss_points_results)
        out = 0

        
    out = out + Dt

swimming_DEM_gid_io.finalize_results()

print("CALCULATIONS FINISHED. THE SIMULATION ENDED SUCCESSFULLY.")
sys.stdout.flush()

for i in drag_file_output_list:
    i.close()
