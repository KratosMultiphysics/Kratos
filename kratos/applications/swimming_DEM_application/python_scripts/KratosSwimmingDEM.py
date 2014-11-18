# This script contains an algorithm that models fluid-particle interaction.
# It combines two parts: a FEM model for the fluid and a DEM model for the particles.
# It has been conceived by adding the DEM part and the interaction on top of an original fluid-only script (see kratos/applications/FluidDynamicswApplication)
# Some parts of the original fluid script have been kept practically untouched and are clearly marked.
# Whenever a minor modification has been made on one of these parts, the corresponding line is indicated with a comment: # MOD.

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    B E G I N S
#_____________________________________________________________________________________________________________________________________

from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import the configuration data as read from the GiD
import ProjectParameters as pp # MOD
import define_output

# setting the domain size for the problem to be solved
domain_size = pp.domain_size

import sys
sys.path.append(pp.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    E N D S
#_____________________________________________________________________________________________________________________________________

import math
import os
current_script_path = os.path.dirname(os.path.abspath(__file__))

from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
import DEM_explicit_solver_var
import DEM_procedures
import swimming_DEM_procedures as swim_proc
import CFD_DEM_coupling
import variables_management as vars_man
import embedded

# listing project parameters (to be put in problem type)
pp.dem                                    = DEM_explicit_solver_var
pp.projection_module_option               = 1
pp.print_particles_results_option         = 0
pp.add_each_hydro_force_option            = 1 # add each of the hydrodynamic forces (drag, lift and virtual mass)
pp.project_at_every_substep_option        = 1
pp.velocity_trap_option                   = 0
pp.inlet_option                           = 1
pp.manually_imposed_drag_law_option       = 0
pp.stationary_problem_option              = 0 # stationary, stop calculating the fluid after it reaches the stationary state
pp.flow_in_porous_medium_option           = 0 # the porosity is an imposed field
pp.flow_in_porous_DEM_medium_option       = 0 # the DEM part is kept static
pp.embedded_option                        = 1 # the embedded domain tools are to be used
pp.make_results_directories_option        = 1 # results are written into a folder (../results) inside the problem folder
pp.body_force_on_fluid_option             = 0
pp.print_debug_info_option                = 0 # print a summary of global physical measures
pp.print_REYNOLDS_NUMBER_option           = 1
pp.print_PRESSURE_GRAD_PROJECTED_option   = 1
pp.print_FLUID_VEL_PROJECTED_option       = 1
pp.print_FLUID_ACCEL_PROJECTED_option     = 1
pp.print_BUOYANCY_option                  = 1
pp.print_DRAG_FORCE_option                = 1
pp.print_VIRTUAL_MASS_FORCE_option        = 1
pp.print_LIFT_FORCE_option                = 1
pp.print_SOLID_FRACTION_option            = 1
pp.print_FLUID_FRACTION_option            = 1
pp.print_FLUID_VISCOSITY_PROJECTED_option = 1
pp.print_FLUID_FRACTION_PROJECTED_option  = 1
pp.print_MESH_VELOCITY1_option            = 1
pp.print_FLUID_FRACTION_GRADIENT_option   = 1
pp.print_BODY_FORCE_option                = 1
pp.print_HYDRODYNAMIC_REACTION_option     = 1
pp.print_HYDRODYNAMIC_FORCE_option        = 1
pp.print_HYDRODYNAMIC_MOMENT_option       = 1
pp.print_PRESSURE_option                  = 1
pp.print_particles_results_cycle          = 1 # number of 'ticks' per printing cycle
pp.debug_tool_cycle                       = 10 # number of 'ticks' per debug computations cycle
pp.similarity_transformation_type         = 0 # no transformation (0), Tsuji (1)
pp.dem_inlet_element_type                 = "SphericSwimmingParticle3D"  # "SphericParticle3D", "SphericSwimmingParticle3D"
pp.fluid_model_type                       = 0 # untouched, velocity incremented by 1/fluid_fraction (0), modified mass conservation only (1)
pp.coupling_level_type                    = pp.dem.project_from_particles_option # one way coupling (0), two way coupling (1)
pp.coupling_scheme_type                   = "UpdatedFluid" # "UpdatedFluid", "UpdatedDEM"
pp.coupling_weighing_type                 = 2 # {fluid_to_DEM, DEM_to_fluid, fluid_fraction} = {lin, lin, imposed} (-1), {lin, const, const} (0), {lin, lin, const} (1), {lin, lin, lin} (2), averaging method (3)
pp.buoyancy_force_type                    = 1 # null buoyancy (0), compute buoyancy (1)  if drag_force_type is 2 buoyancy is always parallel to gravity
pp.drag_force_type                        = 2 # null drag (0), Stokes (1), Weatherford (2), Ganser (3), Ishii (4), Newtonian Regime (5)
pp.virtual_mass_force_type                = 0 # null virtual mass force (0)
pp.lift_force_type                        = pp.dem.consider_lift_force_option # null lift force (0), Saffman (1)
pp.magnus_force_type                      = 0 # null magnus force (0), Rubinow adn Keller (1), Oesterle and Bui Dihn (2)
pp.hydro_torque_type                      = 0 # null hydrodynamic torque (0), Dennis (1)
pp.drag_modifier_type                     = pp.dem.drag_modifier_type # Hayder (2), Chien (3) # problemtype option
pp.drag_porosity_correction_type          = 0 # No correction (0), Richardson and Zaki (1)
pp.interaction_start_time                 = 0
pp.min_fluid_fraction                     = 0.4
pp.initial_drag_force                     = 0.0   # problemtype option
pp.drag_law_slope                         = 0.0   # problemtype option
pp.power_law_tol                          = 0.0
pp.model_over_real_diameter_factor        = 1.0 # not active if similarity_transformation_type = 0
pp.max_pressure_variation_rate_tol        = 1e-3 # for stationary problems, criterion to stop the fluid calculations
pp.time_steps_per_stationarity_step       = 15 # number of fluid time steps between consecutive assessment of stationarity steps
pp.meso_scale_length                      = 0.2 # the radius of the support of the averaging function for homogenization (<=0 for automatic calculation)
pp.shape_factor                           = 0.5 # the density function's maximum over its support's radius (only relevant if coupling_weighing_type == 3)
pp.fluid_domain_volume                    = 2 * math.pi # write down the volume you know it has

# defining and adding imposed porosity fields
pp.fluid_fraction_fields = []
field1 = swim_proc.FluidFractionFieldUtility.LinearField(0.0,
                                                         [0.0, 0.0, 0.0],
                                                         [-1.0, -1.0, 0.15],
                                                         [1.0, 1.0, 0.3])
pp.fluid_fraction_fields.append(field1)

# building lists of variables for which memory is to be allocated
vars_man.ConstructListsOfVariables(pp)
#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    B E G I N S
#_____________________________________________________________________________________________________________________________________

# defining variables to be used
# GID IO IS NOT USING THIS NOW. TO BE REMOVED ONCE THE "PRINT IN POINTS"
# CODE IS NOT USING IT

variables_dictionary = {"PRESSURE"   : PRESSURE,
                        "VELOCITY"   : VELOCITY,
                        "MU"         : MU,         #    MOD.
                        "BUOYANCY"   : BUOYANCY,   #    MOD.
                        "DRAG_FORCE" : DRAG_FORCE,  #    MOD.
                        "LIFT_FORCE" : LIFT_FORCE} #    MOD.

# defining a model part for the fluid part
fluid_model_part = ModelPart("FluidPart")

if "REACTION" in pp.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(REACTION)
if "DISTANCE" in pp.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

#
#
# importing the solvers needed
SolverSettings = pp.FluidSolverConfiguration
solver_module = import_solver(SolverSettings)

#
# importing variables
print('Adding nodal variables to the fluid_model_part')  #     MOD.
sys.stdout.flush()

# caution with breaking up this block (memory allocation)! {
solver_module.AddVariables(fluid_model_part, SolverSettings)
vars_man.AddNodalVariables(fluid_model_part, pp.fluid_vars)  #     MOD.
# }

# introducing input file name
input_file_name = pp.problem_name

# reading the fluid part
model_part_io_fluid = ModelPartIO(input_file_name, True)
model_part_io_fluid.ReadModelPart(fluid_model_part)

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    E N D S
#_____________________________________________________________________________________________________________________________________

# creating utilities from DEM_procedures
DEM_proc = DEM_procedures.Procedures(pp.dem)
DEM_proc.PreProcessModel(pp.dem)

# defining model parts for the balls part and for the DEM-FEM interaction elements

balls_model_part = ModelPart("SolidPart")
clusters_model_part    = ModelPart("Cluster_Part");
rigid_faces_model_part = ModelPart("RigidFace_Part");
DEM_inlet_model_part = ModelPart("DEMInletPart")

import sphere_strategy as DEMSolverStrategy

DEM_proc.AddCommonVariables(balls_model_part, pp.dem)
DEM_proc.AddBallsVariables(balls_model_part, pp.dem)
DEM_proc.AddMpiVariables(balls_model_part)
vars_man.AddNodalVariables(balls_model_part, pp.dem_vars)
DEM_proc.AddCommonVariables(rigid_faces_model_part, pp.dem)
DEM_proc.AddFEMVariables(rigid_faces_model_part, pp.dem)
DEM_proc.AddMpiVariables(rigid_faces_model_part)
DEM_proc.AddCommonVariables(clusters_model_part, pp.dem)
DEM_proc.AddClusterVariables(clusters_model_part, pp.dem)
DEM_proc.AddMpiVariables(clusters_model_part)
DEM_proc.AddCommonVariables(DEM_inlet_model_part, pp.dem)
DEM_proc.AddBallsVariables(DEM_inlet_model_part, pp.dem)
DEM_proc.AddMpiVariables(DEM_inlet_model_part)
vars_man.AddNodalVariables(DEM_inlet_model_part, pp.inlet_vars)

# defining a model part for the mixed part
mixed_model_part = ModelPart("MixedPart")

# reading the balls model part
model_part_io_solid = ModelPartIO(pp.dem.problem_name + "DEM",True)
model_part_io_solid.ReadModelPart(balls_model_part)

# reading the cluster model part
clusters_mp_filename = pp.dem.problem_name + "DEM_Clusters"
model_part_io_clusters = ModelPartIO(clusters_mp_filename)
model_part_io_clusters.ReadModelPart(clusters_model_part)
# reading the fem-dem model part
rigid_face_mp_filename = pp.dem.problem_name + "DEM_FEM_boundary"
model_part_io_solid = ModelPartIO(rigid_face_mp_filename)
model_part_io_solid.ReadModelPart(rigid_faces_model_part)
# reading the dem inlet model part
DEM_Inlet_filename = pp.dem.problem_name + "DEM_Inlet"
model_part_io_demInlet = ModelPartIO(DEM_Inlet_filename)
model_part_io_demInlet.ReadModelPart(DEM_inlet_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
if pp.virtual_mass_force_type > 0:
    buffer_size = 2

else:
    buffer_size = 1

balls_model_part.SetBufferSize(buffer_size)
clusters_model_part.SetBufferSize(buffer_size)
DEM_inlet_model_part.SetBufferSize(buffer_size)

# adding nodal degrees of freedom
DEMSolverStrategy.AddDofs(balls_model_part)
DEMSolverStrategy.AddDofs(clusters_model_part)
DEMSolverStrategy.AddDofs(DEM_inlet_model_part)

# adding extra process info variables
vars_man.AddingDEMProcessInfoVariables(pp, balls_model_part)

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    B E G I N S
#_____________________________________________________________________________________________________________________________________


# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

# adding nodal degrees of freedom
solver_module.AddDofs(fluid_model_part, SolverSettings)

# If Lalplacian form = 2, free all pressure Dofs
# laplacian_form = pp.laplacian_form
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
if(pp.FluidSolverConfiguration.TurbulenceModel == "Spalart-Allmaras"):
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
                   pp.VolumeOutput,
                   pp.GiDPostMode,
                   pp.GiDMultiFileFlag,
                   pp.GiDWriteMeshFlag,
                   pp.GiDWriteConditionsFlag)

if not pp.VolumeOutput:
    cut_list = define_output.DefineCutPlanes()
    gid_io.define_cuts(fluid_model_part, cut_list)

# gid_io.initialize_results(fluid_model_part) # MOD.

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    E N D S
#_____________________________________________________________________________________________________________________________________

import swimming_DEM_gid_output
swimming_DEM_gid_io = swimming_DEM_gid_output.SwimmingDEMGiDOutput(input_file_name,
                                                                   pp.VolumeOutput,
                                                                   pp.GiDPostMode,
                                                                   pp.GiDMultiFileFlag,
                                                                   pp.GiDWriteMeshFlag,
                                                                   pp.GiDWriteConditionsFlag)

swimming_DEM_gid_io.initialize_swimming_DEM_results(balls_model_part, clusters_model_part, rigid_faces_model_part, mixed_model_part)

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    B E G I N S
#_____________________________________________________________________________________________________________________________________

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

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    E N D S
#_____________________________________________________________________________________________________________________________________

# setting fluid's body force to the same as DEM's
if pp.body_force_on_fluid_option:

    for node in fluid_model_part.Nodes:
        node.SetSolutionStepValue(BODY_FORCE_X, 0, pp.dem.GravityX)
        node.SetSolutionStepValue(BODY_FORCE_Y, 0, pp.dem.GravityY)
        node.SetSolutionStepValue(BODY_FORCE_Z, 0, pp.dem.GravityZ)

# coarse-graining: applying changes to the physical properties of the model to adjust for
# the similarity transformation if required (fluid effects only).
swim_proc.ApplySimilarityTransformations(fluid_model_part, pp.similarity_transformation_type, pp.model_over_real_diameter_factor)

# creating a Post Utils object that executes several post-related tasks
post_utils = swim_proc.PostUtils(swimming_DEM_gid_io,
                                               pp,
                                               fluid_model_part,
                                               balls_model_part,
                                               clusters_model_part,
                                               rigid_faces_model_part,
                                               mixed_model_part)

# creating an IOTools object to perform other printing tasks
io_tools = swim_proc.IOTools(pp)

# creating a projection module for the fluid-DEM coupling
h_min = 0.01
n_balls = 1
fluid_volume = 10
pp.n_particles_in_depth = int(math.sqrt(n_balls / fluid_volume)) # only relevant in 2D problems
# creating a physical calculations module to analyse the DEM model_part
dem_physics_calculator = SphericElementGlobalPhysicsCalculator(balls_model_part)

if pp.projection_module_option:

    if pp.meso_scale_length <= 0.0  and balls_model_part.NumberOfElements(0) > 0:
        biggest_size = 2 * dem_physics_calculator.CalculateMaxNodalVariable(balls_model_part, RADIUS)
        pp.meso_scale_length = 20 * biggest_size

    elif balls_model_part.NumberOfElements(0) == 0:
        pp.meso_scale_length = 1.0

    projection_module = CFD_DEM_coupling.ProjectionModule(fluid_model_part, balls_model_part, rigid_faces_model_part, domain_size, pp)
    projection_module.UpdateDatabase(h_min)

# creating a custom functions calculator for the implementation of additional custom functions
custom_functions_tool = swim_proc.FunctionsCalculator(pp)

# creating a stationarity assessment tool
stationarity_tool = swim_proc.StationarityAssessmentTool(pp.max_pressure_variation_rate_tol, custom_functions_tool)

# creating a debug tool
dem_volume_tool = swim_proc.ProjectionDebugUtils(pp.fluid_domain_volume, fluid_model_part, balls_model_part, custom_functions_tool)

# creating a CreatorDestructor object, responsible for any adding or removing of elements during the simulation
creator_destructor = ParticleCreatorDestructor()

# setting up a bounding box for the DEM balls (it is used for erasing remote balls)
DEM_proc.SetBoundingBox(balls_model_part, clusters_model_part, rigid_faces_model_part, creator_destructor)

# creating a Solver object for the DEM part. It contains the sequence of function calls necessary for the evolution of the DEM system at every time step
dem_solver = DEMSolverStrategy.ExplicitStrategy(balls_model_part, rigid_faces_model_part, clusters_model_part, DEM_inlet_model_part, creator_destructor, pp.dem) 
# Initializing the DEM solver (must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part)
dem_solver.Initialize()

# creating a distance calculation process for the embedded technology
# (used to calculate elemental distances defining the structure embedded in the fluid mesh)
if (pp.embedded_option):
    calculate_distance_process = CalculateSignedDistanceTo3DSkinProcess(rigid_faces_model_part, fluid_model_part)
    calculate_distance_process.Execute()

# constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation

if pp.inlet_option:
    
    vars_man.AddingDEMProcessInfoVariables(pp, DEM_inlet_model_part)
    max_node_Id = creator_destructor.FindMaxNodeIdInModelPart(balls_model_part)
    max_FEM_node_Id = creator_destructor.FindMaxNodeIdInModelPart(rigid_faces_model_part)
    max_fluid_node_Id = swim_proc.FindMaxNodeIdInFLuid(fluid_model_part)

    if max_FEM_node_Id > max_node_Id:
        max_node_Id = max_FEM_node_Id
    if max_fluid_node_Id > max_node_Id:
        max_node_Id = max_FEM_node_Id
    
    creator_destructor.SetMaxNodeId(max_node_Id)      

    for properties in DEM_inlet_model_part.Properties:
        DiscontinuumConstitutiveLawString = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME];
        DiscontinuumConstitutiveLaw = globals().get(DiscontinuumConstitutiveLawString)()
        DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties)

    # constructiong the inlet and intializing it (must be done AFTER the balls_model_part Initialize)
    DEM_inlet = DEM_Inlet(DEM_inlet_model_part)
    DEM_inlet.InitializeDEM_Inlet(balls_model_part, creator_destructor, pp.dem_inlet_element_type)

# creating problem directories
directories = ['results']

if pp.make_results_directories_option:
    dir_path_dictionary = io_tools.CreateProblemDirectories(current_script_path, directories)

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    B E G I N S
#_____________________________________________________________________________________________________________________________________

# renumerating IDs if required

# swim_proc.RenumberNodesIdsToAvoidRepeating(fluid_model_part, balls_model_part, rigid_faces_model_part)

# Stepping and time settings

Dt = pp.Dt
Nsteps = pp.nsteps
final_time = pp.max_time
output_time = pp.output_time

time = pp.Start_time
out = Dt
step = 0

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    E N D S
#_____________________________________________________________________________________________________________________________________

if pp.flow_in_porous_medium_option:
    fluid_frac_util = swim_proc.FluidFractionFieldUtility(fluid_model_part, pp.min_fluid_fraction)

    for field in pp.fluid_fraction_fields:
        fluid_frac_util.AppendLinearField(field)

    fluid_frac_util.AddFluidFractionField()

if pp.flow_in_porous_DEM_medium_option:
    swim_proc.FixModelPart(balls_model_part)

# choosing the directory in which we want to work (print to)

if pp.make_results_directories_option:
    os.chdir(dir_path_dictionary['results'])

def yield_DEM_time(current_time, current_time_plus_increment, delta_time):
    current_time += delta_time

    while current_time < current_time_plus_increment:
        yield current_time
        current_time += delta_time

    current_time = current_time_plus_increment
    yield current_time

######################################################################################################################################

#                      I N I T I A L I Z I N G    T I M E    L O O P     ...   ( M I X E D    F L U I D / D E M    B L O C K )

######################################################################################################################################

# setting up loop counters
embedded_counter          = swim_proc.Counter(1, 3, pp.embedded_option)  # MA: because I think DISTANCE,1 (from previous time step) is not calculated correctly for step=1
DEM_to_fluid_counter      = swim_proc.Counter(1, 1, pp.coupling_level_type == 1)
pressure_gradient_counter = swim_proc.Counter(1, 1, pp.projection_module_option)
stationarity_counter      = swim_proc.Counter(pp.time_steps_per_stationarity_step, 1, pp.stationary_problem_option)
print_counter             = swim_proc.Counter(1, 1, out >= output_time)
debug_info_counter        = swim_proc.Counter(pp.debug_tool_cycle , 1, pp.print_debug_info_option)
particles_results_counter = swim_proc.Counter(pp.print_particles_results_cycle, 1, pp.print_particles_results_option)

DEM_step     = 0      # necessary to get a good random insertion of particles   # relevant to the stationarity assessment tool
time_dem     = 0.0
Dt_DEM       = balls_model_part.ProcessInfo.GetValue(DELTA_TIME)
rigid_faces_model_part.ProcessInfo[DELTA_TIME] = Dt_DEM
clusters_model_part.ProcessInfo[DELTA_TIME] = Dt_DEM
stationarity = False

mesh_motion = DEMFEMUtilities()
swim_proc.InitializeVariablesWithNonZeroValues(fluid_model_part, balls_model_part) # all variables are set to 0 by default

while time <= final_time:

    time = time + Dt
    step += 1
    fluid_model_part.CloneTimeStep(time)
    print("\n", "TIME = ", time)
    sys.stdout.flush()

    if pp.coupling_scheme_type == "UpdatedDEM":
        time_final_DEM_substepping = time + Dt

    else:
        time_final_DEM_substepping = time
        time_dem = time - Dt + Dt_DEM

    # walls movement
    mesh_motion.MoveAllMeshes(rigid_faces_model_part, time)

    # calculating elemental distances defining the structure embedded in the fluid mesh
    if pp.embedded_option:
        calculate_distance_process.Execute()

    if embedded_counter.Tick():
        embedded.ApplyEmbeddedBCsToFluid(fluid_model_part)
        embedded.ApplyEmbeddedBCsToBalls(balls_model_part, pp.dem)

    # solving the fluid part

    if step >= 3 and not stationarity:

        print("Solving Fluid... (", fluid_model_part.NumberOfElements(0), " elements )")
        sys.stdout.flush()

        fluid_solver.Solve()

    # assessing stationarity

        if stationarity_counter.Tick():
            print("Assessing Stationarity...")
            stationarity = stationarity_tool.Assess(fluid_model_part)
            sys.stdout.flush()

    # printing if required

    if particles_results_counter.Tick():
        # eliminating remote balls
        
        if pp.dem.BoundingBoxOption == "ON":
            creator_destructor.DestroyParticlesOutsideBoundingBox(balls_model_part)
            
        io_tools.PrintParticlesResults(pp.variables_to_print_in_file, time, balls_model_part)
        graph_printer.PrintGraphs(time)
        PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

    if output_time <= out and pp.coupling_scheme_type == "UpdatedDEM":

        if pp.projection_module_option:
            projection_module.ComputePostProcessResults(balls_model_part.ProcessInfo)

        post_utils.Writeresults(time)
        out = 0

    # solving the DEM part
    pressure_gradient_counter.Deactivate(time < pp.interaction_start_time)

    if pressure_gradient_counter.Tick():
        custom_functions_tool.CalculatePressureGradient(fluid_model_part)

    print("Solving DEM... (", balls_model_part.NumberOfElements(0), " elements)")
    sys.stdout.flush()
    first_dem_iter = True

    for time_dem in yield_DEM_time(time_dem, time_final_DEM_substepping, Dt_DEM):

        DEM_step += 1   # this variable is necessary to get a good random insertion of particles

        balls_model_part.ProcessInfo[TIME_STEPS] = DEM_step

        # applying fluid-to-DEM coupling if required

        if time >= pp.interaction_start_time and pp.projection_module_option and (pp.project_at_every_substep_option or first_dem_iter):

            if pp.coupling_scheme_type == "UpdatedDEM":
                projection_module.ProjectFromNewestFluid()

            else:
                projection_module.ProjectFromFluid((time_final_DEM_substepping - time_dem) / Dt)

        # performing the time integration of the DEM part

        balls_model_part.CloneTimeStep(time_dem)

        if not pp.flow_in_porous_DEM_medium_option: # in porous flow particles remain static
            dem_solver.Solve()

        # adding DEM elements by the inlet

        if pp.inlet_option:
            DEM_inlet.CreateElementsFromInletMesh(balls_model_part, DEM_inlet_model_part, creator_destructor, pp.dem_inlet_element_type)  # After solving, to make sure that neighbours are already set.        
        
        # eliminating remote balls

        creator_destructor.DestroyParticlesOutsideBoundingBox(balls_model_part)
        
        first_dem_iter = False

    # measuring mean velocities in a certain control volume (the 'velocity trap')

    if pp.velocity_trap_option:
        post_utils.ComputeMeanVelocitiesinTrap("Average_Velocity", time_dem)

    # applying DEM-to-fluid coupling

    if DEM_to_fluid_counter.Tick() and time >= pp.interaction_start_time:
        print("Projecting from particles to the fluid...")
        sys.stdout.flush()
        projection_module.ProjectFromParticles()

    # coupling checks (debugging)

    if debug_info_counter.Tick():
        dem_volume_tool.UpdateDataAndPrint(pp.fluid_domain_volume)

    # printing if required

    if particles_results_counter.Tick():
        
        # eliminating remote balls
        
        if pp.dem.BoundingBoxOption == "ON":
            creator_destructor.DestroyParticlesOutsideBoundingBox(balls_model_part)
            
        io_tools.PrintParticlesResults(pp.variables_to_print_in_file, time, balls_model_part)
        graph_printer.PrintGraphs(time)
        PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

    if output_time <= out and pp.coupling_scheme_type == "UpdatedFluid":

        if pp.projection_module_option:
            projection_module.ComputePostProcessResults(balls_model_part.ProcessInfo)

        post_utils.Writeresults(time)
        out = 0

    out = out + Dt

swimming_DEM_gid_io.finalize_results()

print("CALCULATIONS FINISHED. THE SIMULATION ENDED SUCCESSFULLY.")
sys.stdout.flush()

for i in drag_file_output_list:
    i.close()
