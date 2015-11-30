# This script contains an algorithm that models fluid-particle interaction.
# It combines two parts: a FEM model for the fluid and a DEM model for the particles.
# It has been conceived by adding the DEM part and the interaction on top of an original fluid-only script (see kratos/applications/FluidDynamicswApplication)
# Some parts of the original fluid script have been kept practically untouched and are clearly marked.
# Whenever a minor modification has been made on one of these parts, the corresponding line is indicated with a comment: # MOD.

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Python imports
import time as timer
import os
import sys
import math

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

sys.path.insert(0,'')
# DEM Application
import DEM_explicit_solver_var as DEM_parameters
import DEM_procedures

# import the configuration data as read from the GiD
print(os.getcwd())
import ProjectParameters as pp # MOD
import define_output

# setting the domain size for the problem to be solved
domain_size = pp.domain_size

import swimming_dem_parameters
import swimming_DEM_procedures as swim_proc
import CFD_DEM_coupling
import variables_management as vars_man
import embedded

# listing project parameters (to be put in problem type)
pp.swim                                   = swimming_dem_parameters

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI. For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ:
    # Kratos MPI
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *

    # DEM Application MPI
    import DEM_procedures_mpi as DEM_procedures
    import DEM_material_test_script_mpi as DEM_material_test_script
    print("Running under MPI...........")
else:
    # DEM Application
    import DEM_procedures
    import DEM_material_test_script
    print("Running under OpenMP........")

# TO_DO: Ugly fix. Change it. I don't like this to be in the main...
# Strategy object

import sphere_strategy as SolverStrategy

pp.swim.coupling_level_type = DEM_parameters.project_from_particles_option
pp.swim.lift_force_type = DEM_parameters.consider_lift_force_option
pp.swim.drag_modifier_type = DEM_parameters.drag_modifier_type
pp.swim.fluid_domain_volume                    = 2 * math.pi # write down the volume you know it has

##############################################################################
#                                                                            #
#    INITIALIZE                                                              #
#                                                                            #
##############################################################################

# Import utilities from models
procedures    = DEM_procedures.Procedures(DEM_parameters)
demio         = DEM_procedures.DEMIo(DEM_parameters)
report        = DEM_procedures.Report()
parallelutils = DEM_procedures.ParallelUtils()
materialTest  = DEM_procedures.MaterialTest()
 
# Set the print function TO_DO: do this better...
KRATOSprint   = procedures.KRATOSprint

# Preprocess the model
procedures.PreProcessModel(DEM_parameters)

# Prepare modelparts
spheres_model_part    = ModelPart("SpheresPart")
rigid_face_model_part = ModelPart("RigidFace_Part")
mixed_model_part      = ModelPart("Mixed_Part")
cluster_model_part    = ModelPart("Cluster_Part")
DEM_inlet_model_part  = ModelPart("DEMInletPart")
mapping_model_part    = ModelPart("Mappingmodel_part")
contact_model_part    = ""

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
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    E N D S
#_____________________________________________________________________________________________________________________________________

# Constructing a utilities objects
creator_destructor = ParticleCreatorDestructor()
dem_fem_search = DEM_FEM_Search()

# Creating a solver object and set the search strategy

solver = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, dem_fem_search, DEM_parameters)

# Add variables
procedures.AddCommonVariables(spheres_model_part, DEM_parameters)
procedures.AddSpheresVariables(spheres_model_part, DEM_parameters)
procedures.AddMpiVariables(spheres_model_part)
solver.AddAdditionalVariables(spheres_model_part, DEM_parameters)
procedures.AddCommonVariables(cluster_model_part, DEM_parameters)
procedures.AddClusterVariables(cluster_model_part, DEM_parameters)
procedures.AddMpiVariables(cluster_model_part)
procedures.AddCommonVariables(DEM_inlet_model_part, DEM_parameters)
procedures.AddSpheresVariables(DEM_inlet_model_part, DEM_parameters)
solver.AddAdditionalVariables(DEM_inlet_model_part, DEM_parameters)  
procedures.AddCommonVariables(rigid_face_model_part, DEM_parameters)
procedures.AddRigidFaceVariables(rigid_face_model_part, DEM_parameters)
procedures.AddMpiVariables(rigid_face_model_part)
vars_man.AddNodalVariables(spheres_model_part, pp.dem_vars)
vars_man.AddNodalVariables(rigid_face_model_part, pp.rigid_faces_vars)
vars_man.AddNodalVariables(DEM_inlet_model_part, pp.inlet_vars)

# defining a model part for the mixed part
mixed_model_part = ModelPart("MixedPart")


# reading the balls model part
spheres_mp_filename   = DEM_parameters.problem_name + "DEM"
model_part_io_spheres = ModelPartIO(spheres_mp_filename)

# performing the initial partition
[model_part_io_spheres, spheres_model_part, MPICommSetup] = parallelutils.PerformInitialPartition(spheres_model_part, model_part_io_spheres, spheres_mp_filename)
model_part_io_spheres.ReadModelPart(spheres_model_part)

# reading the cluster model part
clusters_mp_filename = DEM_parameters.problem_name + "DEM_Clusters"
model_part_io_clusters = ModelPartIO(clusters_mp_filename)
model_part_io_clusters.ReadModelPart(cluster_model_part)
# reading the fem-dem model part
rigid_face_mp_filename = DEM_parameters.problem_name + "DEM_FEM_boundary"
model_part_io_fem = ModelPartIO(rigid_face_mp_filename)
model_part_io_fem.ReadModelPart(rigid_face_model_part)
# reading the dem inlet model part
DEM_Inlet_filename = DEM_parameters.problem_name + "DEM_Inlet"
model_part_io_demInlet = ModelPartIO(DEM_Inlet_filename)
model_part_io_demInlet.ReadModelPart(DEM_inlet_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
if pp.swim.virtual_mass_force_type  > 0:
    buffer_size = 2

else:
    buffer_size = 1

spheres_model_part.SetBufferSize(buffer_size)
cluster_model_part.SetBufferSize(buffer_size)
DEM_inlet_model_part.SetBufferSize(buffer_size)
rigid_face_model_part.SetBufferSize(1)

# adding nodal degrees of freedom
SolverStrategy.AddDofs(spheres_model_part)
SolverStrategy.AddDofs(cluster_model_part)
SolverStrategy.AddDofs(DEM_inlet_model_part)

# adding extra process info variables
vars_man.AddingDEMProcessInfoVariables(pp, spheres_model_part)

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

swimming_DEM_gid_io.initialize_swimming_DEM_results(spheres_model_part, cluster_model_part, rigid_face_model_part, mixed_model_part)

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
#import point_graph_printer

#output_nodes_list = define_output.DefineOutputPoints()
#graph_printer = point_graph_printer.PrintGraphPrinter(
    #output_nodes_list,
    #fluid_model_part,
    #variables_dictionary,
    #domain_size)

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    E N D S
#_____________________________________________________________________________________________________________________________________

# setting fluid's body force to the same as DEM's
if pp.swim.body_force_on_fluid_option:

    for node in fluid_model_part.Nodes:
        node.SetSolutionStepValue(BODY_FORCE_X, 0, DEM_parameters.GravityX)
        node.SetSolutionStepValue(BODY_FORCE_Y, 0, DEM_parameters.GravityY)
        node.SetSolutionStepValue(BODY_FORCE_Z, 0, DEM_parameters.GravityZ)

# coarse-graining: applying changes to the physical properties of the model to adjust for
# the similarity transformation if required (fluid effects only).
swim_proc.ApplySimilarityTransformations(fluid_model_part, pp.swim.similarity_transformation_type , pp.swim.model_over_real_diameter_factor )

# creating a Post Utils object that executes several post-related tasks
post_utils = swim_proc.PostUtils(swimming_DEM_gid_io,
                                               pp,
                                               fluid_model_part,
                                               spheres_model_part,
                                               cluster_model_part,
                                               rigid_face_model_part,
                                               mixed_model_part)

# creating an IOTools object to perform other printing tasks
io_tools = swim_proc.IOTools(pp)

# creating a projection module for the fluid-DEM coupling
h_min = 0.01
n_balls = 1
fluid_volume = 10
pp.swim.n_particles_in_depth = int(math.sqrt(n_balls / fluid_volume)) # only relevant in 2D problems
# creating a physical calculations module to analyse the DEM model_part
dem_physics_calculator = SphericElementGlobalPhysicsCalculator(spheres_model_part)

if pp.swim.projection_module_option:

    if pp.swim.meso_scale_length  <= 0.0 and spheres_model_part.NumberOfElements(0) > 0:
        biggest_size = 2 * dem_physics_calculator.CalculateMaxNodalVariable(spheres_model_part, RADIUS)
        pp.swim.meso_scale_length  = 20 * biggest_size

    elif spheres_model_part.NumberOfElements(0) == 0:
        pp.swim.meso_scale_length  = 1.0

    projection_module = CFD_DEM_coupling.ProjectionModule(fluid_model_part, spheres_model_part, rigid_face_model_part, domain_size, pp)
    projection_module.UpdateDatabase(h_min)

# creating a custom functions calculator for the implementation of additional custom functions
custom_functions_tool = swim_proc.FunctionsCalculator(pp)

# creating a stationarity assessment tool
stationarity_tool = swim_proc.StationarityAssessmentTool(pp.swim.max_pressure_variation_rate_tol , custom_functions_tool)

# creating a debug tool
dem_volume_tool = swim_proc.ProjectionDebugUtils(pp.swim.fluid_domain_volume, fluid_model_part, spheres_model_part, custom_functions_tool)

# creating a CreatorDestructor object, responsible for any adding or removing of elements during the simulation
creator_destructor = ParticleCreatorDestructor()
dem_fem_search = DEM_FEM_Search()

# setting up a bounding box for the DEM balls (it is used for erasing remote balls)
procedures.SetBoundingBox(spheres_model_part, cluster_model_part, rigid_face_model_part, creator_destructor)

# creating a Solver object for the DEM part. It contains the sequence of function calls necessary for the evolution of the DEM system at every time step
dem_solver = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, dem_fem_search, DEM_parameters)

# Initializing the DEM solver (must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part)
dem_solver.search_strategy = parallelutils.GetSearchStrategy(dem_solver, spheres_model_part)
dem_solver.Initialize()

# creating a distance calculation process for the embedded technology
# (used to calculate elemental distances defining the structure embedded in the fluid mesh)
if (pp.swim.embedded_option):
    calculate_distance_process = CalculateSignedDistanceTo3DSkinProcess(rigid_face_model_part, fluid_model_part)
    calculate_distance_process.Execute()

# constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation

if pp.swim.inlet_option:
    
    vars_man.AddingDEMProcessInfoVariables(pp, DEM_inlet_model_part)
    max_DEM_node_Id = creator_destructor.FindMaxNodeIdInModelPart(spheres_model_part)
    max_FEM_node_Id = creator_destructor.FindMaxNodeIdInModelPart(rigid_face_model_part)
    max_fluid_node_Id = swim_proc.FindMaxNodeIdInFLuid(fluid_model_part)
    max_node_Id = max(max_DEM_node_Id, max_FEM_node_Id, max_fluid_node_Id)
    
    creator_destructor.SetMaxNodeId(max_node_Id)      

    for properties in DEM_inlet_model_part.Properties:
        DiscontinuumConstitutiveLawString = properties[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME];
        DiscontinuumConstitutiveLaw = globals().get(DiscontinuumConstitutiveLawString)()
        DiscontinuumConstitutiveLaw.SetConstitutiveLawInProperties(properties)

    # constructiong the inlet and intializing it (must be done AFTER the spheres_model_part Initialize)
    DEM_inlet = DEM_Inlet(DEM_inlet_model_part)
    DEM_inlet.InitializeDEM_Inlet(spheres_model_part, creator_destructor)

# Creating necessary directories
main_path = os.getcwd()
[post_path,data_and_results,graphs_path,MPI_results] = procedures.CreateDirectories(str(main_path),str(DEM_parameters.problem_name))

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    B E G I N S
#_____________________________________________________________________________________________________________________________________

# renumerating IDs if required

# swim_proc.RenumberNodesIdsToAvoidRepeating(fluid_model_part, spheres_model_part, rigid_face_model_part)

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

if pp.swim.flow_in_porous_medium_option:
    fluid_frac_util = swim_proc.FluidFractionFieldUtility(fluid_model_part, pp.swim.min_fluid_fraction )

    for field in pp.fluid_fraction_fields:
        fluid_frac_util.AppendLinearField(field)

    fluid_frac_util.AddFluidFractionField()

if pp.swim.flow_in_porous_DEM_medium_option:
    swim_proc.FixModelPart(spheres_model_part)

# choosing the directory in which we want to work (print to)

os.chdir(post_path)

def yield_DEM_time(current_time, current_time_plus_increment, delta_time):
    current_time += delta_time

    tolerance = 0.0001
    while current_time < (current_time_plus_increment - tolerance * delta_time):
        yield current_time
        current_time += delta_time

    current_time = current_time_plus_increment
    yield current_time

######################################################################################################################################

#                      I N I T I A L I Z I N G    T I M E    L O O P     ...   ( M I X E D    F L U I D / D E M    B L O C K )

######################################################################################################################################

# setting up loop counters
embedded_counter          = swim_proc.Counter(1, 3, pp.swim.embedded_option)  # MA: because I think DISTANCE,1 (from previous time step) is not calculated correctly for step=1
DEM_to_fluid_counter      = swim_proc.Counter(1, 1, pp.swim.coupling_level_type  == 1)
pressure_gradient_counter = swim_proc.Counter(1, 1, pp.swim.projection_module_option)
stationarity_counter      = swim_proc.Counter(pp.swim.time_steps_per_stationarity_step , 1, pp.swim.stationary_problem_option)
print_counter             = swim_proc.Counter(1, 1, out >= output_time)
debug_info_counter        = swim_proc.Counter(pp.swim.debug_tool_cycle, 1, pp.swim.print_debug_info_option)
particles_results_counter = swim_proc.Counter(pp.swim.print_particles_results_cycle , 1, pp.swim.print_particles_results_option)

#G
#fluid_model_part.AddNodalSolutionStepVariable(BODY_FORCE)

#for node in fluid_model_part.Nodes:
    #node.SetSolutionStepValue(BODY_FORCE_X, 0, 0.0)
    #node.SetSolutionStepValue(BODY_FORCE_Y, 0, 0.0)
    #node.SetSolutionStepValue(BODY_FORCE_Z, 0, - 9.81)
#import meshing_utils

#meshing_utils.ModifyRadiusesGivenPorosity(model_part, porosity, tol)
#Z

DEM_step     = 0      # necessary to get a good random insertion of particles   # relevant to the stationarity assessment tool
time_dem     = 0.0
Dt_DEM       = spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
rigid_face_model_part.ProcessInfo[DELTA_TIME] = Dt_DEM
cluster_model_part.ProcessInfo[DELTA_TIME] = Dt_DEM
stationarity = False

mesh_motion = DEMFEMUtilities()
swim_proc.InitializeVariablesWithNonZeroValues(fluid_model_part, spheres_model_part) # all variables are set to 0 by default

while time <= final_time:

    time = time + Dt
    step += 1
    fluid_model_part.CloneTimeStep(time)
    print("\n", "TIME = ", time)
    sys.stdout.flush()

    if pp.swim.coupling_scheme_type  == "UpdatedDEM":
        time_final_DEM_substepping = time + Dt

    else:
        time_final_DEM_substepping = time

    # walls movement
    mesh_motion.MoveAllMeshes(rigid_face_model_part, time, Dt)

    # calculating elemental distances defining the structure embedded in the fluid mesh
    if pp.swim.embedded_option:
        calculate_distance_process.Execute()

    if embedded_counter.Tick():
        embedded.ApplyEmbeddedBCsToFluid(fluid_model_part)
        embedded.ApplyEmbeddedBCsToBalls(spheres_model_part, DEM_parameters)

    # solving the fluid part

    if step >= 3 and not stationarity:

        print("Solving Fluid... (", fluid_model_part.NumberOfElements(0), "elements )")
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
        
        #if pp.dem.BoundingBoxOption == "ON":
        #    creator_destructor.DestroyParticlesOutsideBoundingBox(spheres_model_part)
            
        io_tools.PrintParticlesResults(pp.variables_to_print_in_file, time, spheres_model_part)
        graph_printer.PrintGraphs(time)
        PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

    if output_time <= out and pp.swim.coupling_scheme_type == "UpdatedDEM":

        if pp.swim.projection_module_option:
            projection_module.ComputePostProcessResults(spheres_model_part.ProcessInfo)

        post_utils.Writeresults(time)
        out = 0

    # solving the DEM part
    pressure_gradient_counter.Deactivate(time < pp.swim.interaction_start_time)

    if pressure_gradient_counter.Tick():
        custom_functions_tool.CalculatePressureGradient(fluid_model_part)

    print("Solving DEM... (", spheres_model_part.NumberOfElements(0), "elements )")
    sys.stdout.flush()
    first_dem_iter = True

    for time_dem in yield_DEM_time(time_dem, time_final_DEM_substepping, Dt_DEM):

        DEM_step += 1   # this variable is necessary to get a good random insertion of particles

        spheres_model_part.ProcessInfo[TIME_STEPS] = DEM_step

        # applying fluid-to-DEM coupling if required

        if time >= pp.swim.interaction_start_time and pp.swim.projection_module_option and (pp.swim.project_at_every_substep_option or first_dem_iter):

            if pp.swim.coupling_scheme_type == "UpdatedDEM":
                projection_module.ProjectFromNewestFluid()

            else:
                projection_module.ProjectFromFluid((time_final_DEM_substepping - time_dem) / Dt)

        # performing the time integration of the DEM part

        spheres_model_part.CloneTimeStep(time_dem)

        if not pp.swim.flow_in_porous_DEM_medium_option: # in porous flow particles remain static
            dem_solver.Solve()

        # adding DEM elements by the inlet

        if pp.swim.inlet_option:
            DEM_inlet.CreateElementsFromInletMesh(spheres_model_part, cluster_model_part, creator_destructor)  # After solving, to make sure that neighbours are already set.        
        
        # eliminating remote balls

        #creator_destructor.DestroyParticlesOutsideBoundingBox(spheres_model_part)
        
        first_dem_iter = False

    # measuring mean velocities in a certain control volume (the 'velocity trap')

    if pp.velocity_trap_option:
        post_utils.ComputeMeanVelocitiesinTrap("Average_Velocity", time_dem)

    # applying DEM-to-fluid coupling

    if DEM_to_fluid_counter.Tick() and time >= pp.swim.interaction_start_time:
        print("Projecting from particles to the fluid...")
        sys.stdout.flush()
        projection_module.ProjectFromParticles()

    # coupling checks (debugging)

    if debug_info_counter.Tick():
        dem_volume_tool.UpdateDataAndPrint(pp.swim.fluid_domain_volume)

    # printing if required

    if particles_results_counter.Tick():
        
        # eliminating remote balls
        
        if pp.dem.BoundingBoxOption == "ON":
            creator_destructor.DestroyParticlesOutsideBoundingBox(spheres_model_part)
            
        io_tools.PrintParticlesResults(pp.variables_to_print_in_file, time, spheres_model_part)
        graph_printer.PrintGraphs(time)
        PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

    if output_time <= out and pp.swim.coupling_scheme_type == "UpdatedFluid":

        if pp.swim.projection_module_option:
            projection_module.ComputePostProcessResults(spheres_model_part.ProcessInfo)

        post_utils.Writeresults(time)
        out = 0

    out = out + Dt

swimming_DEM_gid_io.finalize_results()

print("CALCULATIONS FINISHED. THE SIMULATION ENDED SUCCESSFULLY.")
sys.stdout.flush()

for i in drag_file_output_list:
    i.close()
