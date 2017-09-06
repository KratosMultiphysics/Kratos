# This script contains an algorithm that models fluid-particle interaction.
# It combines two parts: a FEM model for the fluid and a DEM model for the particles.
# It has been conceived by adding the DEM part and the interaction on top of an original fluid-only script (see kratos/applications/FluidDynamicswApplication)
# Some parts of the original fluid script have been kept practically untouched and are clearly marked.
# Whenever a minor modification has been made on one of these parts, the corresponding line is indicated with a comment: # MOD.

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Python imports
import time as timer
init_time = timer.time()
import os
os.system('cls' if os.name == 'nt' else 'clear')
import sys
sys.path.append("/home/gcasas/kratos")
import numpy as np

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.path_to_console_out_file = "console_output.txt"
        self.log = open(self.path_to_console_out_file, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass

def num_type(value):
  try:
    int(str(value))
    return 'int'
  except:
    try:
        float(str(value))
        return 'float'
    except:
        return None


sys.stdout = Logger()
import math

#G
from matplotlib import pyplot as plt
import pylab
#Z
simulation_start_time = timer.clock()

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

sys.path.insert(0,'')
# DEM Application
import DEM_explicit_solver_var as DEM_parameters

# import the configuration data as read from the GiD
import ProjectParameters as pp # MOD
import define_output

# setting the domain size for the problem to be solved
domain_size = pp.domain_size

import swimming_dem_parameters
import swimming_DEM_procedures as swim_proc
import CFD_DEM_coupling
import variables_management as vars_man
import embedded

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

import swimming_sphere_strategy as SolverStrategy

DEM_parameters.fluid_domain_volume                    = 0.5 ** 2 * 2 * math.pi # write down the volume you know it has

##############################################################################
#                                                                            #
#    INITIALIZE                                                              #
#                                                                            #
##############################################################################
n_div = 20
material_derivative_type = 2
laplacian_type = 2
file_name, size_parameter, material_derivative_type, laplacian_type = sys.argv
#G
pp.CFD_DEM = DEM_parameters
pp.CFD_DEM.fluid_already_calculated = 0
pp.CFD_DEM.recovery_echo_level = 1
pp.CFD_DEM.gradient_calculation_type = 5
pp.CFD_DEM.pressure_grad_recovery_type = 1
pp.CFD_DEM.store_full_gradient = 1
pp.CFD_DEM.laplacian_calculation_type = int(laplacian_type)
pp.CFD_DEM.do_search_neighbours = False
pp.CFD_DEM.material_acceleration_calculation_type = int(material_derivative_type)
pp.CFD_DEM.faxen_force_type = 0
pp.CFD_DEM.vorticity_calculation_type = 0
pp.CFD_DEM.print_FLUID_VEL_PROJECTED_RATE_option = 0
pp.CFD_DEM.print_MATERIAL_FLUID_ACCEL_PROJECTED_option = True
pp.CFD_DEM.basset_force_type = 4
pp.CFD_DEM.print_BASSET_FORCE_option = 1
pp.CFD_DEM.basset_force_integration_type = 2
pp.CFD_DEM.n_init_basset_steps = 0
pp.CFD_DEM.time_steps_per_quadrature_step = 1
pp.CFD_DEM.delta_time_quadrature = pp.CFD_DEM.time_steps_per_quadrature_step * pp.CFD_DEM.MaxTimeStep
pp.CFD_DEM.quadrature_order = 2
pp.CFD_DEM.time_window = 0.01
pp.CFD_DEM.number_of_exponentials = 8
pp.CFD_DEM.number_of_quadrature_steps_in_window = int(pp.CFD_DEM.time_window / pp.CFD_DEM.delta_time_quadrature)
pp.CFD_DEM.print_steps_per_plot_step = 1
pp.CFD_DEM.PostCationConcentration = False
pp.CFD_DEM.do_impose_flow_from_field = False
pp.CFD_DEM.print_MATERIAL_ACCELERATION_option = True
pp.CFD_DEM.print_FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED_option = False
number_of_vectors_to_be_kept_in_memory = pp.CFD_DEM.time_window / pp.CFD_DEM.MaxTimeStep * pp.CFD_DEM.time_steps_per_quadrature_step + pp.CFD_DEM.number_of_exponentials
pp.CFD_DEM.print_VORTICITY_option = 1
pp.CFD_DEM.print_MATERIAL_ACCELERATION_option = True
pp.CFD_DEM.print_VELOCITY_GRADIENT_option = 1
print('\nNumber of vectors to be kept in memory: ', number_of_vectors_to_be_kept_in_memory)
# Making the fluid step an exact multiple of the DEM step
pp.Dt = int(pp.Dt / pp.CFD_DEM.MaxTimeStep) * pp.CFD_DEM.MaxTimeStep
# Creating a code for the used input variables
run_code = '_ndiv_' + str(size_parameter) + '_mat_deriv_type_' + str(material_derivative_type) + '_lapl_type_' + str(laplacian_type)
#Z

# Creating swimming DEM procedures
procedures    = DEM_procedures.Procedures(DEM_parameters)

# Creating necessary directories
main_path = os.getcwd()

[post_path, data_and_results, graphs_path, MPI_results] = procedures.CreateDirectories(str(main_path), str(DEM_parameters.problem_name), run_code)

# Import utilities from models

demio         = DEM_procedures.DEMIo(DEM_parameters, post_path)
report        = DEM_procedures.Report()
parallelutils = DEM_procedures.ParallelUtils()
materialTest  = DEM_procedures.MaterialTest()

# Moving to the recently created folder
os.chdir(main_path)
swim_proc.CopyInputFilesIntoFolder(main_path, post_path)

# Set the print function TO_DO: do this better...
KRATOSprint   = procedures.KRATOSprint

# Preprocess the model
procedures.PreProcessModel(DEM_parameters)

# Prepare modelparts
spheres_model_part    = ModelPart("SpheresPart")
rigid_face_model_part = ModelPart("RigidFacePart")
cluster_model_part    = ModelPart("ClusterPart")
DEM_inlet_model_part  = ModelPart("DEMInletPart")
mapping_model_part    = ModelPart("MappingPart")
contact_model_part    = ModelPart("ContactPart")
all_model_parts = DEM_procedures.SetOfModelParts(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, mapping_model_part, contact_model_part)

# defining and adding imposed porosity fields
pp.fluid_fraction_fields = []
field1 = swim_proc.FluidFractionFieldUtility.LinearField(0.0,
                                                         [0.0, 0.0, 0.0],
                                                         [-1.0, -1.0, 0.15],
                                                         [1.0, 1.0, 0.3])
pp.fluid_fraction_fields.append(field1)

# building lists of variables for which memory is to be allocated
# TEMPORARY, HORRIBLE !!!
pp.viscosity_modification_type = 0.0

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
if num_type(size_parameter) == 'int':
    model_part_io_fluid = ModelPartIO(input_file_name.replace('error', 'error_ndiv_' + str(size_parameter)))
elif num_type(size_parameter) == 'float':
    model_part_io_fluid = ModelPartIO(input_file_name.replace('error', 'error_h_' + str(size_parameter)))

model_part_io_fluid.ReadModelPart(fluid_model_part)

#_____________________________________________________________________________________________________________________________________
#
#                               F L U I D    B L O C K    E N D S
#_____________________________________________________________________________________________________________________________________

# Constructing a utilities objects
creator_destructor = ParticleCreatorDestructor()
dem_fem_search = DEM_FEM_Search()

#Getting chosen scheme:
if DEM_parameters["IntegrationScheme"].GetString() == 'Forward_Euler':
    scheme = ForwardEulerScheme()
elif DEM_parameters["IntegrationScheme"].GetString() == 'Symplectic_Euler':
    if pp.CFD_DEM.basset_force_type > 0:
        scheme = SymplecticEulerOldVelocityScheme()
    else:
        scheme = SymplecticEulerScheme()
elif DEM_parameters["IntegrationScheme"].GetString() == 'Taylor_Scheme':
    scheme = TaylorScheme()
elif DEM_parameters["IntegrationScheme"].GetString() == 'Newmark_Beta_Method':
    scheme = NewmarkBetaScheme(0.5, 0.25)
elif DEM_parameters["IntegrationScheme"].GetString() == 'Verlet_Velocity':
    scheme = VerletVelocityScheme()
elif DEM_parameters["IntegrationScheme"].GetString() == 'Hybrid_Bashforth':
    scheme = HybridBashforthScheme()
else:
    KRATOSprint('Error: selected scheme not defined. Please select a different scheme')

if DEM_parameters.ElementType == "SwimmingNanoParticle":
    scheme = TerminalVelocityScheme()

# Creating a solver object and set the search strategy
solver = SolverStrategy.SwimmingStrategy(all_model_parts, creator_destructor, dem_fem_search, scheme, DEM_parameters, procedures)

# Add variables

procedures.AddAllVariablesInAllModelParts(solver, scheme, all_model_parts, DEM_parameters)
vars_man.AddNodalVariables(spheres_model_part, pp.dem_vars)
vars_man.AddNodalVariables(rigid_face_model_part, pp.rigid_faces_vars)
vars_man.AddNodalVariables(DEM_inlet_model_part, pp.inlet_vars)

# If Lalplacian form = 2, free all pressure Dofs
# laplacian_form = pp.laplacian_form
# if(laplacian_form >= 2):
    # for node in fluid_model_part.Nodes:
        # node.Free(PRESSURE)

# adding extra process info variables
vars_man.AddingExtraProcessInfoVariables(pp, fluid_model_part, spheres_model_part)

# defining a model part for the mixed part
mixed_model_part = ModelPart("MixedPart")


# Reading the model_part
spheres_mp_filename   = DEM_parameters.problem_name + "DEM"
model_part_io_spheres = ModelPartIO(spheres_mp_filename)

if (hasattr(DEM_parameters, "do_not_perform_initial_partition") and DEM_parameters.do_not_perform_initial_partition == 1):
    pass
else:
    parallelutils.PerformInitialPartition(model_part_io_spheres)

[model_part_io_spheres, spheres_model_part, MPICommSetup] = parallelutils.SetCommunicator(spheres_model_part, model_part_io_spheres, spheres_mp_filename)

model_part_io_spheres.ReadModelPart(spheres_model_part)

rigidFace_mp_filename = DEM_parameters.problem_name + "DEM_FEM_boundary"
model_part_io_fem = ModelPartIO(rigidFace_mp_filename)
model_part_io_fem.ReadModelPart(rigid_face_model_part)

clusters_mp_filename = DEM_parameters.problem_name + "DEM_Clusters"
model_part_io_clusters = ModelPartIO(clusters_mp_filename)
model_part_io_clusters.ReadModelPart(cluster_model_part)

DEM_Inlet_filename = DEM_parameters.problem_name + "DEM_Inlet"
model_part_io_demInlet = ModelPartIO(DEM_Inlet_filename)
model_part_io_demInlet.ReadModelPart(DEM_inlet_model_part)

# Setting up the buffer size
spheres_model_part.SetBufferSize(1)
cluster_model_part.SetBufferSize(1)
DEM_inlet_model_part.SetBufferSize(1)
rigid_face_model_part.SetBufferSize(1)
fluid_model_part.SetBufferSize(3)

# Adding dofs
solver.AddDofs(spheres_model_part)
solver.AddDofs(cluster_model_part)
solver.AddDofs(DEM_inlet_model_part)
solver_module.AddDofs(fluid_model_part, SolverSettings)
swim_proc.AddExtraDofs(pp, fluid_model_part, spheres_model_part, cluster_model_part, DEM_inlet_model_part)

# copy Y_WALL
for node in fluid_model_part.Nodes:
    y = node.GetSolutionStepValue(Y_WALL, 0)
    node.SetValue(Y_WALL, y)

# Creating necessary directories
main_path = os.getcwd()
[post_path, data_and_results, graphs_path, MPI_results] = procedures.CreateDirectories(str(main_path), str(DEM_parameters.problem_name), run_code)

os.chdir(main_path)

KRATOSprint("\nInitializing Problem...")

# Initialize GiD-IO
demio.AddGlobalVariables()
demio.AddSpheresVariables()
demio.AddFEMBoundaryVariables()
demio.AddClusterVariables()
demio.AddContactVariables()
# MPI
demio.AddMpiVariables()

demio.Configure(DEM_parameters.problem_name,
                DEM_parameters.OutputFileType,
                DEM_parameters.Multifile,
                DEM_parameters.ContactMeshOption)

demio.SetOutputName(DEM_parameters.problem_name)

os.chdir(post_path)

demio.InitializeMesh(all_model_parts)

os.chdir(post_path)

# Perform a partition to balance the problem
parallelutils.Repart(spheres_model_part)
parallelutils.CalculateModelNewIds(spheres_model_part)

os.chdir(post_path)

#Setting up the BoundingBox
bounding_box_time_limits = []
if (DEM_parameters.BoundingBoxOption == "ON"):
    procedures.SetBoundingBox(spheres_model_part, cluster_model_part, rigid_face_model_part, creator_destructor)
    bounding_box_time_limits = [solver.bounding_box_start_time, solver.bounding_box_stop_time]

#Creating a solver object and set the search strategy
#solver                 = SolverStrategy.ExplicitStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, DEM_parameters)
solver.search_strategy = parallelutils.GetSearchStrategy(solver, spheres_model_part)

# Creating the fluid solver
fluid_solver = solver_module.CreateSolver(fluid_model_part, SolverSettings)

Dt_DEM = DEM_parameters.MaxTimeStep

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

solver.CreateCPlusPlusStrategy()
solver.Initialize()    # Possible modifications of DELTA_TIME

fluid_solver.Initialize()
print("fluid solver created")
sys.stdout.flush()

if (DEM_parameters.ContactMeshOption =="ON"):
    contact_model_part = solver.contact_model_part

# constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation
# Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part

if not pp.VolumeOutput:
    cut_list = define_output.DefineCutPlanes()
    gid_io.define_cuts(fluid_model_part, cut_list)

# gid_io.initialize_results(fluid_model_part) # MOD.

import swimming_DEM_gid_output
swimming_DEM_gid_io = swimming_DEM_gid_output.SwimmingDEMGiDOutput(input_file_name,
                                                                   pp.VolumeOutput,
                                                                   pp.GiDPostMode,
                                                                   pp.GiDMultiFileFlag,
                                                                   pp.GiDWriteMeshFlag,
                                                                   pp.GiDWriteConditionsFlag)

swimming_DEM_gid_io.initialize_swimming_DEM_results(spheres_model_part, cluster_model_part, rigid_face_model_part, mixed_model_part)

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

# setting fluid's body force to the same as DEM's
if DEM_parameters.body_force_on_fluid_option:

    for node in fluid_model_part.Nodes:
        node.SetSolutionStepValue(BODY_FORCE_X, 0, DEM_parameters.GravityX)
        node.SetSolutionStepValue(BODY_FORCE_Y, 0, DEM_parameters.GravityY)
        node.SetSolutionStepValue(BODY_FORCE_Z, 0, DEM_parameters.GravityZ)

# coarse-graining: applying changes to the physical properties of the model to adjust for
# the similarity transformation if required (fluid effects only).
swim_proc.ApplySimilarityTransformations(fluid_model_part, DEM_parameters.similarity_transformation_type , DEM_parameters.model_over_real_diameter_factor )

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
DEM_parameters.n_particles_in_depth = int(math.sqrt(n_balls / fluid_volume)) # only relevant in 2D problems
# creating a physical calculations module to analyse the DEM model_part
dem_physics_calculator = SphericElementGlobalPhysicsCalculator(spheres_model_part)

if DEM_parameters["coupling_level_type"].GetInt():

    if DEM_parameters.meso_scale_length <= 0.0 and spheres_model_part.NumberOfElements(0) > 0:
        biggest_size = 2 * dem_physics_calculator.CalculateMaxNodalVariable(spheres_model_part, RADIUS)
        DEM_parameters.meso_scale_length  = 20 * biggest_size

    elif spheres_model_part.NumberOfElements(0) == 0:
        DEM_parameters.meso_scale_length  = 1.0

    #G
    L = 0.1
    U = 0.3
    k = 2.72
    omega = 0*math.pi

    flow_field = CellularFlowField(L, U, k, omega)
    space_time_set = SpaceTimeSet()
    field_utility = FluidFieldUtility(space_time_set, flow_field, 1000.0, 1e-6)
    #Z

    projection_module = CFD_DEM_coupling.ProjectionModule(fluid_model_part, spheres_model_part, rigid_face_model_part, domain_size, pp, field_utility)
    projection_module.UpdateDatabase(h_min)

# creating a custom functions calculator for the implementation of additional custom functions
custom_functions_tool = swim_proc.FunctionsCalculator(pp)

# creating a derivative recovery tool to calculate the necessary derivatives from the fluid solution (gradient, laplacian, material acceleration...)
derivative_recovery_tool = DerivativeRecoveryTool3D(fluid_model_part)

# creating a basset_force tool to perform the operations associated with the calculation of this force along the path of each particle
if pp.CFD_DEM.basset_force_type > 0:
    basset_force_tool = swim_proc.BassetForceTools()

# creating a stationarity assessment tool
stationarity_tool = swim_proc.StationarityAssessmentTool(DEM_parameters.max_pressure_variation_rate_tol , custom_functions_tool)

# creating a debug tool
dem_volume_tool = swim_proc.ProjectionDebugUtils(DEM_parameters.fluid_domain_volume, fluid_model_part, spheres_model_part, custom_functions_tool)

# creating a distance calculation process for the embedded technology
# (used to calculate elemental distances defining the structure embedded in the fluid mesh)
if (DEM_parameters.embedded_option):
    calculate_distance_process = CalculateSignedDistanceTo3DSkinProcess(rigid_face_model_part, fluid_model_part)
    calculate_distance_process.Execute()

# constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation

if (DEM_parameters.dem_inlet_option):
    max_DEM_node_Id = creator_destructor.FindMaxNodeIdInModelPart(spheres_model_part)
    max_elem_Id = creator_destructor.FindMaxElementIdInModelPart(spheres_model_part)
    max_FEM_node_Id = creator_destructor.FindMaxNodeIdInModelPart(rigid_face_model_part)
    max_fluid_node_Id = swim_proc.FindMaxNodeId(fluid_model_part)
    max_node_Id = max(max_DEM_node_Id, max_FEM_node_Id, max_fluid_node_Id, max_elem_Id)
    #max_node_Id = max(max_DEM_node_Id, max_FEM_node_Id, max_fluid_node_Id)

    creator_destructor.SetMaxNodeId(max_node_Id)

    # constructing the inlet and initializing it (must be done AFTER the spheres_model_part Initialize)
    DEM_inlet = DEM_Inlet(DEM_inlet_model_part)
    DEM_inlet.InitializeDEM_Inlet(spheres_model_part, creator_destructor, False)

DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(DEM_parameters, graphs_path, spheres_model_part, rigid_face_model_part)


KRATOSprint("Initialization Complete" + "\n")

step           = 0
time           = pp.Start_time
Dt             = pp.Dt
out            = Dt
Nsteps         = pp.nsteps
final_time     = pp.CFD_DEM.FinalTime
output_time    = pp.CFD_DEM.OutputTimeStep

report.Prepare(timer, DEM_parameters.ControlTime)

first_print = True; index_5 = 1; index_10 = 1; index_50 = 1; control = 0.0

if (DEM_parameters.ModelDataInfo == "ON"):
    os.chdir(data_and_results)
    if (DEM_parameters.ContactMeshOption == "ON"):
        (coordination_number) = procedures.ModelData(spheres_model_part, solver) # Calculates the mean number of neighbours the mean radius, etc..
        KRATOSprint ("Coordination Number: " + str(coordination_number) + "\n")
        os.chdir(main_path)
    else:
        KRATOSprint("Activate Contact Mesh for ModelData information")

if DEM_parameters.flow_in_porous_medium_option:
    fluid_frac_util = swim_proc.FluidFractionFieldUtility(fluid_model_part, DEM_parameters.min_fluid_fraction )

    for field in pp.fluid_fraction_fields:
        fluid_frac_util.AppendLinearField(field)

    fluid_frac_util.AddFluidFractionField()

if DEM_parameters["flow_in_porous_DEM_medium_option"].GetBool():
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

# setting up loop counters: Counter(steps_per_tick_step, initial_step, active_or_inactive_boolean)
embedded_counter             = swim_proc.Counter(1,
                                                 3,
                                                 DEM_parameters.embedded_option)  # MA: because I think DISTANCE,1 (from previous time step) is not calculated correctly for step=1
DEM_to_fluid_counter         = swim_proc.Counter(1,
                                                 1,
                                                 DEM_parameters["coupling_level_type"].GetInt() > 1)
derivative_recovery_counter    = swim_proc.Counter(1,
                                                 4,
                                                 DEM_parameters["coupling_level_type"].GetInt() or pp.CFD_DEM.print_PRESSURE_GRADIENT_option)
stationarity_counter         = swim_proc.Counter(DEM_parameters.time_steps_per_stationarity_step,
                                                 1,
                                                 DEM_parameters.stationary_problem_option)
print_counter                = swim_proc.Counter(1,
                                                 1,
                                                 out >= output_time)
debug_info_counter           = swim_proc.Counter(DEM_parameters.debug_tool_cycle,
                                                 1,
                                                 DEM_parameters.print_debug_info_option)
particles_results_counter    = swim_proc.Counter(DEM_parameters.print_particles_results_cycle,
                                                 1,
                                                 DEM_parameters.print_particles_results_option)
quadrature_counter           = swim_proc.Counter(pp.CFD_DEM.time_steps_per_quadrature_step,
                                                 1,
                                                 pp.CFD_DEM.print_BASSET_FORCE_option)
mat_deriv_averager           = swim_proc.Averager(1, 3)
laplacian_averager           = swim_proc.Averager(1, 3)
#G

# NANO BEGIN
frequence = 1

if DEM_parameters.ElementType == "SwimmingNanoParticle":
    frequence = pp.cation_concentration_frequence

cation_concentration_counter = swim_proc.Counter(frequence, 1, pp.CFD_DEM.drag_force_type == 9)
# NANO END

#fluid_model_part.AddNodalSolutionStepVariable(BODY_FORCE)

#for node in fluid_model_part.Nodes:
    #node.SetSolutionStepValue(BODY_FORCE_X, 0, 0.0)
    #node.SetSolutionStepValue(BODY_FORCE_Y, 0, 0.0)
    #node.SetSolutionStepValue(BODY_FORCE_Z, 0, - 9.81)
#import meshing_utils

#meshing_utils.ModifyRadiusesGivenPorosity(model_part, porosity, tol)
#Z


##############################################################################
#                                                                            #
#    MAIN LOOP                                                               #
#                                                                            #
##############################################################################

DEM_step     = 0      # necessary to get a good random insertion of particles   # relevant to the stationarity assessment tool
time_dem     = 0.0
Dt_DEM       = spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
rigid_face_model_part.ProcessInfo[DELTA_TIME] = Dt_DEM
cluster_model_part.ProcessInfo[DELTA_TIME] = Dt_DEM
stationarity = False

report.total_steps_expected = int(DEM_parameters.FinalTime / Dt_DEM)

KRATOSprint(report.BeginReport(timer))

mesh_motion = DEMFEMUtilities()

# creating a Post Utils object that executes several post-related tasks
post_utils_DEM = DEM_procedures.PostUtils(DEM_parameters, spheres_model_part)

swim_proc.InitializeVariablesWithNonZeroValues(fluid_model_part, spheres_model_part, pp) # otherwise variables are set to 0 by default


# ANALYTICS BEGIN
pp.CFD_DEM.perform_analytics_option = False

if pp.CFD_DEM.perform_analytics_option:
    import analytics
    variables_to_measure = [PRESSURE]
    steps_between_measurements = 100
    gauge = analytics.Gauge(fluid_model_part, Dt, final_time, variables_to_measure, steps_between_measurements)
    point_coors = [0.0, 0.0, 0.01]
    target_node = swim_proc.FindClosestNode(fluid_model_part, point_coors)
    target_id = target_node.Id
    print(target_node.X, target_node.Y, target_node.Z)
    print(target_id)
    def condition(node):
        return node.Id == target_id

    gauge.ConstructArrayOfNodes(condition)
    print(gauge.variables)
    print_analytics_counter = swim_proc.Counter( 5 * steps_between_measurements, 1, 1)
# ANALYTICS END

# NANO BEGIN
if pp.CFD_DEM.drag_force_type == 9:
    alpha = 1000.0

    for node in spheres_model_part.Nodes:
        node.SetSolutionStepValue(CATION_CONCENTRATION, pp.initial_concentration)

# NANO END

#G

import derivative_recovery.derivative_recovery_strategy as derivative_recoverer
recovery = derivative_recoverer.DerivativeRecoveryStrategy(pp, fluid_model_part, derivative_recovery_tool, custom_functions_tool)

number=0
for node in fluid_model_part.Nodes:
    number += 1
#Z

N_steps = int(final_time / Dt_DEM) + 20

if pp.CFD_DEM.basset_force_type > 0:
    basset_force_tool.FillDaitcheVectors(N_steps, pp.CFD_DEM.quadrature_order, pp.CFD_DEM.time_steps_per_quadrature_step)
if pp.CFD_DEM.basset_force_type >= 3 or pp.CFD_DEM.basset_force_type == 1:
    basset_force_tool.FillHinsbergVectors(spheres_model_part, pp.CFD_DEM.number_of_exponentials, pp.CFD_DEM.number_of_quadrature_steps_in_window)

post_utils.Writeresults(time)

mat_deriv_errors = []
laplacian_errors = []
current_mat_deriv_errors = np.zeros(2)
current_laplacian_errors = np.zeros(2)
#fluid_model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 1)

while (time <= final_time):

    time = time + Dt
    step += 1
    fluid_model_part.CloneTimeStep(time)
    print("\n", "TIME = ", time)
    print('ELAPSED TIME = ', timer.time() - init_time)
    sys.stdout.flush()

    if DEM_parameters["coupling_scheme_type"].GetString()  == "UpdatedDEM":
        time_final_DEM_substepping = time + Dt

    else:
        time_final_DEM_substepping = time

    # calculating elemental distances defining the structure embedded in the fluid mesh
    if DEM_parameters.embedded_option:
        calculate_distance_process.Execute()

    if embedded_counter.Tick():
        embedded.ApplyEmbeddedBCsToFluid(fluid_model_part)
        embedded.ApplyEmbeddedBCsToBalls(spheres_model_part, DEM_parameters)

    # solving the fluid part

    if step >= 3 and not stationarity:
        print("Solving Fluid... (", fluid_model_part.NumberOfElements(0), "elements )")
        sys.stdout.flush()

        if not pp.CFD_DEM.drag_force_type == 9:
            #fluid_solver.Solve()
#G
            for node in fluid_model_part.Nodes:
                vel= Vector(3)
                coor= Vector(3)
                coor[0]=node.X
                coor[1]=node.Y
                coor[2]=node.Z
                flow_field.Evaluate(time,coor,vel,0)
                node.SetSolutionStepValue(VELOCITY_X, vel[0])
                node.SetSolutionStepValue(VELOCITY_Y, vel[1])
                node.SetSolutionStepValue(VELOCITY_Z, vel[2])
                #node.SetSolutionStepValue(VELOCITY_X, 7*node.X**2 + 6 * node.Y + 19)
                #node.SetSolutionStepValue(VELOCITY_Y,  9 *node.X**2 - 8 * node.Y**2 +30*node.X)
                #node.SetSolutionStepValue(VELOCITY_Z, 0.0)
#Z


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

    if output_time <= out and DEM_parameters["coupling_scheme_type"].GetString() == "UpdatedDEM":

        if DEM_parameters["coupling_level_type"].GetInt() > 0:
            projection_module.ComputePostProcessResults(spheres_model_part.ProcessInfo)

        post_utils.Writeresults(time)
        out = 0

    # solving the DEM part
    derivative_recovery_counter.Deactivate(time < DEM_parameters.interaction_start_time)

    if derivative_recovery_counter.Tick():
        recovery.Recover()

#G
    error_mat_deriv = 0.
    max_error_mat_deriv = - float('inf')
    error_laplacian = 0.
    max_error_laplacian = - float('inf')

    total_volume = 0.
    mat_deriv_average = Vector(3)
    laplacian_average = Vector(3)
    for k in range(3):
        mat_deriv_average[k] = 0.0
        laplacian_average[k] = 0.0
    norm_mat_deriv_average = 0.
    norm_laplacian_average = 0.

    if step > 3:
        module_mat_deriv = 0.
        module_laplacian = 0.
        for i_node, node in enumerate(fluid_model_part.Nodes):
            calc_mat_deriv = [0.] * 3
            calc_laplacian = [0.] * 3
            mat_deriv= Vector(3)
            laplacian= Vector(3)
            coor= Vector(3)
            coor[0]=node.X
            coor[1]=node.Y
            coor[2]=0
            flow_field.CalculateMaterialAcceleration(time, coor, mat_deriv, 0)
            flow_field.CalculateLaplacian(time, coor, laplacian, 0)
            calc_mat_deriv[0] = node.GetSolutionStepValue(MATERIAL_ACCELERATION_X)
            calc_mat_deriv[1] = node.GetSolutionStepValue(MATERIAL_ACCELERATION_Y)
            calc_mat_deriv[2] = node.GetSolutionStepValue(MATERIAL_ACCELERATION_Z)
            calc_laplacian[0] = node.GetSolutionStepValue(VELOCITY_LAPLACIAN_X)
            calc_laplacian[1] = node.GetSolutionStepValue(VELOCITY_LAPLACIAN_Y)
            calc_laplacian[2] = node.GetSolutionStepValue(VELOCITY_LAPLACIAN_Z)
            module_mat_deriv += math.sqrt(calc_mat_deriv[0] ** 2 + calc_mat_deriv[1] ** 2 + calc_mat_deriv[2] ** 2)
            module_laplacian += math.sqrt(calc_laplacian[0] ** 2 + calc_laplacian[1] ** 2 + calc_laplacian[2] ** 2)
            #module_mat_deriv = max(math.sqrt(mat_deriv[0] ** 2 + mat_deriv[1] ** 2 + mat_deriv[2] ** 2), 1e-8)
            #module_laplacian = max(math.sqrt(laplacian[0] ** 2 + laplacian[1] ** 2 + laplacian[2] ** 2), 1e-8)
            #laplacian[0] = 14
            #laplacian[1] = 2
            #laplacian[2] = 0
            #nodal_volume = node.GetSolutionStepValue(NODAL_AREA)
            #total_volume += nodal_volume
            current_error = swim_proc.NormOfDifference(calc_mat_deriv, mat_deriv)
            error_mat_deriv += current_error
            max_error_mat_deriv = max(max_error_mat_deriv, current_error)
            current_error = swim_proc.NormOfDifference(calc_laplacian, laplacian)
            error_laplacian += current_error
            max_error_laplacian = max(max_error_laplacian, current_error)
            diff_mat_deriv = [calc_mat_deriv[i] - mat_deriv[i] for i in range(len(calc_mat_deriv))]
            diff_laplacian = [calc_laplacian[i] - laplacian[i] for i in range(len(calc_laplacian))]
            #mat_deriv_averager.swim_proc.Norm(diff_mat_deriv)
            #laplacian_averager.swim_proc.Norm(diff_laplacian)
            #for k in range(3):
                #mat_deriv_average[k] += mat_deriv[k]
                #laplacian_average[k] += laplacian[k]
            norm_mat_deriv_average += swim_proc.Norm(mat_deriv)
            norm_laplacian_average += swim_proc.Norm(laplacian)

            node.SetSolutionStepValue(MATERIAL_ACCELERATION_X, calc_mat_deriv[0] - mat_deriv[0])
            node.SetSolutionStepValue(MATERIAL_ACCELERATION_Y, calc_mat_deriv[1] - mat_deriv[1])
            node.SetSolutionStepValue(MATERIAL_ACCELERATION_Z, calc_mat_deriv[2] - mat_deriv[2])
            #node.SetSolutionStepValue(MATERIAL_ACCELERATION_X, mat_deriv[0])
            #node.SetSolutionStepValue(MATERIAL_ACCELERATION_Y, mat_deriv[1])
            #node.SetSolutionStepValue(MATERIAL_ACCELERATION_Z, mat_deriv[2])


            node.SetSolutionStepValue(VELOCITY_LAPLACIAN_X, calc_laplacian[0] - laplacian[0])
            node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Y, calc_laplacian[1] - laplacian[1])
            node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Z, calc_laplacian[2] - laplacian[2])
            #node.SetSolutionStepValue(VELOCITY_LAPLACIAN_X, calc_laplacian_0)
            #node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Y, calc_laplacian_1)
            #node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Z, calc_laplacian_2)
            #node.SetSolutionStepValue(VELOCITY_LAPLACIAN_X, laplacian[0])
            #node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Y, laplacian[1])
            #node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Z, laplacian[2])

        module_mat_deriv /= len(fluid_model_part.Nodes)
        module_laplacian /= len(fluid_model_part.Nodes)
        swim_proc.MultiplyNodalVariableByFactor(fluid_model_part, MATERIAL_ACCELERATION, 1.0 / module_mat_deriv)
        swim_proc.MultiplyNodalVariableByFactor(fluid_model_part, VELOCITY_LAPLACIAN, 1.0 / module_laplacian)

        if norm_mat_deriv_average > 0. and norm_laplacian_average > 0:
            current_mat_deriv_errors[0] = error_mat_deriv / norm_mat_deriv_average
            current_mat_deriv_errors[1] = max_error_mat_deriv / norm_mat_deriv_average * len(fluid_model_part.Nodes)
            current_laplacian_errors[0] = error_laplacian / norm_laplacian_average
            current_laplacian_errors[1] = max_error_laplacian / norm_laplacian_average * len(fluid_model_part.Nodes)
            mat_deriv_errors.append(current_mat_deriv_errors)
            laplacian_errors.append(current_laplacian_errors)
            #print('mat_deriv: min, max, avg, ', mat_deriv_averager.GetCurrentData())
            #print('laplacian: min, max, avg, ', laplacian_averager.GetCurrentData())
            print('rel_error_mat_deriv', error_mat_deriv / norm_mat_deriv_average)
            print('rel_error_laplacian', error_laplacian / norm_laplacian_average)
#Z
    print("Solving DEM... (", spheres_model_part.NumberOfElements(0), "elements )")
    sys.stdout.flush()
    first_dem_iter = True

    for time_dem in yield_DEM_time(time_dem, time_final_DEM_substepping, Dt_DEM):
        DEM_step += 1   # this variable is necessary to get a good random insertion of particles

        spheres_model_part.ProcessInfo[TIME_STEPS]    = DEM_step
        rigid_face_model_part.ProcessInfo[TIME_STEPS] = DEM_step
        cluster_model_part.ProcessInfo[TIME_STEPS]    = DEM_step
        # NANO BEGIN
        if cation_concentration_counter.Tick():
            concentration = pp.final_concentration + (pp.initial_concentration - pp.final_concentration) * math.exp(- alpha * time_dem)

            for node in spheres_model_part.Nodes:
                node.SetSolutionStepValue(CATION_CONCENTRATION, concentration)
        # NANO END
        # applying fluid-to-DEM coupling if required

        if time >= DEM_parameters.interaction_start_time and DEM_parameters["coupling_level_type"].GetInt() and (DEM_parameters.project_at_every_substep_option or first_dem_iter):

            if DEM_parameters["coupling_scheme_type"].GetString() == "UpdatedDEM":
                projection_module.ProjectFromNewestFluid()

            else:
                projection_module.ProjectFromFluid((time_final_DEM_substepping - time_dem) / Dt)

                if DEM_parameters["IntegrationScheme"].GetString() == 'Hybrid_Bashforth':
                    solver.Solve() # only advance in space
                    projection_module.InterpolateVelocity()

                if quadrature_counter.Tick():
                    if pp.CFD_DEM.basset_force_type == 1 or pp.CFD_DEM.basset_force_type >= 3:
                        basset_force_tool.AppendIntegrandsWindow(spheres_model_part)
                    elif pp.CFD_DEM.basset_force_type == 2:
                        basset_force_tool.AppendIntegrands(spheres_model_part)

        # performing the time integration of the DEM part

        spheres_model_part.ProcessInfo[TIME]    = time_dem
        rigid_face_model_part.ProcessInfo[TIME] = time_dem
        cluster_model_part.ProcessInfo[TIME]    = time_dem

        if not DEM_parameters["flow_in_porous_DEM_medium_option"].GetBool(): # in porous flow particles remain static
            #solver.Solve()
            pass

        # Walls movement:
        mesh_motion.MoveAllMeshes(rigid_face_model_part, time, Dt)
        mesh_motion.MoveAllMeshes(spheres_model_part, time, Dt)
        mesh_motion.MoveAllMeshes(DEM_inlet_model_part, time, Dt)

        #### TIME CONTROL ##################################

        # adding DEM elements by the inlet:
        if (DEM_parameters.dem_inlet_option):
            DEM_inlet.CreateElementsFromInletMesh(spheres_model_part, cluster_model_part, creator_destructor)  # After solving, to make sure that neighbours are already set.

        #first_dem_iter = False

    if DEM_parameters.ElementType == "SwimmingNanoParticle":
        print("concentration: " + str(concentration))

    #### PRINTING GRAPHS ####
    os.chdir(graphs_path)
    # measuring mean velocities in a certain control volume (the 'velocity trap')
    if (DEM_parameters.VelocityTrapOption):
        post_utils_DEM.ComputeMeanVelocitiesinTrap("Average_Velocity.txt", time)

    os.chdir(post_path)

    # applying DEM-to-fluid coupling

    if DEM_to_fluid_counter.Tick() and time >= DEM_parameters.interaction_start_time:
        projection_module.ProjectFromParticles()

    # coupling checks (debugging)

    if debug_info_counter.Tick():
        dem_volume_tool.UpdateDataAndPrint(DEM_parameters.fluid_domain_volume)

    # printing if required

    if particles_results_counter.Tick():

        # eliminating remote balls

        #if pp.dem.BoundingBoxOption == "ON":
            #creator_destructor.DestroyParticlesOutsideBoundingBox(spheres_model_part)

        io_tools.PrintParticlesResults(pp.variables_to_print_in_file, time, spheres_model_part)
        graph_printer.PrintGraphs(time)
        PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

    if output_time <= out and DEM_parameters["coupling_scheme_type"].GetString() == "UpdatedFluid":

        if DEM_parameters["coupling_level_type"].GetInt():
            projection_module.ComputePostProcessResults(spheres_model_part.ProcessInfo)

        post_utils.Writeresults(time)
        out = 0

    out = out + Dt

swimming_DEM_gid_io.finalize_results()
print("\n CALCULATIONS FINISHED. THE SIMULATION ENDED SUCCESSFULLY.")
simulation_elapsed_time = timer.clock() - simulation_start_time
print("Elapsed time: " + "%.5f"%(simulation_elapsed_time) + " s ")
print("per fluid time step: " + "%.5f"%(simulation_elapsed_time/ step) + " s ")
print("per DEM time step: " + "%.5f"%(simulation_elapsed_time/ DEM_step) + " s")

if not os.path.exists('../errors_recorded'):
    os.makedirs('../errors_recorded')

import h5py

file_name = main_path + '/errors_recorded/recovery_errors.hdf5'
# with h5py.File(self.file_name, 'r+') as f:
#     f.create_dataset('material_derivative', shape = self.shape, dtype = np.float32)

with h5py.File(file_name) as f:
    mat_deriv_grp = f.require_group('material derivative')
    mat_deriv_mthd_group = mat_deriv_grp.require_group('method = ' + str(pp.CFD_DEM.material_acceleration_calculation_type))
    laplacian_grp = f.require_group('laplacian')
    laplacian_mthd_group = laplacian_grp.require_group('method = ' + str(pp.CFD_DEM.laplacian_calculation_type))

    if num_type(size_parameter) == 'int':
        mesh_grp_mat_deriv = mat_deriv_mthd_group.require_group('regular mesh')
        mesh_grp_laplacian = laplacian_mthd_group.require_group('regular mesh')
        dset_mat_deriv = mesh_grp_mat_deriv.require_dataset('n_div = ' + str(size_parameter), (2,), dtype = np.float64)
        dset_laplacian = mesh_grp_laplacian.require_dataset('n_div = ' + str(size_parameter), (2,), dtype = np.float64)
    elif num_type(size_parameter) == 'float':
        mesh_grp_mat_deriv = mat_deriv_mthd_group.require_group('irregular mesh')
        mesh_grp_laplacian = laplacian_mthd_group.require_group('irregular mesh')
        dset_mat_deriv = mesh_grp_mat_deriv.require_dataset('h = '   + str(size_parameter), (2,), dtype = np.float64)
        dset_laplacian = mesh_grp_laplacian.require_dataset('h = '   + str(size_parameter), (2,), dtype = np.float64)
    dset_mat_deriv.attrs['mesh size'] =  min(float(size_parameter), 1.0 / float(size_parameter))
    dset_mat_deriv[0] = current_mat_deriv_errors[0]
    dset_mat_deriv[1] = current_mat_deriv_errors[1]
    dset_laplacian[0] = current_laplacian_errors[0]
    dset_laplacian[1] = current_laplacian_errors[1]

if num_type(size_parameter) == 'int':
    size_parameter_name = '_ndiv_'
elif num_type(size_parameter) == 'float':
    size_parameter_name = '_h_'
with open('../errors_recorded/mat_deriv_errors' + size_parameter_name + str(size_parameter) + '_type_' + str(pp.CFD_DEM.material_acceleration_calculation_type) + '.txt', 'w') as mat_errors_file:
    for error in mat_deriv_errors:
        mat_errors_file.write(str(error) + '\n')
with open('../errors_recorded/laplacian_errors' + size_parameter_name + str(size_parameter) + '_type_' + str(pp.CFD_DEM.laplacian_calculation_type) + '.txt', 'w') as laplacian_errors_file:
    for error in laplacian_errors:
        laplacian_errors_file.write(str(error) + '\n')
sys.stdout.flush()
os.chdir(main_path)
sys.stdout.path_to_console_out_file
os.rename(sys.stdout.path_to_console_out_file, post_path + '/' + sys.stdout.path_to_console_out_file)
empty = False
add_to_name = ''
count = 0
dir_name = post_path + '_FINISHED_AT_t=' + str(round(time, 1))
#dir_name_it = dir_name
#while os.path.isdir(dir_name_it):
    #dir_name_it = dir_name + '_version_' + str(count + 1)
    #count += 1
if os.path.isdir(dir_name):
    import shutil
    shutil.rmtree(dir_name)
    print(dir_name)
    print(post_path)
os.rename(post_path, dir_name)

for i in drag_file_output_list:
    i.close()
