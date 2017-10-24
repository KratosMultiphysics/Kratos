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
import analytics

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

pp.CFD_DEM = DEM_parameters

##############################################################################
#                                                                            #
#    INITIALIZE                                                              #
#                                                                            #
##############################################################################

#G
file_name, j, k, current_run_path = sys.argv
print(j, k)
M = int(j)
TIME_WINDOW = 5 * pp.CFD_DEM.MaxTimeStep * 2 ** int(k)
sys.path.append(current_run_path)

pp.CFD_DEM.recover_gradient_option = True
pp.CFD_DEM.do_search_neighbours = False
#pp.CFD_DEM.print_PRESSURE_GRADIENT_option = False
DEM_parameters.fluid_domain_volume = 0.04 * math.pi # write down the volume you know it has
pp.CFD_DEM.faxen_terms_type = 0
pp.CFD_DEM.material_acceleration_calculation_type = 0
pp.CFD_DEM.faxen_force_type = 0
pp.CFD_DEM.print_FLUID_VEL_PROJECTED_RATE_option = 0
pp.CFD_DEM.basset_force_type = 4
pp.CFD_DEM.print_BASSET_FORCE_option = 1
pp.CFD_DEM.basset_force_integration_type = 1
pp.CFD_DEM.n_init_basset_steps = 2
pp.CFD_DEM.time_steps_per_quadrature_step = 1
pp.CFD_DEM.delta_time_quadrature = pp.CFD_DEM.time_steps_per_quadrature_step * pp.CFD_DEM.MaxTimeStep
pp.CFD_DEM.quadrature_order = 2
pp.CFD_DEM.time_window = TIME_WINDOW
pp.CFD_DEM.number_of_exponentials = M
pp.CFD_DEM.number_of_quadrature_steps_in_window = int(pp.CFD_DEM.time_window / pp.CFD_DEM.delta_time_quadrature)
pp.CFD_DEM.print_steps_per_plot_step = 1
pp.CFD_DEM.PostCationConcentration = False
pp.CFD_DEM.do_impose_flow_from_field = True
pp.CFD_DEM.print_MATERIAL_ACCELERATION_option = True
pp.CFD_DEM.print_FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED_option = False
number_of_vectors_to_be_kept_in_memory = pp.CFD_DEM.time_window / pp.CFD_DEM.MaxTimeStep * pp.CFD_DEM.time_steps_per_quadrature_step + pp.CFD_DEM.number_of_exponentials 
print('\nNumber of vectors to be kept in memory: ', number_of_vectors_to_be_kept_in_memory)
# Making the fluid step an exact multiple of the DEM step
pp.Dt = int(pp.Dt / pp.CFD_DEM.MaxTimeStep) * pp.CFD_DEM.MaxTimeStep

# Creating a code for the used input variables
run_code = swim_proc.CreateRunCode(pp)
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
rigid_face_model_part = ModelPart("RigidFace_Part")
cluster_model_part    = ModelPart("Cluster_Part")
DEM_inlet_model_part  = ModelPart("DEMInletPart")
mapping_model_part    = ModelPart("Mappingmodel_part")

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
model_part_io_fluid = ModelPartIO(input_file_name)
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
solver = SolverStrategy.SwimmingStrategy(spheres_model_part, rigid_face_model_part, cluster_model_part, DEM_inlet_model_part, creator_destructor, dem_fem_search, scheme, DEM_parameters, procedures)

# Add variables
procedures.AddCommonVariables(spheres_model_part, DEM_parameters)
procedures.AddSpheresVariables(spheres_model_part, DEM_parameters)
procedures.AddMpiVariables(spheres_model_part)
solver.AddAdditionalVariables(spheres_model_part, DEM_parameters)
scheme.AddSpheresVariables(spheres_model_part)
procedures.AddCommonVariables(cluster_model_part, DEM_parameters)
procedures.AddClusterVariables(cluster_model_part, DEM_parameters)
procedures.AddMpiVariables(cluster_model_part)
scheme.AddClustersVariables(cluster_model_part)
procedures.AddCommonVariables(DEM_inlet_model_part, DEM_parameters)
procedures.AddSpheresVariables(DEM_inlet_model_part, DEM_parameters)
solver.AddAdditionalVariables(DEM_inlet_model_part, DEM_parameters)  
scheme.AddSpheresVariables(DEM_inlet_model_part)
procedures.AddCommonVariables(rigid_face_model_part, DEM_parameters)
procedures.AddRigidFaceVariables(rigid_face_model_part, DEM_parameters)
procedures.AddMpiVariables(rigid_face_model_part)
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

demio.InitializeMesh(spheres_model_part,
                     rigid_face_model_part,
                     cluster_model_part,
                     solver.contact_model_part,
                     mapping_model_part)

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

    projection_module = CFD_DEM_coupling.ProjectionModule(fluid_model_part, spheres_model_part, rigid_face_model_part, domain_size, pp)
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
pressure_gradient_counter    = swim_proc.Counter(1, 
                                                 1, 
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
#G

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

swim_proc.InitializeVariablesWithNonZeroValues(fluid_model_part, spheres_model_part, pp) # all variables are set to 0 by default

#G
linear_solver = CGSolver()
scheme = ResidualBasedIncrementalUpdateStaticScheme()
post_process_strategy = ResidualBasedLinearStrategy(fluid_model_part, scheme, linear_solver, False, True, False, False)
#Z

# CHANDELIER BEGIN
import math 
import cmath
import mpmath
import matplotlib.pyplot as plt

import chandelier_parameters as ch_pp
import chandelier as ch
import quadrature as quad
sim = ch.AnalyticSimulator(ch_pp)

post_utils.Writeresults(time)
coors = [None] * 3
exact_vel = [None] * 3 
Dt_DEM_inv = 1.0 / Dt_DEM                    
vel = [0., 0.0, 0.] 
old_vel = [v for v in vel] 
H = [0.] * 3
H_old = [0.] * 3
Delta_H = [0.] * 3 
exact_Delta_H = [0.] * 3
basset_force = [0.] * 3
exact_basset_force = [0.] * 3 
sim.CalculatePosition(coors, 0.0)      
times = [0.0]
radii = [1.0]
integrands = []
particle_mass = 4. / 3 * math.pi * ch_pp.a ** 3 * ch_pp.rho_p
units_coefficient = 6 * ch_pp.a ** 2 * ch_pp.rho_f * math.sqrt(math.pi * ch_pp.nu)                        

# Impose initial velocity to be the terminal velocity
sim.CalculateNonDimensionalVars()
terminal_velocity = sim.NDw0 * ch_pp.R * ch_pp.omega

# NODE HISTORY RESULTS BEGIN
scalar_vars = []
vector_vars = [DISPLACEMENT]

for node in spheres_model_part.Nodes:
    node_to_follow_id = node.Id
    node.SetSolutionStepValue(VELOCITY_Z, terminal_velocity)

results_creator = swim_proc.ResultsFileCreator(spheres_model_part, node_to_follow_id, scalar_vars, vector_vars)
# NODE HISTORY RESULTS END 
# CHANDELLIER END
N_steps = int(final_time / Dt_DEM) + 20

if pp.CFD_DEM.basset_force_type > 0:
    basset_force_tool.FillDaitcheVectors(N_steps, pp.CFD_DEM.quadrature_order, pp.CFD_DEM.time_steps_per_quadrature_step)
if pp.CFD_DEM.basset_force_type >= 3 or pp.CFD_DEM.basset_force_type == 1:
    basset_force_tool.FillHinsbergVectors(spheres_model_part, pp.CFD_DEM.number_of_exponentials, pp.CFD_DEM.number_of_quadrature_steps_in_window)

for node in spheres_model_part.Nodes:
    node.SetSolutionStepValue(VELOCITY_Y, ch_pp.u0)    
    node.SetSolutionStepValue(VELOCITY_Y, ch_pp.v0)
    node.SetSolutionStepValue(VELOCITY_Z, 2. / 9 * 9.8 * ch_pp.a ** 2 / (ch_pp.nu * ch_pp.rho_f) * (ch_pp.rho_f - ch_pp.rho_p))
    node.Fix(VELOCITY_Z)
    node.SetSolutionStepValue(VELOCITY_OLD_X, ch_pp.u0)
    node.SetSolutionStepValue(VELOCITY_OLD_Y, ch_pp.v0)
    node.SetSolutionStepValue(VELOCITY_OLD_Z, 2. / 9 * 9.8 * ch_pp.a ** 2 / (ch_pp.nu * ch_pp.rho_f) * (ch_pp.rho_f - ch_pp.rho_p))
    node.Fix(VELOCITY_OLD_Z)
    node.SetSolutionStepValue(FLUID_VEL_PROJECTED_X, ch_pp.u0)
    node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Y, ch_pp.v0)
    node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Z, 0.0)   
r = 0
x = 0
y = 0

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
            if VELOCITY_LAPLACIAN in pp.fluid_vars:
                print("\nSolving for the Laplacian...")
                sys.stdout.flush()
                fractional_step = fluid_model_part.ProcessInfo[FRACTIONAL_STEP]
                fluid_model_part.ProcessInfo[FRACTIONAL_STEP] = 2

                #post_process_strategy.Solve()

                print("\nFinished solving for the Laplacian.")
                sys.stdout.flush()
                fluid_model_part.ProcessInfo[FRACTIONAL_STEP] = fractional_step
                
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
    pressure_gradient_counter.Deactivate(time < DEM_parameters.interaction_start_time)

    if pressure_gradient_counter.Tick():
        if pp.CFD_DEM.gradient_calculation_type == 2:
            derivative_recovery_tool.RecoverSuperconvergentGradient(fluid_model_part, PRESSURE, PRESSURE_GRADIENT)
        elif pp.CFD_DEM.gradient_calculation_type == 1:            
            custom_functions_tool.CalculatePressureGradient(fluid_model_part)            
        if pp.CFD_DEM.laplacian_calculation_type == 1:
            derivative_recovery_tool.CalculateVectorLaplacian(fluid_model_part, VELOCITY, VELOCITY_LAPLACIAN)
        elif pp.CFD_DEM.laplacian_calculation_type == 2:
            derivative_recovery_tool.RecoverSuperconvergentLaplacian(fluid_model_part, VELOCITY, VELOCITY_LAPLACIAN)
        if pp.CFD_DEM.material_acceleration_calculation_type == 1:
            derivative_recovery_tool.CalculateVectorMaterialDerivative(fluid_model_part, VELOCITY, ACCELERATION, MATERIAL_ACCELERATION)    

    print("Solving DEM... (", spheres_model_part.NumberOfElements(0), "elements )")
    sys.stdout.flush()
    first_dem_iter = True

    for time_dem in yield_DEM_time(time_dem, time_final_DEM_substepping, Dt_DEM):
        DEM_step += 1   # this variable is necessary to get a good random insertion of particles

        spheres_model_part.ProcessInfo[TIME_STEPS]    = DEM_step
        rigid_face_model_part.ProcessInfo[TIME_STEPS] = DEM_step
        cluster_model_part.ProcessInfo[TIME_STEPS]    = DEM_step

        # applying fluid-to-DEM coupling if required

        if time >= DEM_parameters.interaction_start_time and DEM_parameters["coupling_level_type"].GetInt() and (DEM_parameters.project_at_every_substep_option or first_dem_iter):
            if DEM_parameters["coupling_scheme_type"].GetString() == "UpdatedDEM":
                projection_module.ProjectFromNewestFluid()

            else:         
                projection_module.ProjectFromFluid((time_final_DEM_substepping - time_dem) / Dt)             
                
                for node in spheres_model_part.Nodes:
                    x = node.X
                    y = node.Y
                    z = node.Z
                    r = math.sqrt(x ** 2 + y ** 2)
                    omega = ch_pp.omega
                    vx = - omega * y
                    vy =   omega * x
                    ax = - x * omega ** 2
                    ay = - y * omega ** 2             
                    node.SetSolutionStepValue(FLUID_VEL_PROJECTED_X, vx)
                    node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Y, vy)
                    node.SetSolutionStepValue(FLUID_VEL_PROJECTED_Z, 0.0)            
                    node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED_X, ax)
                    node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED_Y, ay)
                    node.SetSolutionStepValue(FLUID_ACCEL_PROJECTED_Z, 0.0)      
                    node.SetSolutionStepValue(MATERIAL_FLUID_ACCEL_PROJECTED_X, ax)
                    node.SetSolutionStepValue(MATERIAL_FLUID_ACCEL_PROJECTED_Y, ay)
                    node.SetSolutionStepValue(MATERIAL_FLUID_ACCEL_PROJECTED_Z, 0.0)       

                    if DEM_parameters["IntegrationScheme"].GetString() == 'Hybrid_Bashforth':
                        solver.Solve() # only advance in space
                        #projection_module.InterpolateVelocity()  
                        x = node.X
                        y = node.Y
                        r = math.sqrt(x ** 2 + y ** 2)
                        new_vx = - omega * y
                        new_vy =   omega * x
                        node.SetSolutionStepValue(SLIP_VELOCITY_X, new_vx)
                        node.SetSolutionStepValue(SLIP_VELOCITY_Y, new_vy)       
                    else:
                        if pp.CFD_DEM.basset_force_type > 0:
                            node.SetSolutionStepValue(SLIP_VELOCITY_X, vx)
                            node.SetSolutionStepValue(SLIP_VELOCITY_Y, vy) 
                          
                    
                    vp_x = node.GetSolutionStepValue(VELOCITY_X)
                    vp_y = node.GetSolutionStepValue(VELOCITY_Y)
                    vp_z = node.GetSolutionStepValue(VELOCITY_Z) 
                    integrands.append([vx - vp_x, vy - vp_y, 0.])                                                                          

                    if quadrature_counter.Tick():
                        if pp.CFD_DEM.basset_force_type == 1 or pp.CFD_DEM.basset_force_type >= 3:
                            basset_force_tool.AppendIntegrandsWindow(spheres_model_part) 
                        elif pp.CFD_DEM.basset_force_type == 2:
                            basset_force_tool.AppendIntegrands(spheres_model_part)     

                    H_old[:] = H[:]
                    if len(times) < pp.CFD_DEM.n_init_basset_steps:                        
                        sim.CalculateBassetForce(basset_force, time_dem * ch_pp.omega)
                        node.SetSolutionStepValue(BASSET_FORCE_X, basset_force[0])
                        node.SetSolutionStepValue(BASSET_FORCE_Y, basset_force[1])    
                        node.SetSolutionStepValue(BASSET_FORCE_Z, 0.)  
                        sim.CalculatePosition(coors, time_dem * ch_pp.omega, exact_vel)    
                        #print("\nEXACT", basset_force) 
                        H[0] = H_old[0] + Dt_DEM / units_coefficient * basset_force[0]
                        H[1] = H_old[1] + Dt_DEM / units_coefficient * basset_force[1]

                        for node in spheres_model_part.Nodes:
                            node.X = coors[0] * ch_pp.R
                            node.Y = coors[1] * ch_pp.R
                            node.Z = coors[2] * ch_pp.R
                            x = node.X
                            y = node.Y
                            r = math.sqrt(x ** 2 + y ** 2)

                            node.SetSolutionStepValue(DISPLACEMENT_X, coors[0] * ch_pp.R - ch_pp.x0)
                            node.SetSolutionStepValue(DISPLACEMENT_Y, coors[1] * ch_pp.R - ch_pp.y0)
                            node.SetSolutionStepValue(DISPLACEMENT_Z, coors[2] * ch_pp.R - ch_pp.z0)
                            vp_x = exact_vel[0] * ch_pp.R * ch_pp.omega
                            vp_y = exact_vel[1] * ch_pp.R * ch_pp.omega
                            vp_z = exact_vel[2] * ch_pp.R * ch_pp.omega
                            #node.SetSolutionStepValue(VELOCITY_X, vp_x)
                            #node.SetSolutionStepValue(VELOCITY_Y, vp_y)
                            #node.SetSolutionStepValue(VELOCITY_Z, vp_z)
                            #node.SetSolutionStepValue(VELOCITY_OLD_X, vp_x)
                            #node.SetSolutionStepValue(VELOCITY_OLD_Y, vp_y)
                            #node.SetSolutionStepValue(VELOCITY_Z, vp_z)                            
                    
                    else:
                        pass
                        #sim.CalculateBassetForce(exact_basset_force, time_dem * ch_pp.omega)
                        #sim.CalculatePosition(coors, time_dem * ch_pp.omega, exact_vel)
                        #vp_x = exact_vel[0] * ch_pp.R / ch_pp.omega
                        #vp_y = exact_vel[1] * ch_pp.R / ch_pp.omega
                        #vp_z = exact_vel[2] * ch_pp.R / ch_pp.omega
                        #H[0] = H_old[0] + Dt_DEM / units_coefficient * exact_basset_force[0]
                        #H[1] = H_old[1] + Dt_DEM / units_coefficient * exact_basset_force[1]          
                        #sqrt_h =  math.sqrt(Dt_DEM)
                        #exact_Delta_H[0] = (H[0] - H_old[0]) / sqrt_h
                        #exact_Delta_H[1] = (H[1] - H_old[1]) / sqrt_h
                        ##Delta_H[0], Delta_H[1], Delta_H[2], present_coefficient = quad.DaitcheTimesAndIntegrands(times, integrands, 1)           
                        #Delta_H = [d / sqrt_h for d in Delta_H]
                        #basset_force[0] = units_coefficient * Dt_DEM_inv * (Delta_H[0])
                        #basset_force[1] = units_coefficient * Dt_DEM_inv * (Delta_H[1])
                        #print("\n integrands", integrands)
                        #print("delta_time", times[-1]-times[-2])
                        #print("DAITCHE", basset_force)                                                
                        #print("EXACT", exact_basset_force)                   
                        #print("\nDelta_H", Delta_H)
                        #print("exact_Delta_H", exact_Delta_H)     
                        #node.SetSolutionStepValue(BASSET_FORCE_X, basset_force[0])
                        #node.SetSolutionStepValue(BASSET_FORCE_Y, basset_force[1])     
                        #node.SetSolutionStepValue(BASSET_FORCE_Z, 0.) 
                        #print("\nINTEGRANDS", integrands)
                        #print("Delta_H",Delta_H)
                        #print("present_coefficient", present_coefficient)
                        #print("my n",len(times) - 2)
                    times.append(time_dem)
                    radii.append(r) 

        # performing the time integration of the DEM part

        spheres_model_part.ProcessInfo[TIME]    = time_dem
        rigid_face_model_part.ProcessInfo[TIME] = time_dem
        cluster_model_part.ProcessInfo[TIME]    = time_dem

        if not DEM_parameters["flow_in_porous_DEM_medium_option"].GetBool(): # in porous flow particles remain static      
            solver.Solve()
            #results_creator.Record(spheres_model_part, node_to_follow_id, time_dem)    
                
        # Walls movement:
        mesh_motion.MoveAllMeshes(rigid_face_model_part, time, Dt)
        mesh_motion.MoveAllMeshes(spheres_model_part, time, Dt)
        mesh_motion.MoveAllMeshes(DEM_inlet_model_part, time, Dt)

        #### TIME CONTROL ##################################
    
        # adding DEM elements by the inlet:
        if (DEM_parameters.dem_inlet_option):
            DEM_inlet.CreateElementsFromInletMesh(spheres_model_part, cluster_model_part, creator_destructor)  # After solving, to make sure that neighbours are already set.              
        
        #first_dem_iter = False

    #### PRINTING GRAPHS ####
    os.chdir(graphs_path)
    # measuring mean velocities in a certain control volume (the 'velocity trap')
    if (DEM_parameters.VelocityTrapOption):
        post_utils_DEM.ComputeMeanVelocitiesinTrap("Average_Velocity.txt", time)

    os.chdir(post_path)

    # applying DEM-to-fluid coupling

    if DEM_to_fluid_counter.Tick() and time >= DEM_parameters.interaction_start_time:
        print("Projecting from particles to the fluid...")
        sys.stdout.flush()
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
results_creator.PrintFile()
print("\n CALCULATIONS FINISHED. THE SIMULATION ENDED SUCCESSFULLY.")
simulation_elapsed_time = timer.clock() - simulation_start_time
print("Elapsed time: " + "%.5f"%(simulation_elapsed_time) + " s ")
print("per fluid time step: " + "%.5f"%(simulation_elapsed_time/ step) + " s ")
print("per DEM time step: " + "%.5f"%(simulation_elapsed_time/ DEM_step) + " s")
sys.stdout.flush()

dt_quad_over_dt = pp.CFD_DEM.delta_time_quadrature / pp.CFD_DEM.MaxTimeStep

with open('radii' + str(int(dt_quad_over_dt)) + '.txt','w') as f: 
    for i in range(len(radii)):
        f.write(str(times[i]) + ' ' + str(radii[i]) + '\n')
os.chdir(main_path)
sys.stdout.path_to_console_out_file
os.rename(sys.stdout.path_to_console_out_file, post_path + '/' + sys.stdout.path_to_console_out_file)
import shutil
folder_name = post_path + '_FINISHED_AT_t=' + str(round(time, 1))
try:
    shutil.rmtree(folder_name)
except OSError:
    pass

os.rename(post_path, folder_name)
final_position_file = open(folder_name + "/final_position", 'w')
final_position_file.write(str(x) + ' ' + str(y))
final_position_file.close()
from shutil import copyfile
copyfile(folder_name + "/final_position", main_path + "/final_position")
for i in drag_file_output_list:
    i.close()
