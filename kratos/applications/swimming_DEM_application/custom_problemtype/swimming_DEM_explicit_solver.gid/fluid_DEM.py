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

#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
import math

from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import DEM_explicit_solver_var as DEMParameters
import DEM_procedures
import swimming_DEM_procedures

# PROJECT PARAMETERS (to be put in problem type)
ProjectParameters.projection_module_option         = 1
ProjectParameters.print_particles_results_option   = 0
ProjectParameters.project_from_particles_option    = 0
ProjectParameters.project_at_every_substep_option  = 0
ProjectParameters.velocity_trap_option             = 0
ProjectParameters.inlet_option                     = 1
ProjectParameters.non_newtonian_option             = 0
ProjectParameters.manually_imposed_drag_law_option = 0
ProjectParameters.similarity_transformation_type   = 0 # no transformation (0), Tsuji (1)
ProjectParameters.dem_inlet_element_type           = "SphericSwimmingParticle3D"  # "SphericParticle3D", "SphericSwimmingParticle3D"
ProjectParameters.coupling_scheme_type             = "UpdatedFluid" # "UpdatedFluid", "UpdatedDEM"
ProjectParameters.coupling_weighing_type           = 2 # {fluid_to_DEM, DEM_to_fluid, Solid_fraction} = {lin, const, const} (0), {lin, lin, const} (1), {lin, lin, lin} (2)
ProjectParameters.buoyancy_force_type              = 1 # null buoyancy (0), standard (1)
ProjectParameters.drag_force_type                  = 1 # null drag (0), standard (1), Weatherford (2), Ganser (3)
ProjectParameters.virtual_mass_force_type          = 0 # null virtual mass force (0)
ProjectParameters.lift_force_type                  = 0 # null lift force (0)
ProjectParameters.drag_modifier_type               = 3 # Hayder (2), Chien (3)
ProjectParameters.interaction_start_time           = 0.00
ProjectParameters.gravity_x                        = 0.0
ProjectParameters.gravity_y                        = 0.0
ProjectParameters.gravity_z                        = 0.0 #- 9.81
ProjectParameters.smoothing_parameter_m            = 0.035
ProjectParameters.yield_stress_value               = 0.0
ProjectParameters.max_solid_fraction               = 0.6
ProjectParameters.gel_strength                     = 0.0
ProjectParameters.power_law_n                      = 0.0
ProjectParameters.power_law_k                      = 0.0
ProjectParameters.initial_drag_force               = 0.0
ProjectParameters.drag_law_slope                   = 0.0
ProjectParameters.power_law_tol                    = 0.0
ProjectParameters.model_over_real_diameter_factor  = 2.0 # not active if similarity_transformation_type = 0

# variables to be printed
ProjectParameters.dem_nodal_results                = ["RADIUS", "FLUID_VEL_PROJECTED", "DRAG_FORCE", "BUOYANCY", "PRESSURE_GRAD_PROJECTED", "REYNOLDS_NUMBER"]
ProjectParameters.mixed_nodal_results              = ["VELOCITY", "DISPLACEMENT"]
ProjectParameters.variables_to_print_in_file       = ["DRAG_FORCE", "BUOYANCY", "VELOCITY"]

# changes on PROJECT PARAMETERS for the sake of consistency
ProjectParameters.nodal_results.append("SOLID_FRACTION")
ProjectParameters.nodal_results.append("MESH_VELOCITY1")
ProjectParameters.nodal_results.append("BODY_FORCE")
ProjectParameters.nodal_results.append("DRAG_REACTION")

ProjectParameters.project_from_particles_option *= ProjectParameters.projection_module_option
ProjectParameters.project_at_every_substep_option *= ProjectParameters.projection_module_option

DEMParameters.GravityX                       = ProjectParameters.gravity_x
DEMParameters.GravityY                       = ProjectParameters.gravity_y
DEMParameters.GravityZ                       = ProjectParameters.gravity_z

for var in ProjectParameters.mixed_nodal_results:

    if var in ProjectParameters.nodal_results:
        ProjectParameters.nodal_results.remove(var)

# extra nodal variables (related to the coupling) to be added to the model parts (memory will be allocated for them)
fluid_variables_to_add = [PRESSURE_GRADIENT,
                          AUX_DOUBLE_VAR,
                          DRAG_REACTION,
                          SOLID_FRACTION,
                          MESH_VELOCITY1]

balls_variables_to_add = [FLUID_VEL_PROJECTED,
                          FLUID_DENSITY_PROJECTED,
                          PRESSURE_GRAD_PROJECTED,
                          FLUID_VISCOSITY_PROJECTED,
                          DRAG_FORCE,
                          BUOYANCY,
                          SOLID_FRACTION_PROJECTED,
                          REYNOLDS_NUMBER]

fem_dem_variables_to_add = [VELOCITY,
                            DISPLACEMENT]

DEM_inlet_variables_to_add = balls_variables_to_add

# constructing a DEM_procedures object
dem_procedures = DEM_procedures.Procedures(DEMParameters)
##AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

# defining variables to be used
# GID IO IS NOT USING THIS NOW. TO BE REMOVED ONCE THE "PRINT IN POINTS"
# CODE IS NOT USING IT

variables_dictionary = {"PRESSURE"   : PRESSURE,
                        "VELOCITY"   : VELOCITY,
                        "MU"         : MU,         #SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
                        "BUOYANCY"   : BUOYANCY,   #SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
                        "DRAG_FORCE" : DRAG_FORCE} #SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
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
print 'Adding nodal variables to the fluid_model_part' #(memory allocation) #SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

solver_module.AddVariables(fluid_model_part, SolverSettings)
swimming_DEM_procedures.AddNodalVariables(fluid_model_part, fluid_variables_to_add) #SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

# introducing input file name
input_file_name = ProjectParameters.problem_name

# reading the fluid part
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS


# defining model parts for the balls part and for the DEM-FEM interaction elements

balls_model_part = ModelPart("SolidPart")
fem_dem_model_part = ModelPart("RigidFace_Part");

print 'Adding nodal variables to the balls_model_part' #(memory allocation)

swimming_DEM_procedures.AddNodalVariables(balls_model_part, balls_variables_to_add)

print 'Adding nodal variables to the dem_fem_wall_model_part' #(memory allocation)

swimming_DEM_procedures.AddNodalVariables(fem_dem_model_part, fem_dem_variables_to_add)

# defining a model part for the mixed part
mixed_model_part = ModelPart("MixedPart")

import sphere_strategy as SolverStrategy
SolverStrategy.AddVariables(balls_model_part, DEMParameters)

# reading the balls model part
model_part_io_solid = ModelPartIO(DEMParameters.problem_name)
model_part_io_solid.ReadModelPart(balls_model_part)
rigidFace_mp_filename = DEMParameters.problem_name + "_FEM_boundary"
model_part_io_solid = ModelPartIO(rigidFace_mp_filename)
model_part_io_solid.ReadModelPart(fem_dem_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
balls_model_part.SetBufferSize(1)

# adding nodal degrees of freedom
SolverStrategy.AddDofs(balls_model_part)

# adding extra process info variables
balls_model_part.ProcessInfo.SetValue(BUOYANCY_FORCE_TYPE, ProjectParameters.buoyancy_force_type)
balls_model_part.ProcessInfo.SetValue(DRAG_FORCE_TYPE, ProjectParameters.drag_force_type)
balls_model_part.ProcessInfo.SetValue(VIRTUAL_MASS_FORCE_TYPE, ProjectParameters.virtual_mass_force_type)
balls_model_part.ProcessInfo.SetValue(LIFT_FORCE_TYPE, ProjectParameters.lift_force_type)
balls_model_part.ProcessInfo.SetValue(NON_NEWTONIAN_OPTION, ProjectParameters.non_newtonian_option)
balls_model_part.ProcessInfo.SetValue(MANUALLY_IMPOSED_DRAG_LAW_OPTION, ProjectParameters.manually_imposed_drag_law_option)
balls_model_part.ProcessInfo.SetValue(DRAG_MODIFIER_TYPE, ProjectParameters.drag_modifier_type)
balls_model_part.ProcessInfo.SetValue(GEL_STRENGTH, ProjectParameters.gel_strength)
balls_model_part.ProcessInfo.SetValue(POWER_LAW_N, ProjectParameters.power_law_n)
balls_model_part.ProcessInfo.SetValue(POWER_LAW_K, ProjectParameters.power_law_k)
balls_model_part.ProcessInfo.SetValue(INIT_DRAG_FORCE, ProjectParameters.initial_drag_force)
balls_model_part.ProcessInfo.SetValue(DRAG_LAW_SLOPE, ProjectParameters.drag_law_slope)
balls_model_part.ProcessInfo.SetValue(POWER_LAW_TOLERANCE, ProjectParameters.power_law_tol)

fluid_model_part.ProcessInfo.SetValue(YIELD_STRESS, ProjectParameters.yield_stress_value)
fluid_model_part.ProcessInfo.SetValue(M, ProjectParameters.smoothing_parameter_m)
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

# adding nodal degrees of freedom
solver_module.AddDofs(fluid_model_part, SolverSettings)

# If Lalplacian form = 2, free all pressure Dofs
#laplacian_form = ProjectParameters.laplacian_form
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
        ##get the nodes of the wall for SA.
        nodes = fluid_model_part.GetNodes(i)
        for node in nodes:
            fluid_solver.wall_nodes.append(node)
            node.SetSolutionStepValue(TURBULENT_VISCOSITY, 0, 0.0)
            node.Fix(TURBULENT_VISCOSITY)


fluid_solver.Initialize()
print "fluid solver created"

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
swimming_DEM_gid_io = swimming_DEM_gid_output.SwimmingDEMGiDOutput(input_file_name,
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

print drag_file_output_list

def PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time):
    i = 0
    for it in drag_list:
        print it[0]
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

#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
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
    projection_module = swimming_DEM_procedures.ProjectionModule(fluid_model_part, balls_model_part, fem_dem_model_part, domain_size, ProjectParameters.max_solid_fraction, ProjectParameters.coupling_weighing_type, n_particles_in_depth)
    projection_module.UpdateDatabase(h_min)
    interaction_calculator = CustomFunctionsCalculator()

# TEMPORARY!! setting the initial kinematic viscosity, should be done from problemtype
for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(BODY_FORCE_X, 0, ProjectParameters.gravity_x)
    node.SetSolutionStepValue(BODY_FORCE_Y, 0, ProjectParameters.gravity_y)
    node.SetSolutionStepValue(BODY_FORCE_Z, 0, ProjectParameters.gravity_z)

# creating a CreatorDestructor object, encharged of any adding or removing of elements during the simulation
creator_destructor = ParticleCreatorDestructor()
creator_destructor.SetMaxNodeId(max_fluid_node_Id)

# creating a Solver object for the DEM part. It contains the sequence of function calls necessary for the evolution of the DEM system at every time step
dem_solver = SolverStrategy.ExplicitStrategy(balls_model_part, fem_dem_model_part, creator_destructor, DEMParameters)

# constructing a model part for the DEM inlet. it contains the DEM elements to be released during the simulation

if (ProjectParameters.inlet_option):
    DEM_inlet_model_part = ModelPart("DEMInletPart")
    DEM_Inlet_filename = DEMParameters.problem_name + "_Inlet"
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
    DEM_inlet_model_part.ProcessInfo.SetValue(VIRTUAL_MASS_FORCE_TYPE, ProjectParameters.virtual_mass_force_type)
    DEM_inlet_model_part.ProcessInfo.SetValue(LIFT_FORCE_TYPE, ProjectParameters.lift_force_type)
    DEM_inlet_model_part.ProcessInfo.SetValue(NON_NEWTONIAN_OPTION, ProjectParameters.non_newtonian_option)
    DEM_inlet_model_part.ProcessInfo.SetValue(MANUALLY_IMPOSED_DRAG_LAW_OPTION, ProjectParameters.manually_imposed_drag_law_option)
    DEM_inlet_model_part.ProcessInfo.SetValue(DRAG_MODIFIER_TYPE, ProjectParameters.drag_modifier_type)
    DEM_inlet_model_part.ProcessInfo.SetValue(GEL_STRENGTH, ProjectParameters.gel_strength)
    DEM_inlet_model_part.ProcessInfo.SetValue(POWER_LAW_N, ProjectParameters.power_law_n)
    DEM_inlet_model_part.ProcessInfo.SetValue(POWER_LAW_K, ProjectParameters.power_law_k)
    DEM_inlet_model_part.ProcessInfo.SetValue(INIT_DRAG_FORCE, ProjectParameters.initial_drag_force)
    DEM_inlet_model_part.ProcessInfo.SetValue(DRAG_LAW_SLOPE, ProjectParameters.drag_law_slope)
    DEM_inlet_model_part.ProcessInfo.SetValue(POWER_LAW_TOLERANCE, ProjectParameters.power_law_tol)

    # constructiong the inlet and intializing it
    DEM_inlet = DEM_Inlet(DEM_inlet_model_part)
    DEM_inlet.InitializeDEM_Inlet(balls_model_part, creator_destructor)

def yield_DEM_time(current_time, current_time_plus_increment, delta_time):
    current_time += delta_time

    while current_time < current_time_plus_increment - delta_time:
        yield current_time
        current_time += delta_time

    current_time = current_time_plus_increment
    yield current_time
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

# renumerating IDs if required

#swimming_DEM_procedures.RenumberNodesIdsToAvoidRepeating(fluid_model_part, balls_model_part, fem_dem_model_part)

# Stepping and time settings

Dt = ProjectParameters.Dt
Nsteps = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = Dt
step = 0

#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
time_dem = 0.0
DEM_step = 0 # this variable is necessary to get a good random insertion of particles
dem_solver.Initialize()
Dt_DEM = DEMParameters.MaxTimeStep
#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

while(time <= final_time):
    
    time = time + Dt
    step = step + 1
    fluid_model_part.CloneTimeStep(time)
    print "\n", "TIME = ", time

    if (ProjectParameters.coupling_scheme_type == "UpdatedDEM"):
        time_final_DEM_substepping = time + Dt

    else:
        time_final_DEM_substepping = time

    # solving the fluid part

    if (step >= 3):
        print "Solving Fluid... (", fluid_model_part.NumberOfElements(0), " elements )"
        fluid_solver.Solve()

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

    print "Solving DEM... (", balls_model_part.NumberOfElements(0), " elements)"

    if (time >= ProjectParameters.interaction_start_time and ProjectParameters.projection_module_option):
        #print "Calculating Pressure Gradient ..."
        interaction_calculator.CalculatePressureGradient(fluid_model_part)                 
            
    for time_dem in yield_DEM_time(time_dem, time_final_DEM_substepping, Dt_DEM):
        
        DEM_step = DEM_step + 1   # this variable is necessary to get a good random insertion of particles
        
        balls_model_part.ProcessInfo[TIME_STEPS] = DEM_step
        
        # applying fluid-to-DEM coupling
                
        if (time >= ProjectParameters.interaction_start_time and ProjectParameters.project_at_every_substep_option):
            
            if (ProjectParameters.coupling_scheme_type == "UpdatedDEM"):
                projection_module.ProjectFromNewestFluid()

            else:
                projection_module.ProjectFromFluid((time_final_DEM_substepping - time_dem) / Dt)

        # performing the time integration of the DEM part

        balls_model_part.CloneTimeStep(time_dem)
        dem_solver.Solve()

        # adding DEM elements by the inlet

        if (ProjectParameters.inlet_option):
            DEM_inlet.CreateElementsFromInletMesh(balls_model_part, DEM_inlet_model_part, creator_destructor, ProjectParameters.dem_inlet_element_type) #After solving, to make sure that neighbours are already set.

        # measuring mean velocities in a certain control volume (the 'velocity trap')

    if (ProjectParameters.velocity_trap_option):
        post_utils.ComputeMeanVelocitiesinTrap("Average_Velocity", time_dem)

    # applying DEM-to-fluid coupling

    if (time >= ProjectParameters.interaction_start_time and ProjectParameters.project_from_particles_option):
        print "Project from particles to the fluid"
        projection_module.ProjectFromParticles()

    # printing if required

    if (ProjectParameters.print_particles_results_option):
        io_tools.PrintParticlesResults(ProjectParameters.variables_to_print_in_file, time, balls_model_part)
        graph_printer.PrintGraphs(time)
        PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)
        
    if (output_time <= out and ProjectParameters.coupling_scheme_type == "UpdatedFluid"):
        
        if (ProjectParameters.projection_module_option):
            projection_module.ComputePostProcessResults(balls_model_part.ProcessInfo) 
            
        print ""
        print "*******************  PRINTING RESULTS FOR GID  ***************************" 

        if (ProjectParameters.GiDMultiFileFlag == "Multiples"):
            mixed_model_part.Elements.clear()
            mixed_model_part.Nodes.clear()
            post_utilities.AddModelPartToModelPart(mixed_model_part, fluid_model_part)
            post_utilities.AddModelPartToModelPart(mixed_model_part, balls_model_part)
            post_utilities.AddModelPartToModelPart(mixed_model_part, fem_dem_model_part)

        swimming_DEM_gid_io.write_swimming_DEM_results(time, fluid_model_part, balls_model_part, fem_dem_model_part, mixed_model_part, ProjectParameters.nodal_results, ProjectParameters.dem_nodal_results, ProjectParameters.mixed_nodal_results, ProjectParameters.gauss_points_results)
        out = 0
# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

#       gid_io.write_results(  # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#           time,               # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#           fluid_model_part,    # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#           ProjectParameters.nodal_results,    # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#           ProjectParameters.gauss_points_results)    # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

        out = 0

    out = out + Dt

#gid_io.finalize_results()  # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
swimming_DEM_gid_io.finalize_results()
#
print "CALCULATIONS FINISHED. THE SIMULATION ENDED SUCCESSFULLY."
#
for i in drag_file_output_list:
    i.close()
