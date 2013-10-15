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
import DEM_explicit_solver_var as DEM_parameters
import DEM_procedures as DEMProc
import swimming_DEM_procedures as SwimProc

# PROJECT PARAMETERS TO BE PUT IN PROBLEMTYPE
ProjectParameters.ProjectionModuleOption      = 1
ProjectParameters.PrintParticlesResultsOption = 0
ProjectParameters.ProjectFromParticlesOption  = 1
ProjectParameters.ProjectAtEverySubStepOption = 1
ProjectParameters.CreateParticlesOption       = 0
ProjectParameters.CalculatePorosity           = 1
ProjectParameters.Interaction_start_time      = 0.01
ProjectParameters.gravity_x                   = 0.00000e+00
ProjectParameters.gravity_y                   = 0.0
ProjectParameters.gravity_z                   = -9.81000e+00
DEM_parameters.GravityX                       = ProjectParameters.gravity_x
DEM_parameters.GravityY                       = ProjectParameters.gravity_y
DEM_parameters.GravityZ                       = ProjectParameters.gravity_z
ProjectParameters.plastic_viscosity           = 0.000014 #It is divided by the density
ProjectParameters.smoothing_parameter_m       = 0.035
ProjectParameters.yield_stress_value          = 10.0
ProjectParameters.DEM_nodal_results           = ["RADIUS", "FLUID_VEL_PROJECTED", "DRAG_FORCE", "BUOYANCY", "PRESSURE_GRAD_PROJECTED"]
ProjectParameters.mixed_nodal_results         = ["VELOCITY", "DISPLACEMENT"]
ProjectParameters.nodal_results.append("SOLID_FRACTION")

for var in ProjectParameters.mixed_nodal_results:

    if var in ProjectParameters.nodal_results:
        ProjectParameters.nodal_results.remove(var)

# Constructing a DEM_procedures object
DEM_proc = DEMProc.Procedures(DEM_parameters)
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

# defining a model part for the fluid
fluid_model_part = ModelPart("FluidPart")

if "REACTION" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(REACTION)
if "DISTANCE" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#if "MU" in ProjectParameters.nodal_results:
#    fluid_model_part.AddNodalSolutionStepVariable(MU)
##AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#
#
# importing the solvers needed
SolverSettings = ProjectParameters.FluidSolverConfiguration
solver_module = import_solver(SolverSettings)

#
# importing variables
solver_module.AddVariables(fluid_model_part, SolverSettings)

# introducing input file name
input_file_name = ProjectParameters.problem_name

# reading the fluid part
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
fluid_model_part.AddNodalSolutionStepVariable(PRESSURE_GRADIENT)
fluid_model_part.AddNodalSolutionStepVariable(AUX_DOUBLE_VAR)
fluid_model_part.AddNodalSolutionStepVariable(DRAG_REACTION)
fluid_model_part.AddNodalSolutionStepVariable(SOLID_FRACTION)

# Defining a model part for the balls part
my_timer = Timer()

balls_model_part = ModelPart("SolidPart")

# HYDRODYNAMICS
balls_model_part.AddNodalSolutionStepVariable(FLUID_VEL_PROJECTED)
balls_model_part.AddNodalSolutionStepVariable(FLUID_DENSITY_PROJECTED)
balls_model_part.AddNodalSolutionStepVariable(PRESSURE_GRAD_PROJECTED)
balls_model_part.AddNodalSolutionStepVariable(FLUID_VISCOSITY_PROJECTED)

# FORCES
balls_model_part.AddNodalSolutionStepVariable(DRAG_FORCE)
balls_model_part.AddNodalSolutionStepVariable(BUOYANCY)

# Defining a model part for the mixed part
mixed_model_part = ModelPart("MixedPart")

import sphere_strategy as SolverStrategy
SolverStrategy.AddVariables(balls_model_part, DEM_parameters)

# reading the balls model part
model_part_io_solid = ModelPartIO(DEM_parameters.problem_name)
model_part_io_solid.ReadModelPart(balls_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
balls_model_part.SetBufferSize(3)

# adding nodal degrees of freedom
SolverStrategy.AddDofs(balls_model_part)
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
#
#

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
import swimming_DEM_gid_output as gid_output #S
swimming_DEM_gid_io = gid_output.SwimmingDEMGiDOutput(input_file_name, #SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
                   ProjectParameters.VolumeOutput,
                   ProjectParameters.GiDPostMode,
                   ProjectParameters.GiDMultiFileFlag,
                   ProjectParameters.GiDWriteMeshFlag,
                   ProjectParameters.GiDWriteConditionsFlag)

swimming_DEM_gid_io.initialize_swimming_DEM_results(balls_model_part, mixed_model_part)
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

# Stepping and time settings
Dt = ProjectParameters.Dt
Nsteps = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0
#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
full_Dt = Dt
initial_Dt = 0.001 * full_Dt

# creating projection modeule for the fluid-DEM coupling and initial projection
h_min = 0.01
n_balls = 1
fluid_volume = 10
n_particles_in_depth = int(math.sqrt(n_balls / fluid_volume))

if (ProjectParameters.ProjectionModuleOption):
    projection_module = SwimProc.ProjectionModule(fluid_model_part, balls_model_part, domain_size, n_particles_in_depth)
    projection_module.UpdateDatabase(h_min)

fluid_model_part.ProcessInfo.SetValue(YIELD_STRESS, ProjectParameters.yield_stress_value )
fluid_model_part.ProcessInfo.SetValue(M, ProjectParameters.smoothing_parameter_m)

for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(VISCOSITY, 0, ProjectParameters.plastic_viscosity)
    node.SetSolutionStepValue(BODY_FORCE_X, 0, ProjectParameters.gravity_x)
    node.SetSolutionStepValue(BODY_FORCE_Y, 0, ProjectParameters.gravity_y)
    node.SetSolutionStepValue(BODY_FORCE_Z, 0, ProjectParameters.gravity_z)

creator_destructor = ParticleCreatorDestructor()

#balls_model_part.SetBufferSize(3)
DEM_solver = SolverStrategy.ExplicitStrategy(balls_model_part, creator_destructor, DEM_parameters)
DEM_proc.GiDSolverTransfer(balls_model_part, DEM_solver, DEM_parameters)



# Defining a model part for the DEM inlet  #############################################################MA

DEM_inlet_model_part = ModelPart("DEMInletPart") #MA
DEM_Inlet_filename=DEM_parameters.problem_name+"_Inlet" #MA
#DEM_Inlet_filename="1DEM_Inlet" #MA
SolverStrategy.AddVariables(DEM_inlet_model_part, DEM_parameters) #MA
# HYDRODYNAMICS
DEM_inlet_model_part.AddNodalSolutionStepVariable(FLUID_VEL_PROJECTED) #MA
DEM_inlet_model_part.AddNodalSolutionStepVariable(FLUID_DENSITY_PROJECTED) #MA
DEM_inlet_model_part.AddNodalSolutionStepVariable(PRESSURE_GRAD_PROJECTED) #MA
DEM_inlet_model_part.AddNodalSolutionStepVariable(FLUID_VISCOSITY_PROJECTED) #MA
# FORCES
DEM_inlet_model_part.AddNodalSolutionStepVariable(DRAG_FORCE) #MA
DEM_inlet_model_part.AddNodalSolutionStepVariable(BUOYANCY) #MA


model_part_io_demInlet = ModelPartIO(DEM_Inlet_filename) #MA
model_part_io_demInlet.ReadModelPart(DEM_inlet_model_part) #MA
DEM_inlet_model_part.SetBufferSize(2) #MA
SolverStrategy.AddDofs(DEM_inlet_model_part)
DEM_inlet_parameters = DEM_inlet_model_part.Properties  #MA
DEM_inlet = DEM_Inlet(DEM_inlet_model_part) #MA Constructor for an Inlet
DEM_inlet.InitializeDEM_Inlet(balls_model_part,creator_destructor) #MA
DEM_step = 0   # MA #this variable is necessary to get a good random insertion of particles
########################################################################################################MA




DEM_solver.Initialize()

def my_time_dem(time_dem_initial, time_final, Dt_DEM):

    while time_dem_initial < time_final:
        yield time_dem_initial
        time_dem_initial += Dt_DEM

    yield time_final

#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

while(time <= final_time):
    
    if(step < 3):
        Dt = initial_Dt

    else:
        Dt = full_Dt
        
    print "\n" , "TIME = " , time 
     
    step     = step + 1
    time_dem = time
    time     = time + Dt
    fluid_model_part.CloneTimeStep(time)

    if(step < 3):
        balls_model_part.CloneTimeStep(time)

    if(step >= 3):
        
        print "Solving Fluid... (" , fluid_model_part.NumberOfElements(0) , " elements)"
        fluid_solver.Solve()
#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
        if (time > ProjectParameters.Interaction_start_time):
            time_final = time
            Dt_DEM     = DEM_parameters.MaxTimeStep

            if (ProjectParameters.ProjectionModuleOption):
                interaction_calculator = CustomFunctionsCalculator()
                interaction_calculator.PressureGradientCalculator(fluid_model_part)
                projection_module.ProjectFromFluid((time_dem + Dt - time) / Dt)
                    
            print "Solving DEM...(" , balls_model_part.NumberOfElements(0) , " elements)"        
            for time_dem in my_time_dem(time_dem, time_final, Dt_DEM):
				
                #print "DEM time = " , time_dem
                DEM_step = DEM_step + 1 #MA
                balls_model_part.ProcessInfo[TIME_STEPS] = DEM_step #MA
                
                if (ProjectParameters.ProjectionModuleOption and ProjectParameters.ProjectAtEverySubStepOption):
                    projection_module.ProjectFromFluid((time_dem + Dt - time) / Dt)
                
                balls_model_part.CloneTimeStep(time_dem)
                DEM_solver.Solve()
                DEM_inlet.CreateElementsFromInletMesh( balls_model_part, DEM_inlet_model_part, creator_destructor, "SphericSwimmingParticle3D" ) #MA #After solving, to make sure that neighbours are already set.
    # options are: "SphericParticle3D" , "SphericSwimmingParticle3D" 

            if (ProjectParameters.ProjectionModuleOption):

                if (ProjectParameters.ProjectFromParticlesOption):
                    projection_module.ProjectFromParticles()

        if (ProjectParameters.PrintParticlesResultsOption):
            print_particles_results = PrintParticlesResults("DRAG_FORCE", time, balls_model_part)
            print_particles_results = PrintParticlesResults("BUOYANCY", time, balls_model_part)
            print_particles_results = PrintParticlesResults("VELOCITY", time, balls_model_part)

#AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

        graph_printer.PrintGraphs(time)
        PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)

    if(output_time <= out):
        ParticleUtils2D().VisualizationModelPart(mixed_model_part, fluid_model_part, balls_model_part)
#       gid_io.write_results(  # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#           time,               # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#           fluid_model_part,    # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#           ProjectParameters.nodal_results,    # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#           ProjectParameters.gauss_points_results)    # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

# SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
        swimming_DEM_gid_io.write_swimming_DEM_results(
            time,
            fluid_model_part,
            balls_model_part,
            mixed_model_part,
            ProjectParameters.nodal_results,
            ProjectParameters.DEM_nodal_results,
            ProjectParameters.mixed_nodal_results,
            ProjectParameters.gauss_points_results)
# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

        out = 0

    out = out + Dt

#gid_io.finalize_results()  # SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
swimming_DEM_.finalize_results()

for i in drag_file_output_list:
    i.close()
