from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Defining variables so that we don't have to import them #
# Setting parameters
SolverType = "linear"
output_step = 1
Emissivity = 1.0
Convection_coefficient = 8.0
vel_type = 0
max_vel = 1.0
# Now setting the ProjectParameters
class ProjectParameters:
    Dt = 0.1
    max_time = 30.0
	# setting the domain size for the problem to be solved
    domain_size = 3
    LATENT_HEAT = 300e3
    SOLID_TEMPERATURE = 600.0
    FLUID_TEMPERATURE = 610.0
    AMBIENT_TEMPERATURE = 30
    DENSITY = 3000.0
    SPECIFIC_HEAT = 1000.0
    CONDUCTIVITY = 150.0
    HTC = 10
    MODULUS = 2.7355725535e-03
    body_force_x = 0.0000000000e+00
    body_force_y = -9.8100000000e+00
    body_force_z = 0.0000000000e+00
    wall_law_y = 0.005
    use_mass_correction = True
    reduction_on_failure = 0.3
    stabdt_pressure_factor = 1.0
    stabdt_convection_factor = 0.01
    tau2_factor = 1.0
    edge_factor = 5.0
    SmagorinskyConstant = 0.75
    edge_detection_angle = 45.0
    Alpha = 1.3256270216e-04

# setting the domain size for the problem to be solved

#
#
# ATTENTION: here the order is important
#
#

## Setting the path
# Declaring Problem name and path
import sys
import os
problem_name = "cube_benchmark"
current_directory=os.getcwd()
problem_path = current_directory
os.chdir("..\\..\\..\\..\\..")
kratos_path = os.getcwd()
os.chdir(current_directory)
#os.chdir(current_directory)


# including kratos path
kratos_libs_path =os.path.join(kratos_path,'libs')  # kratos_root/libs
kratos_applications_path = os.path.join(kratos_path,'applications')  #'\\applications'  # kratos_root/applications
kratos_benchmarking_path=os.path.join(kratos_path,'kratos\\benchmarking') #'\\kratos\\benchmarking'
kratos_kratos_path=os.path.join(kratos_path,'kratos') #'\\kratos\\benchmarking'

sys.path.append(kratos_path)
sys.path.append(kratos_libs_path)
sys.path.append(kratos_applications_path)
sys.path.append(kratos_benchmarking_path)
sys.path.append(kratos_kratos_path)
sys.path.append(problem_path)



# from now on the order is not anymore crucial

#
# Importing Kratos and extensions                                #
#

# importing Kratos main library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ConvectionDiffusionApplication import *

# We import the benchmarking utilities
import benchmarking

def BenchmarkTest(time,node1,node2):
    benchmarking.Output(time,"Time",1e-7)
    benchmarking.Output(node1.GetSolutionStepValue(TEMPERATURE),"TEST NODE 1 TEMPERATURE=",None, 0.005)
    benchmarking.Output(node2.GetSolutionStepValue(TEMPERATURE),"TEST NODE 2 TEMPERATURE=",None, 0.005)


def VelocityLaw(x, y, z, vel_type):
#
# Function defining the velocities inside the domain ##
# It returns a value from 0 to 1, and work in adim.  ##
# coordinates                                        ##
#
    v = Vector(3)
    v[0] = 0
    v[1] = 0
    v[2] = 0
    if vel_type == 1:
        v[0] = (1 - (z) * (z)) * (1 - (y) * (y))
    return v


# Creating mesh, we use the mdpa generator cube_mesher in   ##
# kratos/FluidDynamicsApplication/python_script/cube_mesher ##
#
from cube_mesher import *
h = 0.05
l = 0.05
nx = 1
ny = 1
nz = 1
box = box_data(0, -h / 2, -h / 2, l, h / 2, h / 2, nx, ny, nz)
# Defining element type. SUPGConvDiff3d is in ThermomechanicalApplication
elemtype = "SUPGConvDiffPhaseChange3D"  # "SUPGConvDiff3D" #"ConvDiff3D"#"VMS3D"
condtype = "EnvironmentContact3D"  # "Condition3D"#"WallCondition3D"
# Now we write the mdpa in the file, using cube_mesher functions
input_file_name = problem_path + "/" + problem_name
filename = input_file_name + ".mdpa"
with open(filename, "w") as mdpa:
    write_header(mdpa)
    generate_nodes(mdpa, box)
    generate_elements(mdpa, box, elemtype)
    generate_conditions(mdpa, box, condtype)
    generate_mesh_groups(mdpa, box)

#
# Creating a ModelPart         ##
#
# defining a model part for the fluid and one for the structure
model_part = ModelPart("FluidPart")

# We import the thermal solver
import thermal_solver_levelset

# Define variables and settings
config = ConvectionDiffusionSettings()
config.SetDensityVariable(DENSITY)
config.SetDiffusionVariable(CONDUCTIVITY)
config.SetUnknownVariable(TEMPERATURE)
config.SetVolumeSourceVariable(HEAT_FLUX)
config.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
config.SetMeshVelocityVariable(MESH_VELOCITY)
config.SetConvectionVariable(VELOCITY)
config.SetTransferCoefficientVariable(HTC)

#
# Now we load the mdpa file the we have just created ##
#

# introducing input file name
input_file_name = problem_name

# reading the fluid part
gid_mode = GiDPostMode.GiD_PostBinary
multifile = MultiFileFlag.MultipleFiles
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteConditions  # WriteElementsOnly
gid_io = GidIO(input_file_name, gid_mode,
               multifile, deformed_mesh_flag, write_conditions)

#
# Now we import the Solver                              ##
#
# We Add the Variables to be needed to the ModelPart
thermal_solver_levelset.AddVariables(model_part, config)

model_part.AddNodalSolutionStepVariable(CONDUCTIVITY)
model_part.AddNodalSolutionStepVariable(TEMPERATURE)
model_part.AddNodalSolutionStepVariable(HEAT_FLUX)
model_part.AddNodalSolutionStepVariable(FACE_HEAT_FLUX)
model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT)
model_part.AddNodalSolutionStepVariable(TEMPERATURES)
model_part.AddNodalSolutionStepVariable(NODAL_MASS)
model_part.AddNodalSolutionStepVariable(FRONT_MEETING)
model_part.AddNodalSolutionStepVariable(MACRO_POROSITY)
model_part.AddNodalSolutionStepVariable(SHRINKAGE_POROSITY)
model_part.AddNodalSolutionStepVariable(FILLTIME)
model_part.AddNodalSolutionStepVariable(DP_ALPHA1)
model_part.AddNodalSolutionStepVariable(REACTION)
model_part.AddNodalSolutionStepVariable(ENTHALPY)

# NOW WE READ THE MODELPART
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
model_part.SetBufferSize(3)

# Adding dofs to the modelpart
thermal_solver_levelset.AddDofs(model_part, config)


# model_part.Properties[0][CONVECTION_COEFFICIENT] = Convection_coefficient
model_part.ProcessInfo.SetValue(
    AMBIENT_TEMPERATURE, ProjectParameters.AMBIENT_TEMPERATURE)
model_part.ProcessInfo.SetValue(
    FLUID_TEMPERATURE, ProjectParameters.FLUID_TEMPERATURE)
model_part.ProcessInfo.SetValue(
    SOLID_TEMPERATURE, ProjectParameters.SOLID_TEMPERATURE)
model_part.ProcessInfo.SetValue(LATENT_HEAT, ProjectParameters.LATENT_HEAT)
model_part.ProcessInfo.SetValue(DENSITY, ProjectParameters.DENSITY)
model_part.ProcessInfo.SetValue(SPECIFIC_HEAT, ProjectParameters.SPECIFIC_HEAT)
model_part.ProcessInfo.SetValue(HTC, ProjectParameters.HTC)
model_part.ProcessInfo.SetValue(K0, ProjectParameters.MODULUS)
model_part.ProcessInfo.SetValue(DELTA_TIME, ProjectParameters.Dt)
# model_part.ProcessInfo.SetValue(IS_GRAVITY_FILLING, is_gravity_filling)

#import real_tables as input_tables
import constant_tables as input_tables
input_tables.InitializeTables(model_part)

##start everything at 800 and let's see how it gets cold
for node in model_part.Nodes:
    node.SetSolutionStepValue(TEMPERATURE,0,800)
    
# Generating the Solver
solver = thermal_solver_levelset.Solver(model_part, ProjectParameters.domain_size, config)
solver.MaxNewtonRapshonIterations = 15
solver.echo_level = 1
print("Before initialize")
solver.Initialize()
print("Thermal solver created :)")

# settings to be changed
full_Dt = ProjectParameters.Dt
final_time = ProjectParameters.max_time

def ComputeBDFcoefficients(model_part):
    Dt = model_part.ProcessInfo[DELTA_TIME]
    OldDt = Dt #model_part.ProcessInfo.GetPreviousTimeStepInfo(1)[DELTA_TIME]

    Rho = OldDt / Dt
    TimeCoeff = 1.0 / (Dt * Rho * Rho + Dt * Rho)

    BDFcoeffs = Vector(3)
    BDFcoeffs[0] = TimeCoeff * (Rho * Rho + 2.0 * Rho)
    BDFcoeffs[1] = -TimeCoeff * (Rho * Rho + 2.0 * Rho + 1.0)
    BDFcoeffs[2] = TimeCoeff
    model_part.ProcessInfo.SetValue(BDF_COEFFICIENTS, BDFcoeffs)


for node in model_part.Nodes:
    x = node.X / l
    y = 2 * node.Y / h
    z = 2 * node.Z / h
# Back Face (X min)
for node in model_part.GetMesh(2).Nodes:
    a = 0

# Side Faces
sidefaces = []
for ii in range(3, 7):
    for cond in model_part.GetMesh(ii).Conditions:
        sidefaces.append(cond)
# mesh to be printed
mesh_name = 0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh(model_part.GetMesh())
gid_io.FinalizeMesh()
gid_io.InitializeResults(mesh_name, (model_part).GetMesh())

solver.solver.Check()
out = 0
time = 0.0
step = 0
print("***************************")
print("*** Starting  Time Loop ***")
print("***************************")

for node in model_part.Nodes:
   node.SetSolutionStepValue(DISTANCE,0,-1.0)

#Now we create the objects containing nodes to be checked
node1=model_part.Nodes[1]
node2=model_part.Nodes[8]

while(time < final_time):
    time = time + full_Dt
    print("Time = " + str(time) + " out of " + str(final_time))
    model_part.CloneTimeStep(time)
    if(step >= 3):
        ComputeBDFcoefficients(model_part)
        solver.Solve()

    if(out == output_step):
        gid_io.WriteNodalResults(TEMPERATURE, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(CONDUCTIVITY, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(SPECIFIC_HEAT, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DENSITY, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(VELOCITY, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(ENTHALPY,model_part.Nodes,time,0)
        gid_io.WriteNodalResults(HEAT_FLUX, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(FACE_HEAT_FLUX, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(SOLIDFRACTION, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(DISTANCE, model_part.Nodes, time, 0)
        gid_io.Flush()
        out = 0
    
    out = out + 1
    step = step + 1
    BenchmarkTest(time,node1,node2)

gid_io.FinalizeResults()
