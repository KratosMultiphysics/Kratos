from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

######################################################################################
######################################################################################
######################################################################################

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

# Parsing the parameters
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = KratosMultiphysics.Parameters( parameter_file.read())

## Defining variables ----------------------------------------------------------------------------------------

domain_size         = ProjectParameters["problem_data"]["domain_size"].GetInt()
problem_name        = ProjectParameters["problem_data"]["problem_name"].GetString()
echo_level          = ProjectParameters["solver_settings"]["echo_level"].GetInt()
buffer_size         = ProjectParameters["solver_settings"]["buffer_size"].GetInt()
end_time            = ProjectParameters["problem_data"]["end_time"].GetDouble()
time                = ProjectParameters["problem_data"]["start_time"].GetDouble()
gravity             = ProjectParameters["problem_data"]["gravity"].GetDouble()
time_scale          = ProjectParameters["problem_data"]["time_scale"].GetString()
water_height_scale  = ProjectParameters["problem_data"]["water_height_scale"].GetString()

# Time_unit_converter = 1.0
if   time_scale == "seconds":
    time_unit_converter =     1
elif time_scale == "minutes":
    time_unit_converter =    60
elif time_scale == "hours":
    time_unit_converter =  3600
elif time_scale == "days":
    time_unit_converter = 86400
else:
    raise Exception("unknown time scale")

# Water height unit converter
if   water_height_scale == "meters":
    water_height_unit_converter = 1.0
elif water_height_scale == "millimeters":
    water_height_unit_converter = 0.001
else:
    raise Exception("unknown water height scale")

## Model part ------------------------------------------------------------------------------------------------

# Defining the model part
main_model_part = KratosMultiphysics.ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
main_model_part.ProcessInfo.SetValue(KratosMultiphysics.GRAVITY_Z, gravity * time_unit_converter**2)
main_model_part.ProcessInfo.SetValue(KratosShallow.TIME_UNIT_CONVERTER, time_unit_converter)
main_model_part.ProcessInfo.SetValue(KratosShallow.WATER_HEIGHT_UNIT_CONVERTER, water_height_unit_converter)

# Construct the solver (main settings methods are located in the solver_module)
solver_module = __import__(ProjectParameters["solver_settings"]["solver_type"].GetString())
solver = solver_module.CreateSolver(main_model_part, ProjectParameters["solver_settings"])

# Add problem variables
solver.AddVariables()

# Read model_part (note: the buffer_size is set here)
solver.ImportModelPart()

# Add degrees of freedom
solver.AddDofs()

# Initialize GiD I/O
from gid_output_process import GiDOutputProcess
gid_output = GiDOutputProcess(main_model_part,
                              ProjectParameters["problem_data"]["problem_name"].GetString(),
                              ProjectParameters["output_configuration"])
gid_output.ExecuteInitialize()

# Create the Kratos model
ShallowModel = KratosMultiphysics.Model()
ShallowModel.AddModelPart(main_model_part)

# Print model_part and properties
if(echo_level > 1):
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)


## Initialize ------------------------------------------------------------------------------------------------

# Construct processes to be applied
import process_factory
# "list_of_processes" contains all the processes already constructed (boundary conditions, initial conditions and gravity)
# Note 1: bathymetry is firstly constructed. Initial conditions might need its information.
# Note 2: initial conditions are constructed before BCs. Otherwise, they may overwrite the BCs information.
list_of_processes  = process_factory.KratosProcessFactory(ShallowModel).ConstructListOfProcesses( ProjectParameters["bathymetry_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(ShallowModel).ConstructListOfProcesses( ProjectParameters["initial_conditions_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(ShallowModel).ConstructListOfProcesses( ProjectParameters["boundary_conditions_process_list"] )

# Print list of constructed processes
if(echo_level>1):
    for process in list_of_processes:
        print(process)

# Initialize processes
for process in list_of_processes:
    process.ExecuteInitialize()

# Initialize the solver
solver.Initialize()

# Gid ExecuteBeforeSolutionLoop
gid_output.ExecuteBeforeSolutionLoop()


# ExecuteBeforeSolutionLoop
for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()


## Temporal loop ---------------------------------------------------------------------------------------------

step = 0

while(time <= end_time):

    # Update temporal variables
    delta_time = solver.ComputeDeltaTime()
    main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, delta_time)
    step += 1
    time += delta_time
    main_model_part.CloneTimeStep(time)

    print("STEP = ", step)
    print("TIME = ", time)

    # Update imposed conditions
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()

    gid_output.ExecuteInitializeSolutionStep()

    if step > 1:
        solver.Solve()

    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    gid_output.ExecuteFinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()

    # Write GiD results
    if gid_output.IsOutputStep():
        gid_output.PrintOutput()

    for process in list_of_processes:
        process.ExecuteAfterOutputStep()

## Finalize --------------------------------------------------------------------------------------------------

for process in list_of_processes:
    process.ExecuteFinalize()

# Finalizing output files
gid_output.ExecuteFinalize()
