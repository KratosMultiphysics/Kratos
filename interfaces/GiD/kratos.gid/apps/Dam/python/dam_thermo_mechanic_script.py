from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Time monitoring
from time import *
print (ctime())
start_time = clock()

## Necessary modules -----------------------------------------------------------------------------------------

# Including kratos path
from KratosMultiphysics import *
# Including Applications path
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
from KratosMultiphysics.PoromechanicsApplication import *
from KratosMultiphysics.DamApplication import *

##############################################################################
############################# PARSING PARAMETERS #############################
##############################################################################

#import define_output
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

########################## MODEL PART DEFINITIONS ########################## 

main_model_part = ModelPart(ProjectParameters["general_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["general_data"]["domain_size"].GetInt())

Model = {ProjectParameters["general_data"]["model_part_name"].GetString() : main_model_part}

########################## IMPORTING THE SOLVER ########################## 

# Thermal Solver
thermal_solver_module = __import__("dam_eulerian_convection_diffusion_solver")
thermal_solver = thermal_solver_module.CreateSolver(main_model_part, ProjectParameters)

# Mechanical Solver
solver_module = __import__("dam_new_mechanical_solver")  # Up to now is the unique one
solver = solver_module.CreateSolver(main_model_part, ProjectParameters)

# Add variables
thermal_solver.AddVariables()
solver.AddVariables()

# Import the model part
solver.ImportModelPart()

# Add degrees of freedom
thermal_solver.AddDofs()
solver.AddDofs()

########################## IMPORTING NECESSARY MODULES ########################## 

# Import process modules
from gid_output_process import GiDOutputProcess
import cleaning_utility
import process_factory

############################### PREVIOUS DEFINITIONS ############################# 

# Number of threads
parallel=OpenMPUtils()
parallel.SetNumThreads(int(ProjectParameters["general_data"]["NumberofThreads"].GetInt()))

# Problem parameters
problem_path = os.getcwd()
problem_name = ProjectParameters["general_data"]["problem_name"].GetString()
domain_size = ProjectParameters["general_data"]["domain_size"].GetInt()
type_of_problem = ProjectParameters["general_data"]["type_of_problem"].GetString()
delta_time = ProjectParameters["general_data"]["delta_time"].GetDouble()
ending_time = ProjectParameters["general_data"]["ending_time"].GetDouble()
time_converter = ProjectParameters["general_data"]["time_scale"].GetString()
current_step = 0
current_time = 0.0
current_id = 1
tol = delta_time*1.0e-10
echo_level = ProjectParameters["mechanical_settings"]["echo_level"].GetInt()
buffer_size = ProjectParameters["mechanical_settings"]["buffer_size"].GetInt()

# Time Units Converter
if(time_converter=="Months"):               # Factor to pass from months to seconds
    time_unit_converter = 2592000.0
elif(time_converter=="Days"):               # Factor to pass from days to seconds
    time_unit_converter = 86400.0
elif(time_converter=="Hours"):              # Factor to pass from hours to seconds
    time_unit_converter = 3600.0
else:                                       # No changes
    time_unit_converter = 1.0               

main_model_part.ProcessInfo[TIME_UNIT_CONVERTER] = time_unit_converter

#Thermal Parameters
if (type_of_problem == "Thermo-Mechanical"):
    
    #Type of reference temperature
    main_model_part.ProcessInfo[REFERENCE_TEMPERATURE] = ProjectParameters["diffusion_settings"]["reference_temperature"].GetDouble()
    
    #if (reference_temperature =="Reservoir_Information"):
        #main_model_part.ProcessInfo[REFERENCE_TEMPERATURE] = 10.0  #TODO: for working
        ##model_part.ProcessInfo[REFERENCE_TEMPERATURE] = model_part.GetTable(18).GetNearestValue(0.0)      To start computations at time = 0  / Table 5 = Reference Temperature Values
    #else:
        #main_model_part.ProcessInfo[REFERENCE_TEMPERATURE] = 10.0 # TODO: Here we have to solve the thermal problem until arrives a stationary    
    
    #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
    thermal_scheme = ProjectParameters["diffusion_settings"]["temporal_scheme"].GetString()
    
    if(thermal_scheme=="Backward-Euler"):
        main_model_part.ProcessInfo[THETA] = 1.0
    elif(thermal_scheme=="Crank-Nicolson"):   
        main_model_part.ProcessInfo[THETA] = 0.5
    else:
        main_model_part.ProcessInfo[THETA] = 0.0
    
# Update time variables
delta_time = delta_time * time_unit_converter
ending_time = ending_time * time_unit_converter


# Build sub_model_parts
for i in range(ProjectParameters["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["processes_sub_model_part_list"][i].GetString()
    Model.update({part_name : main_model_part.GetSubModelPart(part_name)})


# Print control
if(echo_level > 1):
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)


########################## INITIALIZE ########################## 

# Construct the processes to be applied
list_of_processes = process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["nodal_processes_sub_model_part_list"] )
list_of_processes += process_factory.KratosProcessFactory(Model).ConstructListOfProcesses( ProjectParameters["load_processes_sub_model_part_list"] )


# Print list of constructed processes
if(echo_level>1):
    for process in list_of_processes:
        print(process)


# Initialize processes
for process in list_of_processes:
    process.ExecuteInitialize()


# Set ProcessInfo variables and fill the previous steps of the buffer with the initial conditions
current_time = current_time-(buffer_size-1)*delta_time
main_model_part.ProcessInfo[TIME] = current_time #current_time and TIME = 0 after filling the buffer


for step in range(buffer_size-1):
    current_time = current_time + delta_time
    main_model_part.CloneTimeStep(current_time)


# Initialize GiD I/O
computing_model_part = solver.GetComputingModelPart()
gid_output = GiDOutputProcess(computing_model_part,problem_name,ProjectParameters["output_configuration"])
gid_output.ExecuteInitialize()

#Initializa the mechanical solver
solver.Initialize()
thermal_solver.Initialize()


#Execute before solution
for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()


# Set the results 
gid_output.ExecuteBeforeSolutionLoop()

################################# TEMPORAL LOOP ####3########################## 

while( (current_time+tol) < ending_time ):
    
    # Update temporal variables
    delta_time = main_model_part.ProcessInfo[DELTA_TIME]
    current_time = current_time + delta_time
    current_step = current_step + 1
    main_model_part.CloneTimeStep(current_time)
    print("--------------------------------------------------")
    print("STEP",current_step," - TIME","%.5f" % current_time)
    
    for process in list_of_processes:
        process.ExecuteInitializeSolutionStep()

    gid_output.ExecuteInitializeSolutionStep()
    
    
    # Solve thermal step
    clock_time = clock()
    thermal_solver.Solve()
    print("Thermal Solving Time = ","%.5f" % (clock() - clock_time)," seconds")
        
    # Solve mechanical step
    clock_time = clock()
    solver.Solve()
    print("Mechanical Solving Time = ","%.5f" % (clock() - clock_time)," seconds")


    gid_output.ExecuteFinalizeSolutionStep()
    
    
    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()
    
    
    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()
    
    
    # Write GiD results
    if gid_output.IsOutputStep():
        gid_output.PrintOutput()
    
    
    for process in list_of_processes:
        process.ExecuteAfterOutputStep()


############################### FINALIZE ############################### 

# Finalizing output files
gid_output.ExecuteFinalize()

for process in list_of_processes:
    process.ExecuteFinalize()
        
# Finalizing strategy
solver.Clear()

# Time control
print (ctime())
print("Analysis Completed, Elapsed Time = ", clock() - start_time," seconds")
