from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Time control
import time
print (time.ctime())
start_time = time.clock()


## Necessary modules -----------------------------------------------------------------------------------------

# Including kratos path
from KratosMultiphysics import *
# Including Applications path
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
from KratosMultiphysics.PoromechanicsApplication import *
from KratosMultiphysics.DamApplication import *

#import define_output
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

#Import solver constructors
import dam_mechanical_solver as mechanical_solver
import eulerian_convection_diffusion_solver as diffusion_solver

# Import utilities
import conditions_utility
import constitutive_law_utility
import gid_print_utility
import cleaning_utility


## Previous definitions --------------------------------------------------------------------------------------

# Number of threads
parallel=OpenMPUtils()
parallel.SetNumThreads(int(ProjectParameters["general_data"]["NumberofThreads"].GetInt()))

# Problem parameters
problem_path = os.getcwd()
problem_name = ProjectParameters["general_data"]["problem_name"].GetString()
delta_time = ProjectParameters["general_data"]["delta_time"].GetDouble()
ending_time = ProjectParameters["general_data"]["ending_time"].GetDouble()
time_converter = ProjectParameters["general_data"]["time_scale"].GetString()
evolution_type = ProjectParameters["general_data"]["evolution_type"].GetString()
current_step = 0
current_time = 0.0
current_id = 1
tol = delta_time*1.0e-10

# Time Units Converter
if(time_converter=="Months"):               # Factor to pass from months to seconds
    time_unit_converter = 2592000.0
elif(time_converter=="Days"):               # Factor to pass from days to seconds
    time_unit_converter = 86400.0
elif(time_converter=="Hours"):              # Factor to pass from hours to seconds
    time_unit_converter = 3600.0
else:
    time_unit_converter = 1.0               # No changes

# Update time variables
delta_time = delta_time * time_unit_converter
ending_time = ending_time * time_unit_converter

# List of variables to write
nodal_res = ProjectParameters["output_configuration"]["result_file_configuration"]["nodal_results"]
gp_res = ProjectParameters["output_configuration"]["result_file_configuration"]["gauss_point_results"]

## Model part ------------------------------------------------------------------------------------------------

# Definition of model part
model_part = ModelPart("SolidDomain")

# Setting thermal variables
diffusion_solver.AddVariables(model_part,ProjectParameters["diffusion_settings"]) 

# Set mechanical variables
mechanical_solver.AddVariables(model_part)

# Reading model part
model_part_io = ModelPartIO(problem_name)
model_part_io.ReadModelPart(model_part)

# Set buffer size
buffer_size = 2
model_part.SetBufferSize(buffer_size)

# Set thermal degrees of freedom
# TODO: fix the problems with convection-diffusion solver
#diffusion_solver.AddDofs(model_part)  

# Set mechanical degrees of freedom
mechanical_solver.AddDofs(model_part)

# Set ProcessInfo variables and fill the previous steps of the buffer with the initial conditions
current_time = -(buffer_size-1)*delta_time
model_part.ProcessInfo[TIME] = current_time #current_time and TIME = 0 after filling the buffer

#Type of reference temperature
reference_temperature = ProjectParameters["diffusion_settings"]["reference_temperature"]
if (reference_temperature =="Reservoir Information"):
    model_part.ProcessInfo[REFERENCE_TEMPERATURE] = 10.0  #TODO: for working
#model_part.ProcessInfo[REFERENCE_TEMPERATURE] = model_part.GetTable(18).GetNearestValue(0.0)      To start computations at time = 0  / Table 5 = Reference Temperature Values
else:
    model_part.ProcessInfo[REFERENCE_TEMPERATURE] = 10.0 # TODO: Here we have to solve the thermal problem until arrives a stationary    

for step in range(buffer_size-1):
    current_time = current_time + delta_time
    model_part.CloneTimeStep(current_time)


## Initialize ------------------------------------------------------------------------------------------------

# Definition of utilities
# TODO: Problems with inputs in gid_output_util, now we dont have Conditions options (think about it) 
#gid_output_util = gid_print_utility.GidPrintUtility(ProjectParameters["output_configuration"],problem_name,current_time,ending_time,delta_time,time_unit_converter)
cleaning_util = cleaning_utility.CleaningUtility(problem_path)
#conditions_util = conditions_utility.ConditionsUtility(delta_time,ProjectParameters.ConditionsOptions, model_part, time_unit_converter,evolution_type)

# Erasing previous results files
cleaning_util.CleanPreviousFiles()

# Set constitutive laws
constitutive_law_utility.SetConstitutiveLaw(model_part)

# Define and initialize the diffusion solver
# TODO: Problems with the solver, domain_size is not found. 
#thermal_diffusion_solver = diffusion_solver.CreateSolver(model_part, ProjectParameters)

#Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
thermal_scheme = ProjectParameters["diffusion_settings"]["temporal_scheme"]
if(thermal_scheme=="Backward-Euler"):
    model_part.ProcessInfo[THETA] = 1.0
elif(thermal_scheme=="Crank-Nicolson"):   
    model_part.ProcessInfo[THETA] = 0.5
else:
    model_part.ProcessInfo[THETA] = 0.0
    
thermal_diffusion_solver.Initialize()
thermal_diffusion_solver.SetEchoLevel(0)

# Define and initialize the mechanical solver
solid_mechanics_solver = mechanical_solver.CreateSolver(model_part, ProjectParameters.MechanicalSolverSettings)
solid_mechanics_solver.Initialize()

# Initialize imposed conditions
conditions_util.Initialize(model_part)

# Initializing new results
gid_output_util.initialize_results(model_part, current_id) #For single post file
gid_output_util.write_results(model_part, nodal_res, gp_res, current_time, current_step, current_id)


## Temporal loop ---------------------------------------------------------------------------------------------

while( (current_time+tol) < ending_time ):

    # Update temporal variables
    current_time = current_time + delta_time
    current_step = current_step + 1
    model_part.CloneTimeStep(current_time)
    print("--------------------------------------------------")
    print("STEP",current_step," - TIME","%.5f" % current_time)
    
    # Update imposed conditions
    conditions_util.UpdateImposedConditions(model_part,current_step)

    # Solve thermal step
    clock_time = time.clock()
    thermal_diffusion_solver.Solve()
    print("Thermal Solving Time = ","%.5f" % (time.clock() - clock_time)," seconds")
    
    # Solve mechanical step
    clock_time = time.clock()
    solid_mechanics_solver.Solve()
    print("Mechanical Solving Time = ","%.5f" % (time.clock() - clock_time)," seconds")

    # Write GiD results
    execute_write = gid_output_util.CheckWriteResults(current_time)
    if(execute_write):
        current_id = current_id + 1
        clock_time = time.clock()
        gid_output_util.write_results(model_part, nodal_res, gp_res, current_time, current_step, current_id)
        print("Writing Time = ","%.5f" % (time.clock() - clock_time)," seconds")


## Finalize --------------------------------------------------------------------------------------------------

# Finalizing mechanical strategy
solid_mechanics_solver.Finalize()

# Finalizing output files
gid_output_util.finalize_results()

# Time control
print (time.ctime())
print("Analysis Completed, Elapsed Time = ", time.clock() - start_time," seconds")
