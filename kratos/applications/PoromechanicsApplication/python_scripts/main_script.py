from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Time control
import time
print (time.ctime())
start_time = time.clock()


## Necessary modules -----------------------------------------------------------------------------------------

# Including kratos path
from KratosMultiphysics import *
# Including Applications path
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.PoromechanicsApplication import *

# Import parameters
import ProjectParameters

#Import solver constructor
import U_Pw_solver

# Import utilities
import poromechanics_conditions
import poromechanics_constitutivelaw
import poromechanics_gid_output
import poromechanics_cleaning_utility


## Previous definitions --------------------------------------------------------------------------------------

# Number of threads
parallel=OpenMPUtils()
parallel.SetNumThreads(int(ProjectParameters.NumberofThreads))

# Problem parameters
problem_name = os.path.join(str(ProjectParameters.problem_path),str(ProjectParameters.problem_name))
delta_time = ProjectParameters.delta_time
ending_time = ProjectParameters.ending_time
current_step = 0
current_time = 0.0
current_id = 0
tol = delta_time*1e-10

# List of variables to write
nodal_res = ProjectParameters.nodal_results
gp_res = ProjectParameters.gauss_points_results


## Model part ------------------------------------------------------------------------------------------------

# Definition of model part
model_part = ModelPart("PorousDomain")

# Set problem variables
U_Pw_solver.AddVariables(model_part)

# Reading model part 
model_part_io = ModelPartIO(problem_name)
model_part_io.ReadModelPart(model_part)

# Set buffer size
buffer_size = 2
model_part.SetBufferSize(buffer_size)

# Set degrees of freedom
U_Pw_solver.AddDofs(model_part)

# Set TIME and DELTA_TIME and fill the previous steps of the buffer with the initial conditions
current_time = -(buffer_size-1)*delta_time
model_part.ProcessInfo[TIME] = current_time #current_time and TIME = 0 after filling the buffer
for step in range(buffer_size-1):
    current_time = current_time + delta_time
    model_part.CloneTimeStep(current_time)

# Print control
#print(model_part)
#print(model_part.Properties[1])


## Initialize ------------------------------------------------------------------------------------------------

# Definition and initilaization of utilities
gid_output_util = poromechanics_gid_output.PoromechanicsGidOutput(ProjectParameters.GidOutputConfiguration,problem_name,current_time,ending_time,delta_time,nodal_res, gp_res)
cleaning_util = poromechanics_cleaning_utility.PoromechanicsCleaningUtility(ProjectParameters.problem_path)
cleaning_util.CleanPreviousFiles()
conditions_util = poromechanics_conditions.PoromechanicsConditions(delta_time,ProjectParameters.ConditionsOptions)
conditions_util.Initialize(model_part)

# Set constitutive laws
poromechanics_constitutivelaw.SetConstitutiveLaw(model_part)

# Define and initialize the main solver
main_step_solver = U_Pw_solver.CreateSolver(model_part, ProjectParameters.SolverSettings)
main_step_solver.Initialize()

# Initializing new results
gid_output_util.initialize_results(model_part, current_id) #For single post file
print("WRITING RESULTS: ID", current_id, " - STEP", current_step, " - TIME", "%.5f" % current_time)
gid_output_util.write_results(model_part, current_time, current_id)


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
    
    # Solve step
    main_step_solver.Solve()
    
    # Write GiD results
    execute_write = gid_output_util.CheckWriteResults(current_time)
    if(execute_write):
        current_id = current_id + 1
        print("WRITING RESULTS: ID", current_id, " - STEP", current_step, " - TIME", "%.5f" % current_time)
        gid_output_util.write_results(model_part, current_time, current_id)


## Finalize --------------------------------------------------------------------------------------------------

# Finalizing strategy
main_step_solver.Finalize()

# Finalizing output files
gid_output_util.finalize_results()

# Time control
print (time.ctime())
print("Analysis Completed, Elapsed Time = ", time.clock() - start_time," seconds")
