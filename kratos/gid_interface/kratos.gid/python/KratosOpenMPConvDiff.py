from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#
#
import ProjectParameters
import define_output
#
#
# setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size
#
#
import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

# defining a model part for the fluid
model_part = ModelPart("FluidPart");  

#
#
# importing the solvers needed
if(ProjectParameters.Linear==True):
    SolverSettings = ProjectParameters.SolverSettings1
else:
    SolverSettings = ProjectParameters.SolverSettings2

solver_module = import_solver(SolverSettings)
#
#
# importing variables
solver_module.AddVariables(model_part, SolverSettings)

#introducing input file name
input_file_name = ProjectParameters.problem_name


# reading the fluid part
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
model_part.SetBufferSize(3)


#define dofs to be stored
solver_module.AddDofs(model_part,SolverSettings)
#
#
# Creating the fluid solver
conv_diff_solver = solver_module.CreateSolver( model_part, SolverSettings)
conv_diff_solver.Initialize()

print("conv_diff solver created")

# initialize GiD  I/O
from gid_output import GiDOutput
gid_io = GiDOutput(input_file_name, ProjectParameters.VolumeOutput, ProjectParameters.GiDPostMode,ProjectParameters.GiDMultiFileFlag,
                   ProjectParameters.GiDWriteMeshFlag,ProjectParameters.GiDWriteConditionsFlag)

if not ProjectParameters.VolumeOutput:
    cut_list = define_output.DefineCutPlanes()
    gid_io.define_cuts(model_part, cut_list)

gid_io.initialize_results(model_part)




#settings to be changed
Dt = ProjectParameters.Dt 
Nsteps  = ProjectParameters.nsteps
full_Dt = Dt 
initial_Dt = 0.01 * full_Dt #0.05 #0.01
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time
output_step = ProjectParameters.output_step
time = ProjectParameters.Start_time
out = 0
step = 0

Stationary=ProjectParameters.Stationary
if(Stationary==True):
    
    model_part.ProcessInfo.SetValue(STATIONARY,1)
    time = 0.01
    model_part.CloneTimeStep(time)	
    	
    conv_diff_solver.Solve()
    gid_io.write_results(time,model_part,ProjectParameters.nodal_results,ProjectParameters.gauss_points_results)	

else:
    model_part.ProcessInfo.SetValue(STATIONARY,0)
    for step in range(0,Nsteps):
        time = Dt*step
        model_part.CloneTimeStep(time)
    
        if(step > 3):
            conv_diff_solver.Solve()
            
	
        if(output_time <= out):
            gid_io.write_results(time,model_part,ProjectParameters.nodal_results,ProjectParameters.gauss_points_results)
            out = 0

        out = out + Dt

    gid_io.finalize_results()    
