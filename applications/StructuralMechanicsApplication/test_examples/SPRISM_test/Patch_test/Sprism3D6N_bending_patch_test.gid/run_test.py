from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# import the configuration data as read from the GiD
import ProjectParameters

from numpy import *

import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

# time control starts
from time import *
print(ctime())
# measure process time
t0p=clock()

def StartTimeMeasuring():
    # measure process time
    time_ip=clock()
    return time_ip

def StopTimeMeasuring(time_ip, process):
    # measure process time
    time_fp=clock()
    print(" ", process, " [ spent time=", time_fp - time_ip, "] ")

def SetParallelSize(num_threads):
    parallel=OpenMPUtils()
    print("Num Threads=", num_threads)
    parallel.SetNumThreads(int(num_threads))

# defining the number of threads:
num_threads=ProjectParameters.NumberofThreads
SetParallelSize(num_threads)

# defining the model parts
model_part = ModelPart("StructurePart")
model_part.AddNodalSolutionStepVariable(ALPHA_EAS)

# importing the solvers needed
SolverSettings=ProjectParameters.SolidSolverConfiguration
solver_module=__import__(SolverSettings.solver_type)

# setting the domain size for the problem to be solved
domain_size=SolverSettings.domain_size

# importing variables
solver_module.AddVariables(model_part, SolverSettings)

# defining the type, the name and the path of the problem:
input_file_name=ProjectParameters.problem_name
problem_path=ProjectParameters.problem_path

# reading the solid part
model_part_io=ModelPartIO(input_file_name)
model_part_io.ReadModelPart(model_part)
solver_module.AddDofs(model_part, ProjectParameters)

# Creating the solid solver, set the constitutive law
import constitutive_law_python_utility as constitutive_law_utils
constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(model_part, domain_size);
constitutive_law.Initialize();

# --- PRINT CONTROL ---#
print(model_part)
print(model_part.Properties[1])

# Find neighbours if required
sprism_neighbour_search = SprismNeighbours(model_part)
sprism_neighbour_search.Execute()

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
model_part.SetBufferSize(3)

# importing the solver files
solid_solver = solver_module.CreateSolver(model_part, SolverSettings)

solid_solver.Initialize()
(solid_solver).SetEchoLevel(SolverSettings.echo_level);
solid_solver.Initialize()

# initialize GiD  I/O
#gid_mode=GiDPostMode.GiD_PostBinary
gid_mode=GiDPostMode.GiD_PostAscii
multifile=MultiFileFlag.MultipleFiles
deformed_mesh_flag=WriteDeformedMeshFlag.WriteUndeformed
write_conditions=WriteConditionsFlag.WriteElementsOnly
gid_io=GidIO(input_file_name, gid_mode, multifile, deformed_mesh_flag, write_conditions)

mesh_name=0.0
gid_io.InitializeMesh(mesh_name)
gid_io.WriteMesh((model_part).GetMesh())
gid_io.FinalizeMesh()

# Stepping and time settings
Start_time=ProjectParameters.Start_time
Dt=ProjectParameters.Dt
Nsteps=ProjectParameters.nsteps
final_time=ProjectParameters.max_time
output_time=ProjectParameters.output_time

out=0
step=0
time=Start_time

gid_io.InitializeResults(mesh_name, (model_part).GetMesh())

for node in model_part.Nodes:
    if (node.X >2.40000e-01 -1.0e-5) | (node.Y > 1.20000e-01 -1.0e-5) | (node.X < 1.0e-5) | (node.Y < 1.0e-5):
        node.Fix(DISPLACEMENT_X)
        node.SetSolutionStepValue(DISPLACEMENT_X, - 1.0e-7 *(node.Z - 0.0005)*(node.X+node.Y/2))
        node.Fix(DISPLACEMENT_Y)
        node.SetSolutionStepValue(DISPLACEMENT_Y, - 1.0e-7 *(node.Z - 0.0005)*(node.Y+node.X/2))
        node.Fix(DISPLACEMENT_Z)
        node.SetSolutionStepValue(DISPLACEMENT_Z, 0.5 * 1.0e-7 *(node.X**2+node.X*node.Y+node.Y**2))

while(time <= final_time): 
      
    model_part.CloneTimeStep(time)

    time+=Dt
    step+=1
    out += Dt
    
    print("STEP=", step)
    print("TIME=", time)
    
    if step>=2:
      
      solid_solver.Solve()

      if(output_time <= out):
        gid_io.WriteNodalResults(DISPLACEMENT, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(REACTION, model_part.Nodes, time, 0)
        gid_io.WriteNodalResults(POINT_LOAD, model_part.Nodes, time, 0)
        gid_io.PrintOnGaussPoints(PK2_STRESS_TENSOR, model_part, time)
        gid_io.PrintOnGaussPoints(CAUCHY_STRESS_TENSOR, model_part, time)
        gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_TENSOR, model_part, time)
        gid_io.PrintOnGaussPoints(HENCKY_STRAIN_TENSOR, model_part, time)

gid_io.FinalizeResults()

print("Analysis Finalized ")

# measure process time
tfp=clock()
print(ctime())
print("Analysis Completed  [Process Time=", tfp - t0p, "] ")