

# import the configuration data as read from the GiD
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import ProjectParameters
import math

#Definition of functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Move nodes with translation

def DefineInterface(model_part):
  for node in model_part.GetNodes(2):
    node.Set(INTERFACE,True)
    

def MoveNodes(model_part,time):
  for node in model_part.Nodes:
    if(node.Is(INTERFACE)):
      disp = Vector(3)
      disp[0] = 0.0    #0.05*math.sin((time-0.02)*100)
      disp[1] = 0.0
      disp[2] = 0.03*math.sin((time-0.02)*100) 
      node.SetSolutionStepValue(MESH_DISPLACEMENT,0,disp)
      node.Fix(MESH_DISPLACEMENT_X)
      node.Fix(MESH_DISPLACEMENT_Y)
      node.Fix(MESH_DISPLACEMENT_Z)


#Apply displacement boundary conditions

def ApplyDisplacementConditions(model_part):
  for node in model_part.Nodes:
    if(node.X < 0.0001 or node.X > 1.999):
      node.Fix(MESH_DISPLACEMENT_X)
      node.Fix(MESH_DISPLACEMENT_Y)
      node.Fix(MESH_DISPLACEMENT_Z)
    if(node.Y > 0.4999 or node.Y < 0.0001):
      node.Fix(MESH_DISPLACEMENT_X)
      node.Fix(MESH_DISPLACEMENT_Y)
      node.Fix(MESH_DISPLACEMENT_Z)
    if(node.Is(INTERFACE)):
      node.Fix(MESH_DISPLACEMENT_X)
      node.Fix(MESH_DISPLACEMENT_Y)
      node.Fix(MESH_DISPLACEMENT_Z)

#Set number of threads

def SetParallelSize(num_threads):
    parallel = OpenMPUtils()
    print("Num Threads = ", num_threads)
    parallel.SetNumThreads(int(num_threads))
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# setting the domain size for the problem to be solved

domain_size = ProjectParameters.domain_size

import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ALEApplication import *



# defining a model part for the fluid
model_part = ModelPart("Modelpart")


# Importing the solvers needed (example fluid simulation)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

num_threads = ProjectParameters.NumberofThreads
SetParallelSize(num_threads)
MeshSolverType = ProjectParameters.MeshSolverType

#Meshsolvers
if(MeshSolverType == "Laplacian"):
    import mesh_solver_laplacian as mesh_solver
    mesh_solver.AddVariables(model_part)
elif(MeshSolverType == "StructuralSimilarity"):
    import mesh_solver_structural_similarity as mesh_solver
    mesh_solver.AddVariables(model_part)
else:
    raise NameError("solver type not supported: options are StructuralSimilarity - Laplacian")
    

    
    
#Imports   
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# introducing input file name
input_file_name = ProjectParameters.problem_name

# reading the fluid part
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
model_part.SetBufferSize(3)

mesh_solver.AddDofs(model_part)




#Creating the mesh solver
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print ("Creating mesh solver")

reform_dofs_at_each_step = False
Tolerance = 10e-5
MaxIter = 40

if(MeshSolverType == "Laplacian"):
  mesh_sol = mesh_solver.MeshSolverLaplacian(model_part,reform_dofs_at_each_step)
  mesh_sol.time_order = 2
  mesh_sol.Initialize()
elif(MeshSolverType == "StructuralSimilarity"):
  mesh_sol = mesh_solver.MeshSolverStructuralSimilarity(model_part,reform_dofs_at_each_step, False)
  mesh_sol.time_order = 2
  mesh_sol.Initialize()
else:
  raise NameError("Selected solver has not been imported !")
  

# initialize GiD  I/O
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Initializing GID I/O")

from gid_output import GiDOutput
gid_io = GiDOutput(input_file_name,
                   ProjectParameters.VolumeOutput,
                   ProjectParameters.GiDPostMode,
                   ProjectParameters.GiDMultiFileFlag,
                   ProjectParameters.GiDWriteMeshFlag,
                   ProjectParameters.GiDWriteConditionsFlag)

if not ProjectParameters.VolumeOutput:
    cut_list = define_output.DefineCutPlanes()
    gid_io.define_cuts(model_part, cut_list)

    
gid_io.initialize_results(model_part)


# Stepping and time settings
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Dt = ProjectParameters.Dt
Nsteps = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0


#Apply boundary conditions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Defining interface and boundary conditions")

DefineInterface(model_part)
ApplyDisplacementConditions(model_part)


#Solution loop
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Solve problem")

while(time <= final_time):

    time = time + Dt
    step = step + 1
    model_part.CloneTimeStep(time)

    print ("STEP = ", step)
    print ("TIME = ", round(time,6))

    if (step >= 3):
      MoveNodes(model_part,time)
      mesh_sol.Solve()
      mesh_sol.MoveNodes()
      
    if(output_time <= out):
      gid_io.write_results(
            time,
            model_part,
            ProjectParameters.nodal_results,
            ProjectParameters.gauss_points_results)
      out = 0
        
    
    out = out + Dt

gid_io.finalize_results()

print("Solution finished")



