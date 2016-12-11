

# import the configuration data as read from the GiD
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

from __future__ import print_function, absolute_import, division

import ProjectParameters
from math import sin,cos

#Definition of functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

movementList = ["TRANSLATION_X", "TRANSLATION_Y", "BENDING", "ROTATION"]
# =============================================================
movement = 0             # set 0-3
amplificationFactor = 1  # factor to amplify the movement
solveMesh = True         # set to "True" to solve the mesh
# =============================================================

#Currently not needed but still defined

def DefineInterface(model_part):
  for node in model_part.GetNodes(2):
    node.Set(INTERFACE,True)
    

#Apply displacement boundary conditions

def ApplyDisplacementConditions(model_part):
  for node in model_part.Nodes:
    if(node.X < 0.0001 or node.X > 1.999):
      node.Fix(MESH_DISPLACEMENT_X)
      node.Fix(MESH_DISPLACEMENT_Z)
    if(node.Y > 0.4999 or node.Y < 0.0001):
      node.Fix(MESH_DISPLACEMENT_X)
      node.Fix(MESH_DISPLACEMENT_Y)
      node.Fix(MESH_DISPLACEMENT_Z)
    if(node.Is(INTERFACE)):
      node.Fix(MESH_DISPLACEMENT_X)
      node.Fix(MESH_DISPLACEMENT_Y)
      node.Fix(MESH_DISPLACEMENT_Z)
      
def DisplacementToMesh(fluid_interface, time, movement, ampFac):
    time *= ampFac
    for node in fluid_interface:
        if movement == "TRANSLATION_X":
            valueX = 0.5 * time
            valueY = 0
        elif movement == "TRANSLATION_Y":
            valueX = 0
            valueY = 0.2 * time
        elif movement == "BENDING":
            valueX = 0
            valueY = 3 * node.X0 * node.X0 * time
        elif movement == "ROTATION":
            xOld = node.X0
            yOld = node.Y0
            valueX = (cos(time*.6)*xOld + sin(time*.6)*yOld) - xOld
            valueY = (-sin(time*.6)*xOld + cos(time*.6)*yOld) - yOld
        else:
            wait = input("Wrong type of movement specified, please correct input")
        # set the prescribed values
        node.SetSolutionStepValue(MESH_DISPLACEMENT_X,0,valueX)
        node.SetSolutionStepValue(MESH_DISPLACEMENT_Y,0,valueY)
      

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
  mesh_sol = mesh_solver.MeshSolverLaplacian(model_part,domain_size,reform_dofs_at_each_step)
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
move = movementList[movement]

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
      DisplacementToMesh(model_part.GetNodes(2), time, move, amplificationFactor)
      if solveMesh == True:
        mesh_sol.Solve()
      gid_io.write_results(
            time,
            model_part,
            ProjectParameters.nodal_results,
            ProjectParameters.gauss_points_results)

gid_io.finalize_results()

print("Solution finished")



