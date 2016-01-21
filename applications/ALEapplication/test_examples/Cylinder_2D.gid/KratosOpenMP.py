

# import the configuration data as read from the GiD
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import ProjectParameters
import define_output
import math



#Definition of functions
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Move nodes with translation

def DefineInterface(model_part):
  for node in model_part.GetNodes(6):
    node.Set(INTERFACE,True)
    

def MoveNodes(model_part,time):
  for node in model_part.Nodes:
    if(node.Is(INTERFACE)):
      disp = Vector(3)
      disp[0] = 0.0
      disp[1] = 0.05*math.sin((time-0.02)*100)
      disp[2] = 0.0
      node.SetSolutionStepValue(DISPLACEMENT,0,disp)
      node.Fix(DISPLACEMENT_X)
      node.Fix(DISPLACEMENT_Y)
      node.Fix(DISPLACEMENT_Z)


#Apply displacement boundary conditions

def ApplyDisplacementConditions(model_part):
  for node in model_part.Nodes:
    if(node.X < 0.0001 or node.X > 0.999):
      node.Fix(DISPLACEMENT_X)
      node.Fix(DISPLACEMENT_Z)
    if(node.Y > 0.4999 or node.Y < 0.0001):
      node.Fix(DISPLACEMENT_X)
      node.Fix(DISPLACEMENT_Y)
      node.Fix(DISPLACEMENT_Z)
    if(node.Is(INTERFACE)):
      node.Fix(DISPLACEMENT_X)
      node.Fix(DISPLACEMENT_Y)
      node.Fix(DISPLACEMENT_Z)

# --SET NUMBER OF THREADS --#################


def SetParallelSize(num_threads):
    parallel = OpenMPUtils()
    print("Num Threads = ", num_threads)
    parallel.SetNumThreads(int(num_threads))
# --SET NUMBER OF THREADS --#################

	  

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# setting the domain size for the problem to be solved
domain_size = ProjectParameters.domain_size

import sys
sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ALEApplication import *



#defining variables to be used
# GID IO IS NOT USING THIS NOW. TO BE REMOVED ONCE THE "PRINT IN POINTS"
# CODE IS NOT USING IT

variables_dictionary = {"PRESSURE": PRESSURE,
                        "VELOCITY": VELOCITY,
                        "DISPLACEMENT": DISPLACEMENT,
                        "MESH_VELOCITY": MESH_VELOCITY,
                        }

# defining a model part for the fluid
fluid_model_part = ModelPart("FluidPart")

if "MESH_VELOCITY" in ProjectParameters.nodal_results:
    fluid_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)



# Importing the solvers needed
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

SolverType = ProjectParameters.SolverType		#Parameters coming from ProjectParameters.mdpa

num_threads = ProjectParameters.NumberofThreads
SetParallelSize(num_threads)
MeshSolverType = ProjectParameters.MeshSolverType	#Parameters coming from ProjectParameters.mdpa


#Fluidsolvers

if(SolverType == "FractionalStep"):
    import vms_fractional_step_solver as solver
    solver.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian"):
    import vms_monolithic_solver as solver
    solver.AddVariables(fluid_model_part)
elif(SolverType == "monolithic_solver_eulerian_compressible"):
    import monolithic_solver_eulerian_compressible as solver
    solver.AddVariables(fluid_model_part)
else:
    raise NameError("solver type not supported: options are FractionalStep  - monolithic_solver_eulerian - monolithic_solver_eulerian_compressible")


SolverSettings = ProjectParameters.FluidSolverConfiguration
solver_module = import_solver(SolverSettings)



#Meshsolvers

if(MeshSolverType == "Laplacian"):
    import mesh_solver_laplacian as mesh_solver
    mesh_solver.AddVariables(fluid_model_part)
elif(MeshSolverType == "StructuralSimilarity"):
    import mesh_solver_structural_similarity as mesh_solver
    mesh_solver.AddVariables(fluid_model_part)
elif(MeshSolverType == "StructuralSimilarityNonlinear"):
    import mesh_solver_structural_similarity_nonlinear as mesh_solver
    mesh_solver.AddVariables(fluid_model_part)
elif(MeshSolverType == "LaplacianComponentwise"):
    import mesh_solver_laplacian_componentwise as mesh_solver
    mesh_solver.AddVariables(fluid_model_part)
else:
    raise NameError("solver type not supported: options are LaplacianComponentwise  - StructuralSimilarity - StructuralSimilarityNonlinear")
    

    
    
#Imports   
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Importing variables
solver_module.AddVariables(fluid_model_part, SolverSettings)

# introducing input file name
input_file_name = ProjectParameters.problem_name

# reading the fluid part
model_part_io_fluid = ModelPartIO(input_file_name)
model_part_io_fluid.ReadModelPart(fluid_model_part)

# setting up the buffer size: SHOULD BE DONE AFTER READING!!!
fluid_model_part.SetBufferSize(3)

solver_module.AddDofs(fluid_model_part, SolverSettings)
mesh_solver.AddDofs(fluid_model_part)

# copy Y_WALL
for node in fluid_model_part.Nodes:
    y = node.GetSolutionStepValue(Y_WALL, 0)
    node.SetValue(Y_WALL, y)


# Creating the fluid solver
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print ("creating fluid solver")
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
print ("fluid solver created")


#Creating the mesh solver
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print ("creating mesh solver")
reform_dofs_at_each_step = False
Tolerance = 10e-5
MaxIter = 40

if(MeshSolverType == "Laplacian"):
  mesh_sol = mesh_solver.MeshSolverLaplacian(fluid_model_part,domain_size,reform_dofs_at_each_step)
  mesh_sol.time_order = fluid_solver.time_order
  mesh_sol.Initialize()
elif(MeshSolverType == "StructuralSimilarity"):
  mesh_sol = mesh_solver.MeshSolverStructuralSimilarity(fluid_model_part,reform_dofs_at_each_step)
  mesh_sol.time_order = fluid_solver.time_order
  mesh_sol.Initialize()
elif(MeshSolverType == "StructuralSimilarityNonlinear"):
  mesh_sol = mesh_solver.MeshSolverStructuralSimilarityNonlinear(fluid_model_part,domain_size,reform_dofs_at_each_step,Tolerance,MaxIter)
  mesh_sol.time_order = fluid_solver.time_order
  mesh_sol.Initialize()
elif(MeshSolverType == "LaplacianComponentwise"):
  mesh_sol = mesh_solver.MeshSolverLaplacianComponentwise(fluid_model_part,domain_size,reform_dofs_at_each_step)
  mesh_sol.time_order = fluid_solver.time_order
  mesh_sol.Initialize()
else:
  raise NameError("Selected solver has not been imported !")
  
  
print ("mesh solver created")

# initialize GiD  I/O
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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

    
gid_io.initialize_results(fluid_model_part)



# define the drag computation list
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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



def PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time):
    i = 0
    for it in drag_list:
        print (it[0])
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
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import point_graph_printer

output_nodes_list = define_output.DefineOutputPoints()
graph_printer = point_graph_printer.PrintGraphPrinter(
    output_nodes_list,
    fluid_model_part,
    variables_dictionary,
    domain_size)


# Stepping and time settings
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Dt = ProjectParameters.Dt
Nsteps = ProjectParameters.nsteps
final_time = ProjectParameters.max_time
output_time = ProjectParameters.output_time

time = ProjectParameters.Start_time
out = 0
step = 0


#Apply boundary conditions
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DefineInterface(fluid_model_part)
ApplyDisplacementConditions(fluid_model_part)



#Solution loop
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

while(time <= final_time):

    time = time + Dt
    step = step + 1
    fluid_model_part.CloneTimeStep(time)

    print ("STEP = ", step)
    print ("TIME = ", round(time,6))

    #if(step >= 3 and step < 50):
      #fluid_solver.Solve()
      #graph_printer.PrintGraphs(time)
      #PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)
    if (step >= 3):
      MoveNodes(fluid_model_part,time)
      mesh_sol.Solve()
      mesh_sol.MoveNodes()
      #fluid_solver.Solve()
      graph_printer.PrintGraphs(time)
      PrintDrag(drag_list, drag_file_output_list, fluid_model_part, time)
      
    if(output_time <= out):
      gid_io.write_results(
            time,
            fluid_model_part,
            ProjectParameters.nodal_results,
            ProjectParameters.gauss_points_results)
      out = 0
        
    
    out = out + Dt

gid_io.finalize_results()


for i in drag_file_output_list:
    i.close()

