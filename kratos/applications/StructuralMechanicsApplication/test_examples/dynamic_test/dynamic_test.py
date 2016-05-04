# Import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

# Modelpart for the solid
model_part = ModelPart("StructuralPart")
model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
model_part.AddNodalSolutionStepVariable(VELOCITY)
model_part.AddNodalSolutionStepVariable(ACCELERATION)
model_part.AddNodalSolutionStepVariable(REACTION) 
model_part.AddNodalSolutionStepVariable(NODAL_MASS) 
model_part.AddNodalSolutionStepVariable(NODAL_STIFFNESS) 

# Reading the fluid part
model_part_io = ModelPartIO("dynamic_test")
model_part_io.ReadModelPart(model_part)

# Setting up the buffer size
model_part.SetBufferSize(3)

# Adding dofs
for node in model_part.Nodes:
  node.AddDof(DISPLACEMENT_X, REACTION_X)
  node.AddDof(DISPLACEMENT_Y, REACTION_Y)
  node.AddDof(DISPLACEMENT_Z, REACTION_Z)
  
max_iters = 30

linear_solver = SkylineLUFactorizationSolver()

# Create solver
builder_and_solver = ResidualBasedBuilderAndSolver(linear_solver)
mechanical_scheme = ResidualBasedBossakDisplacementScheme(0.0, 1)
mechanical_convergence_criterion = ResidualCriteria(1e-4, 1e-4)

mechanical_solver = ResidualBasedNewtonRaphsonStrategy(
            model_part,
            mechanical_scheme,
            linear_solver,
            mechanical_convergence_criterion,
            builder_and_solver,
            max_iters,
            False ,
            True,
            True)

mechanical_solver.SetEchoLevel(0)
mechanical_solver.Initialize()

step = 0
time = 0.0

Dt = 1.0e-3
final_time = 1.0

for node in model_part.Nodes:
  node.SetSolutionStepValue(DISPLACEMENT_X, 0, 0.1)

while(time <= final_time):

  time += Dt
  step += 1
  
  model_part.CloneTimeStep(time)
  
  mechanical_solver.Solve()
  
  for node in model_part.Nodes:
    disp = node.GetSolutionStepValue(DISPLACEMENT_X, 0)
  
  print("Time: ",time, "disp: ", disp )
  