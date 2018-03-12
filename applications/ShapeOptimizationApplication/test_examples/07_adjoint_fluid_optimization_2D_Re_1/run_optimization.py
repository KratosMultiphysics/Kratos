from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# For time measures
import time as timer

# ======================================================================================================================================
# Model part & solver
# ======================================================================================================================================

#import define_output
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

#set echo level
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()

#defining the model_part
main__optimization_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main__optimization_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available in kratos core
Model_optimization = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main__optimization_model_part}

adjoint_optimization = __import__("adjoint_optimization")
adjointAnalyzer = adjoint_optimization.kratosAdjointFluidAnalyzer( main__optimization_model_part, ProjectParameters )

# Create an optimizer 
# Note that internally variables related to the optimizer are added to the model part
optimizerFactory = __import__("optimizer_factory")
optimizer = optimizerFactory.CreateOptimizer( main__optimization_model_part, ProjectParameters["optimization_settings"] )

# Create solver for handling mesh-motion
mesh_solver_module = __import__(ProjectParameters["mesh_solver_settings"]["solver_type"].GetString())
mesh_solver = mesh_solver_module.CreateSolver(main__optimization_model_part, ProjectParameters["mesh_solver_settings"])
mesh_solver.AddVariables()


# ======================================================================================================================================
# Analyzer
# ======================================================================================================================================

adjointAnalyzer.SetMeshSolver(mesh_solver)
optimizer.importModelPart()
mesh_solver.AddDofs()

# Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
## Get the list of the submodel part in the object Model
for i in range(ProjectParameters["adjoint_fluid_analyzer"]["sub_model_part_list"].size()):
    part_name = ProjectParameters["adjoint_fluid_analyzer"]["sub_model_part_list"][i].GetString()
    if( main__optimization_model_part.HasSubModelPart(part_name) ):
        Model_optimization.update({part_name: main__optimization_model_part.GetSubModelPart(part_name)})

adjointAnalyzer.Initialize(Model_optimization)

# ======================================================================================================================================
# Optimization
# ======================================================================================================================================

optimizer.importAnalyzer( adjointAnalyzer )
adjointAnalyzer.importOptimizer( optimizer )
optimizer.optimize()
adjointAnalyzer.finalizeSolutionLoop()

# ======================================================================================================================================
