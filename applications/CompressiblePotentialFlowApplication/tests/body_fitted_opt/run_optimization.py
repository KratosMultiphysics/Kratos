import KratosMultiphysics as km
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.ShapeOptimizationApplication as kso
import KratosMultiphysics.ShapeOptimizationApplication.optimizer_factory as optimizer_factory
km.Logger.GetDefaultOutput().SetSeverity(km.Logger.Severity.WARNING)

with open("optimization_parameters.json",'r') as parameter_file:
    parameters_optimization = km.Parameters(parameter_file.read())

# Create optimizer and perform optimization
model = km.Model()
optimizer = optimizer_factory.CreateOptimizer(parameters_optimization["optimization_settings"], model)
optimizer.Optimize()
