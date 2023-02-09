import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.OptimizationApplication.optimization_analysis import (
    OptimizationAnalysis, )

with open("./optimization_parameters.json", "r") as file_input:
    parameters = Kratos.Parameters(file_input.read())

model = Kratos.Model()
analysis = OptimizationAnalysis(model, parameters)
analysis.Run()
