import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
import csv, os

try:
    import nlopt
    nlopt_available = True
except ImportError:
    nlopt_available = False

with open("optimization_parameters.json", "r") as file_input:
    parameters = Kratos.Parameters(file_input.read())
model = Kratos.Model()
analysis = OptimizationAnalysis(model, parameters)
analysis.Run()
