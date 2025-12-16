import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.surrogate_modelling import FsiSurrogate


with open("optimization_parameters.json", "r") as file_input:
    parameters = Kratos.Parameters(file_input.read())
models = {'fluid': Kratos.Model(), 'structure': Kratos.Model()}
model = Kratos.Model()
analysis = FsiSurrogate(models, parameters)
analysis.Run()
