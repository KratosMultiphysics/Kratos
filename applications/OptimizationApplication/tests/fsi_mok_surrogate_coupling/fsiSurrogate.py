import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.surrogate_modelling import FsiSurrogate


with open("optimization_parameters.json", "r") as file_input:
    parameters = Kratos.Parameters(file_input.read())
models = {'fluid': Kratos.Model(), 'structure': Kratos.Model()}
model = Kratos.Model()
analysis = FsiSurrogate(models, parameters)
analysis.Run()


#2*0.06067*(1-cos(pi*t/10))*y*(1-y)