import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication

from KratosMultiphysics.ShallowWaterApplication.shallow_water_analysis import ShallowWaterAnalysis

with open("ProjectParameters.json",'r') as parameter_file:
    parameters = KratosMultiphysics.Parameters(parameter_file.read())

model = KratosMultiphysics.Model()

simulation = ShallowWaterAnalysis(model, parameters)

simulation.Run()
