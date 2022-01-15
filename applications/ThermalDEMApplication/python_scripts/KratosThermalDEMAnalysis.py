from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.ThermalDEMApplication import *
from KratosMultiphysics.ThermalDEMApplication.thermal_dem_analysis import ThermalDEMAnalysis

class ThermalDEMAnalysisWithFlush(ThermalDEMAnalysis):
    def __init__(self,model,parameters):
	    super(ThermalDEMAnalysisWithFlush,self).__init__(model,parameters)

if __name__ == "__main__":
    model = Model()
    with open('ProjectParametersDEM.json','r') as parameter_file:
        parameters = Parameters(parameter_file.read())
    ThermalDEMAnalysisWithFlush(model,parameters).Run()
