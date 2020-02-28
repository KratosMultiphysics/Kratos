from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication.convection_diffusion_analysis_rom import ConvectionDiffusionAnalysisROM

import numpy as np

class TestConvectionDiffusionTransientROM(ConvectionDiffusionAnalysisROM):

    def __init__(self,model,project_parameters):        
        super(TestConvectionDiffusionTransientROM,self).__init__(model,project_parameters)
        self.selected_time_step_solution_container = []
    
    def FinalizeSolutionStep(self):
        super(TestConvectionDiffusionTransientROM,self).FinalizeSolutionStep()        
        ArrayOfTemperatures = []        
        if self.time==500 or self.time==1200 or self.time==2500 or self.time==3000 or self.time==3600:
            for node in self._solver.GetComputingModelPart().Nodes:       
                ArrayOfTemperatures.append(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
            self.selected_time_step_solution_container.append(ArrayOfTemperatures)  

    def EvaluateQuantityOfInterest(self):
       ################################################################################################
       # Functions evaluating the QoI of the problem: Numpy array with selected snapshots of solution #
       #                                    and nodal area                                            #
       ################################################################################################
        SnapshotMatrix = np.zeros((len(self.selected_time_step_solution_container[0]), len(self.selected_time_step_solution_container)))
        for i in range(len(self.selected_time_step_solution_container)):
            Snapshot_i= np.array(self.selected_time_step_solution_container[i])
            SnapshotMatrix[:,i] = Snapshot_i.transpose()
        return SnapshotMatrix

    def EvaluateQuantityOfInterest2(self):
        computing_model_part = self.model["ThermalModelPart"]
        dimension = self._GetSolver().settings["domain_size"].GetInt()
        area_calculator = KratosMultiphysics.CalculateNodalAreaProcess(computing_model_part , dimension)
        area_calculator.Execute()        

        ArrayOfAreas = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            ArrayOfAreas.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
        return ArrayOfAreas

############################################################################################################

if __name__ == "__main__":

    with open("ProjectParametersROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = TestConvectionDiffusionTransientROM(model,parameters)
    simulation.Run()
