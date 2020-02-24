from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication.structural_mechanics_analysis_rom import StructuralMechanicsAnalysisROM

import numpy as np
import json

class TestStructuralMechanicsDynamicROM(StructuralMechanicsAnalysisROM):

    def __init__(self,model,project_parameters):
        super(TestStructuralMechanicsDynamicROM,self).__init__(model,project_parameters)
        self.selected_time_step_solution_container = []
    
    def FinalizeSolutionStep(self):
        super(TestStructuralMechanicsDynamicROM,self).FinalizeSolutionStep()        
        ArrayOfDisplacements = []
        if np.isclose(self.time,2) or np.isclose(self.time,4) or np.isclose(self.time,6) or np.isclose(self.time,8) or np.isclose(self.time,10):
            for node in self._GetSolver().GetComputingModelPart().Nodes:            
                ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0))
                ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0))
            self.selected_time_step_solution_container.append(ArrayOfDisplacements)

    def EvaluateQuantityOfInterest(self):
       ##############################################################################################
       # Functions evaluating the QoI of the problem: Array of displacement at every node on the mesh #
       #                    and nodal area                                                           #
       ##############################################################################################
        SnapshotMatrix = np.zeros((len(self.selected_time_step_solution_container[0]), len(self.selected_time_step_solution_container)))
        for i in range(len(self.selected_time_step_solution_container)):
            Snapshot_i= np.array(self.selected_time_step_solution_container[i])
            SnapshotMatrix[:,i] = Snapshot_i.transpose()
        return SnapshotMatrix       

    def EvaluateQuantityOfInterest2(self):
        computing_model_part = self._solver.GetComputingModelPart()
        dimension = self._GetSolver().settings["domain_size"].GetInt()
        area_calculator = KratosMultiphysics.CalculateNodalAreaProcess(computing_model_part , dimension)
        area_calculator.Execute()        

        ArrayOfAreas = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            ArrayOfAreas.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))
            ArrayOfAreas.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA))            
        return ArrayOfAreas

if __name__ == "__main__":    
    with open("ProjectParametersROM.json",'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())
    model = KratosMultiphysics.Model()      
    Simulation = TestStructuralMechanicsDynamicROM(model,parameters)
    Simulation.Run()