from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.RomApplication import structural_mechanics_dynamic_implicit_rom_solver

import numpy as np
import json
"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

class StructDynamicROM(StructuralMechanicsAnalysis):

    def __init__(self,model,project_parameters):
        super(StructDynamicROM,self).__init__(model,project_parameters)
        self.selected_time_step_solution_container = []

    def _CreateSolver(self):        
        return structural_mechanics_dynamic_implicit_rom_solver.ROMSolver(self.model, self.project_parameters["solver_settings"])
    
    def FinalizeSolutionStep(self):
        super(StructDynamicROM,self).FinalizeSolutionStep()        
        ArrayOfDisplacements = []
        if np.isclose(self.time,2) or np.isclose(self.time,4) or np.isclose(self.time,6) or np.isclose(self.time,8) or np.isclose(self.time,10):
            for node in self._GetSolver().GetComputingModelPart().Nodes:            
                ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0))
                ArrayOfDisplacements.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0))
            self.selected_time_step_solution_container.append(ArrayOfDisplacements)

    def ModifyInitialGeometry(self):
        """Here is the place where the BASIS_ROM and the AUX_ID are imposed to each node"""
        super(StructDynamicROM,self).ModifyInitialGeometry()
        
        computing_model_part = self._solver.GetComputingModelPart()

        with open('NodalModes.json') as f:
            data=json.load(f)
            counter = 0
            rom_dofs=self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()
            Dimensions = len(data['1.0'])  # dimensions are defined here (a vector or a matrix is imported to each node)
            for node in computing_model_part.Nodes:
                aux = KratosMultiphysics.Matrix(Dimensions, rom_dofs) 
                for j in range(Dimensions):
                    Counter=str(1.0*node.Id)
                    if Dimensions >1:
                        for i in range(rom_dofs):                         
                            aux[j,i] = data[Counter][j][i]
                    else:
                        for i in range(rom_dofs): 
                            aux[j,i] = data[Counter][i]
                node.SetValue(romapp.ROM_BASIS, aux ) # ROM basis
                node.SetValue(romapp.AUX_ID, counter) # Aux ID
                counter+=1

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
    Simulation = StructDynamicROM(model,parameters)
    Simulation.Run()