from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom as solver_wrapper
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import json

class StructuralMechanicsAnalysisROM(StructuralMechanicsAnalysis):

    def __init__(self,model,project_parameters):
        super(StructuralMechanicsAnalysisROM,self).__init__(model,project_parameters)
    
    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return solver_wrapper.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def _GetSimulationName(self):
        return "::[ROM Simulation]:: "
        
    def ModifyInitialGeometry(self):
        """Here is the place where the BASIS_ROM and the AUX_ID are imposed to each node"""
        super(StructuralMechanicsAnalysisROM,self).ModifyInitialGeometry()
        
        computing_model_part = self._solver.GetComputingModelPart()
        import pdb
        with open('NodalModes.json') as f:
            data=json.load(f)
            counter = 0
            rom_dofs=self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()
            Dimensions = len(data['1.0'])  #Checking for the number of DOFs (a vector or a matrix is imported to each node)
            for node in computing_model_part.Nodes:
                aux = KratosMultiphysics.Matrix(Dimensions, rom_dofs) 
                for j in range(Dimensions):
                    Counter=str(1.0*node.Id)
                    for i in range(rom_dofs):                         
                        aux[j,i] = data[Counter][j][i]
                node.SetValue(romapp.ROM_BASIS, aux ) # ROM basis
                node.SetValue(romapp.AUX_ID, counter) # Aux ID
                counter+=1
    


