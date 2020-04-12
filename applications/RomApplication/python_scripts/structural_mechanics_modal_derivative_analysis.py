from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.RomApplication as RomApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

from KratosMultiphysics.StructuralMechanicsApplication import python_solvers_wrapper_structural as structural_solvers

import json

class StructuralMechanicsModalDerivativeAnalysis(StructuralMechanicsAnalysis):

def __init__(self,model,project_parameters):
        super(StructuralMechanicsAnalysisROM,self).__init__(model,project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return structural_solvers.CreateSolver(self.model, self.project_parameters)

    def _GetSimulationName(self):
        return "::[Modal Derivative Simulation]:: "
    
    def ModifyInitialGeometry(self):
        """Here is the place where the BASIS_ROM and the AUX_ID are imposed to each node"""
        super(ModalDerivativeSolver,self).ModifyInitialGeometry()
        computing_model_part = self._solver.GetComputingModelPart()
        with open('RomParameters.json') as f:
            data = json.load(f)
            eigenvalues = data["eigenvalues"]
            number_of_eigenvalues = len(eigenvalues)
            kratos_eigenvalues = KratosMultiphysics.Vector(number_of_eigenvalues)
            computing_model_part.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR] = kratos_eigenvalues
            nodal_dofs = len(data["rom_settings"]["nodal_unknowns"])
            nodal_modes = data["nodal_modes"]
            counter = 0
            initial_rom_dofs = self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()
            modified_rom_dofs = initial_rom_dofs * ( initial_rom_dofs + 1 )
            for node in computing_model_part.Nodes:
                aux = KratosMultiphysics.Matrix(nodal_dofs, modified_rom_dofs)
                for i in range(nodal_dofs):
                    Counter=str(node.Id)
                    for j in range(initial_rom_dofs):
                        aux[i,j] = nodal_modes[Counter][i][j]
                node.SetValue(romapp.ROM_BASIS, aux ) # ROM basis
                node.SetValue(romapp.AUX_ID, counter) # Aux ID
                counter+=1
