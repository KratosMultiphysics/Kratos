from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.RomApplication as RomApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom as solver_wrapper

import json

class StructuralMechanicsModalDerivativeAnalysis(StructuralMechanicsAnalysis):

    def __init__(self,model,project_parameters):
        super(StructuralMechanicsModalDerivativeAnalysis,self).__init__(model,project_parameters)

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        return solver_wrapper.CreateSolver(self.model, self.project_parameters)

    def _GetSimulationName(self):
        return "::[Modal Derivative Simulation]:: "
    
    def ModifyInitialGeometry(self):
        """Here is the place where the BASIS_ROM and the AUX_ID are imposed to each node"""
        super(StructuralMechanicsModalDerivativeAnalysis,self).ModifyInitialGeometry()
        computing_model_part = self._solver.GetComputingModelPart()

        derivative_type = self.project_parameters["solver_settings"]["derivative_type"].GetString()
        derivative_type_flag = None
        if derivative_type == "static":
            derivative_type_flag = False
        elif derivative_type == "dynamic":
            derivative_type_flag = True
        else:
            err_msg  = '\"derivative_type\" can only be \"static\" or \"dynamic\"'
            raise Exception(err_msg)
        
        derivative_parameter_type = self.project_parameters["solver_settings"]["derivative_parameter"].GetString()
        derivative_parameter_type_flag = True
        if derivative_parameter_type != "modal_coordinates":
            derivative_parameter_type_flag = False

        rom_parameters_filename = self.project_parameters["solver_settings"]["rom_parameters_filename"].GetString()
        with open(rom_parameters_filename) as rom_parameters_file:
            data = json.load(rom_parameters_file)

            number_of_initial_rom_dofs = data["rom_settings"]["number_of_rom_dofs"]
            
            if derivative_type_flag:
                eigenvalues = data["eigenvalues"]
            else:
                eigenvalues = [0]*number_of_initial_rom_dofs
                            
            number_of_extended_rom_dofs = None
            if derivative_type_flag and derivative_parameter_type_flag:
                number_of_extended_rom_dofs = int(number_of_initial_rom_dofs * ( number_of_initial_rom_dofs + 1 ))
            elif not derivative_type_flag and derivative_parameter_type_flag: 
                number_of_extended_rom_dofs = int(number_of_initial_rom_dofs + number_of_initial_rom_dofs * ( number_of_initial_rom_dofs + 1 ) / 2)
            elif not derivative_parameter_type_flag:
                number_of_extended_rom_dofs = int(2*number_of_initial_rom_dofs)

            kratos_eigenvalues = KratosMultiphysics.Vector(number_of_initial_rom_dofs)
            for i in range(number_of_initial_rom_dofs):
                kratos_eigenvalues[i] = eigenvalues[i]
            computing_model_part.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR] = kratos_eigenvalues
                        
            counter = 0
            nodal_dofs = len(data["rom_settings"]["nodal_unknowns"])
            nodal_modes = data["nodal_modes"]
            for node in computing_model_part.Nodes:
                aux = KratosMultiphysics.Matrix(nodal_dofs, number_of_extended_rom_dofs)
                for i in range(nodal_dofs):
                    Counter=str(node.Id)
                    for j in range(number_of_initial_rom_dofs):
                        aux[i,j] = nodal_modes[Counter][i][j]
                node.SetValue(RomApplication.ROM_BASIS, aux ) # ROM basis
                node.SetValue(RomApplication.AUX_ID, counter) # Aux ID
                counter+=1
