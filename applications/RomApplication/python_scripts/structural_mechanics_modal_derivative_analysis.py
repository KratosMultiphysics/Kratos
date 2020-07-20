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
        """ Create the Solver (and create and import the ModelPart if it is not already in the model) """
        ## Solver construction
        return solver_wrapper.CreateSolver(self.model, self.project_parameters)

    def _GetSimulationName(self):
        return "::[Modal Derivative Simulation]:: "
    
    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """

        # Creating output
        super(StructuralMechanicsModalDerivativeAnalysis, self).OutputSolutionStep()

        self.WriteNodalModes()
    
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

        derivative_parameter = self.project_parameters["solver_settings"]["derivative_parameter"].GetString()
        derivative_parameter_type_flag = True
        if derivative_parameter != "modal_coordinates":
            derivative_parameter_type_flag = False

        sub_model_parts_list = self.project_parameters["solver_settings"]["sub_model_parts_list"].GetStringArray()
        num_sub_model_parts = len(sub_model_parts_list)
        if not derivative_parameter_type_flag and num_sub_model_parts == 0:
            err_msg  = '"sub_model_parts_list" is empty for "derivative_parameter": "'+derivative_parameter+'".'
            err_msg  += ' Please provide at least one SubModelPart'
            raise Exception(err_msg)

        for sub_model_part_name in sub_model_parts_list:
            if not computing_model_part.HasSubModelPart(sub_model_part_name):
                err_msg  = '"'+sub_model_part_name+'" is not a SubModelPart of "'+computing_model_part.Name+'"'
                raise Exception(err_msg)
        
        rom_parameters_filename = self.project_parameters["solver_settings"]["rom_parameters_filename"].GetString()
        with open(rom_parameters_filename) as rom_parameters_file:
            rom_parameters = json.load(rom_parameters_file)

            number_of_initial_rom_dofs = rom_parameters["rom_settings"]["number_of_rom_dofs"]
            
            number_of_extended_rom_dofs = None
            if derivative_type_flag and derivative_parameter_type_flag:
                number_of_extended_rom_dofs = int(number_of_initial_rom_dofs * ( number_of_initial_rom_dofs + 1 ))
            elif not derivative_type_flag and derivative_parameter_type_flag: 
                number_of_extended_rom_dofs = int(number_of_initial_rom_dofs + number_of_initial_rom_dofs * ( number_of_initial_rom_dofs + 1 ) / 2)
            elif not derivative_parameter_type_flag:
                number_of_extended_rom_dofs = int((1+num_sub_model_parts)*number_of_initial_rom_dofs)

            eigenvalues = rom_parameters["eigenvalues"]
            kratos_eigenvalues = KratosMultiphysics.Vector(number_of_initial_rom_dofs)
            for i in range(number_of_initial_rom_dofs):
                kratos_eigenvalues[i] = eigenvalues[i]
            computing_model_part.ProcessInfo[StructuralMechanicsApplication.EIGENVALUE_VECTOR] = kratos_eigenvalues
                        
            counter = 0
            nodal_dofs = len(rom_parameters["rom_settings"]["nodal_unknowns"])
            nodal_modes = rom_parameters["nodal_modes"]
            for node in computing_model_part.Nodes:
                aux = KratosMultiphysics.Matrix(nodal_dofs, number_of_extended_rom_dofs)
                for i in range(nodal_dofs):
                    Counter=str(node.Id)
                    for j in range(number_of_initial_rom_dofs):
                        aux[i,j] = nodal_modes[Counter][i][j]
                node.SetValue(RomApplication.ROM_BASIS, aux ) # ROM basis
                node.SetValue(RomApplication.AUX_ID, counter) # Aux ID
                counter+=1

    def WriteNodalModes(self):

        # Iniate nodal modes dictionary
        rom_parameters = {}
        rom_parameters["rom_settings"] = {}
        rom_parameters["rom_settings"]["nodal_unknowns"] = ["ROTATION_X", "ROTATION_Y", "ROTATION_Z", "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z"]
        rom_parameters["nodal_modes"] = {}
        
        eigenvalues = self.model.GetModelPart('Structure').ProcessInfo[RomApplication.EIGENVALUE_VECTOR]
        rom_parameters["eigenvalues"] = []
        for eigenvalue in eigenvalues:
            rom_parameters["eigenvalues"].append(eigenvalue)

        # Loop over nodes
        for node in self.model.GetModelPart('Structure').GetNodes():

            # Create Node Id for NodalMode
            rom_parameters["nodal_modes"][str(node.Id)] = []

            # Get nodal eigenvector matrix
            nodal_modes_mtx = node.GetValue(RomApplication.ROM_BASIS)
            num_nodal_dofs = nodal_modes_mtx.Size1()
            num_nodal_modes = nodal_modes_mtx.Size2()
            
            # Loop over nodal DOFs
            for iDOF in range(0, num_nodal_dofs):
                
                # Initiate DOF modes
                dof_modes = []

                # Loop over dof modes                
                for iMode in range(0, num_nodal_modes):
                    dof_modes.append(nodal_modes_mtx[iDOF, iMode])

                # Set nodal modes
                rom_parameters["nodal_modes"][str(node.Id)].append(dof_modes)

        keys_nodal_modes = list(rom_parameters["nodal_modes"])
        num_modes = len(rom_parameters["nodal_modes"][keys_nodal_modes[0]][0])
        rom_parameters["rom_settings"]["number_of_rom_dofs"] = num_modes
        
        with open('RomParameters.json', 'w') as f:
            json.dump(rom_parameters,f, indent=2)

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Nodal modes printed in JSON format")
