import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom as solver_wrapper
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import json
import numpy as np

class StructuralMechanicsAnalysisROM(StructuralMechanicsAnalysis):

    def __init__(self,model,project_parameters, hyper_reduction_element_selector = None):
        super().__init__(model,project_parameters)
        if hyper_reduction_element_selector != None :
            if hyper_reduction_element_selector == "EmpiricalCubature":
                self.hyper_reduction_element_selector = EmpiricalCubatureMethod()
                self.time_step_residual_matrix_container = []
            else:
                err_msg =  "The requested element selection method \"" + hyper_reduction_element_selector + "\" is not in the rom application\n"
                err_msg += "Available options are: \"EmpiricalCubature\""
                raise Exception(err_msg)
        else:
            self.hyper_reduction_element_selector = None

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        with open('RomParameters.json') as rom_parameters:
            rom_settings = KratosMultiphysics.Parameters(rom_parameters.read())
            self.project_parameters["solver_settings"].AddValue("rom_settings", rom_settings["rom_settings"])
        return solver_wrapper.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def _GetSimulationName(self):
        return "::[ROM Simulation]:: "

    def ModifyInitialGeometry(self):
        """Here is the place where the BASIS_ROM and the AUX_ID are imposed to each node"""
        super().ModifyInitialGeometry()
        computing_model_part = self._solver.GetComputingModelPart()
        with open('RomParameters.json') as f:
            data = json.load(f)
            nodal_dofs = len(data["rom_settings"]["nodal_unknowns"])
            nodal_modes = data["nodal_modes"]
            counter = 0
            rom_dofs= self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()
            for node in computing_model_part.Nodes:
                aux = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
                for j in range(nodal_dofs):
                    Counter=str(node.Id)
                    for i in range(rom_dofs):
                        aux[j,i] = nodal_modes[Counter][j][i]
                node.SetValue(romapp.ROM_BASIS, aux ) # ROM basis
                node.SetValue(romapp.AUX_ID, counter) # Aux ID
                counter+=1


    def ModifyAfterSolverInitialize(self):
        super().ModifyAfterSolverInitialize()
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                self.ResidualUtilityObject = romapp.RomResidualsUtility(self._GetSolver().GetComputingModelPart(), self.project_parameters["solver_settings"]["rom_settings"], KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                print('\n\n\n\nGenerating matrix of residuals')
                ResMat = self.ResidualUtilityObject.GetResiduals()
                NP_ResMat = np.array(ResMat, copy=False)
                self.time_step_residual_matrix_container.append(NP_ResMat)

    def Finalize(self):
        super().FinalizeSolutionStep()
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                OriginalNumberOfElements = self._GetSolver().GetComputingModelPart().NumberOfElements()
                ModelPartName = self._GetSolver().settings["model_import_settings"]["input_filename"].GetString()
                self. hyper_reduction_element_selector.SetUp(self.time_step_residual_matrix_container, OriginalNumberOfElements, ModelPartName)
                self.hyper_reduction_element_selector.Run()





