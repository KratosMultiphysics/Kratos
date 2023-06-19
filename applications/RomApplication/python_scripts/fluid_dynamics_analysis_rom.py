import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.RomApplication as romapp
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom as solver_wrapper
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

import json
import numpy as np

class FluidDynamicsAnalysisROM(FluidDynamicsAnalysis):

    def __init__(self,model,project_parameters, path = '.', hyper_reduction_element_selector = None):
        KratosMultiphysics.Logger.PrintWarning('\x1b[1;31m[DEPRECATED CLASS] \x1b[0m',"\'FluidDynamicsAnalysisROM\'", "class is deprecated. Use the \'RomAnalysis\' one instead.")
        self.path = path
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
        with open(self.path +'/RomParameters.json') as rom_parameters:
            rom_settings = KratosMultiphysics.Parameters(rom_parameters.read())
            self.project_parameters["solver_settings"].AddValue("rom_settings", rom_settings["rom_settings"])
        return solver_wrapper.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def _GetSimulationName(self):
        return "::[ROM Simulation]:: "

    def ModifyAfterSolverInitialize(self):
        """Here is where the ROM_BASIS is imposed to each node"""
        super().ModifyAfterSolverInitialize()
        computing_model_part = self._solver.GetComputingModelPart()
        with open(self.path + '/RomParameters.json') as f:
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
                counter+=1
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                self.ResidualUtilityObject = romapp.RomResidualsUtility(self._GetSolver().GetComputingModelPart(), self.project_parameters["solver_settings"]["rom_settings"], self._GetSolver()._GetScheme())


    def FinalizeSolutionStep(self):
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                print('\n\n\n\nGenerating matrix of residuals')
                ResMat = self.ResidualUtilityObject.GetResiduals()
                NP_ResMat = np.array(ResMat, copy=False)
                self.time_step_residual_matrix_container.append(NP_ResMat)
        super().FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                self.residuals_snapshots = self._ObtainBasis()

    def _ObtainBasis(self):
        ### Building the Snapshot matrix ####
        for i in range (len(self.time_step_residual_matrix_container)):
            if i == 0:
                SnapshotMatrix = self.time_step_residual_matrix_container[i]
            else:
                SnapshotMatrix = np.c_[SnapshotMatrix, self.time_step_residual_matrix_container[i]]
        return SnapshotMatrix