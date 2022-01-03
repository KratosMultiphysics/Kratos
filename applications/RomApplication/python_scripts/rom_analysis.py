import json
import importlib
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom
from KratosMultiphysics.RomApplication import new_python_solvers_wrapper_rom
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod

def CreateRomAnalysisInstance(cls, global_model, parameters, hyper_reduction_element_selector = None):
    class RomAnalysis(cls):

        def __init__(self,global_model, parameters, hyper_reduction_element_selector = None):
            super().__init__(global_model, parameters)

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

        def _CreateSolver(self):
            """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """

            # Get the ROM settings from the RomParameters.json input file and set them in the "solver_settings" of the solver introducing the physics
            with open('RomParameters.json') as rom_parameters:
                rom_settings = KratosMultiphysics.Parameters(rom_parameters.read())
                self.project_parameters["solver_settings"].AddValue("rom_settings", rom_settings["rom_settings"])

            # Create the ROM solver
            return new_python_solvers_wrapper_rom.CreateSolver(
                self.model,
                self.project_parameters)

        def _GetSimulationName(self):
            return "::[ROM Simulation]:: "

        def ModifyAfterSolverInitialize(self):
            """Here is where the ROM_BASIS is imposed to each node"""
            super().ModifyAfterSolverInitialize()

            # Get the model part where the ROM is to be applied
            computing_model_part = self._GetSolver().GetComputingModelPart()

            # Set ROM basis
            with open('RomParameters.json') as f:
                # Get the ROM data from RomParameters.json
                data = json.load(f)
                nodal_modes = data["nodal_modes"]
                nodal_dofs = len(self.project_parameters["solver_settings"]["rom_settings"]["nodal_unknowns"].GetStringArray())
                rom_dofs = self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()

                # Set the nodal ROM basis
                aux = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
                for node in computing_model_part.Nodes:
                    node_id = str(node.Id)
                    for j in range(nodal_dofs):
                        for i in range(rom_dofs):
                            aux[j,i] = nodal_modes[node_id][j][i]
                    node.SetValue(KratosROM.ROM_BASIS, aux)

            # Hyper-reduction
            if self.hyper_reduction_element_selector:
                if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                    self.ResidualUtilityObject = KratosROM.RomResidualsUtility(
                        computing_model_part,
                        self.project_parameters["solver_settings"]["rom_settings"],
                        self._GetSolver()._GetScheme())

        def FinalizeSolutionStep(self):
            if self.hyper_reduction_element_selector:
                if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                    KratosMultiphysics.Logger.PrintInfo("RomAnalysis","Generating matrix of residuals.")
                    res_mat = self.ResidualUtilityObject.GetResiduals()
                    np_res_mat = np.array(res_mat, copy=False)
                    self.time_step_residual_matrix_container.append(np_res_mat)

            super().FinalizeSolutionStep()

        def Finalize(self):
            super().Finalize()

            if self.hyper_reduction_element_selector:
                if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                    original_number_of_elements = self._GetSolver().GetComputingModelPart().NumberOfElements()
                    input_filename = self._GetSolver().settings["model_import_settings"]["input_filename"].GetString()
                    self. hyper_reduction_element_selector.SetUp(
                        self.time_step_residual_matrix_container,
                        original_number_of_elements,
                        input_filename)
                    self.hyper_reduction_element_selector.Run()

    return RomAnalysis(global_model, parameters, hyper_reduction_element_selector)

def CreateHRomAnalysisInstance(cls, global_model, parameters, hyper_reduction_element_selector = None):
    # Create an standard ROM analysis instance to get the type
    # This is required as RomAnalysis is defined inside CreateRomAnalysisInstance function
    rom_analysis = CreateRomAnalysisInstance(cls, global_model, parameters, hyper_reduction_element_selector)

    # Extend the ModifyAfterSolverInitialize method to set the element and condition hyper-reduction weights
    class HRomAnalysis(type(rom_analysis)):

        def ModifyAfterSolverInitialize(self):
            super().ModifyAfterSolverInitialize()
            computing_model_part = self._GetSolver().GetComputingModelPart()

            with open('ElementsAndWeights.json') as f:
                HR_data = json.load(f)
                for key in HR_data["Elements"].keys():
                    computing_model_part.GetElement(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, HR_data["Elements"][key])
                for key in HR_data["Conditions"].keys():
                    computing_model_part.GetCondition(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, HR_data["Conditions"][key])

    return HRomAnalysis(global_model, parameters, hyper_reduction_element_selector)

if __name__ == "__main__":

    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

    global_model = KratosMultiphysics.Model()
    simulation = CreateRomAnalysisInstance(analysis_stage_class, global_model, parameters)
    # simulation = CreateHRomAnalysisInstance(analysis_stage_class, global_model, parameters)
    simulation.Run()
