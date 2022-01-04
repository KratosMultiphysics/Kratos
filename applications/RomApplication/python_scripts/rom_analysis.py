import json
import importlib
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics.RomApplication import new_python_solvers_wrapper_rom
from KratosMultiphysics.RomApplication.hrom_training_utility import HRomTrainingUtility

def CreateRomAnalysisInstance(cls, global_model, parameters):
    class RomAnalysis(cls):

        def __init__(self,global_model, parameters):
            super().__init__(global_model, parameters)

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

                # Check for HROM stages
                run_hrom = data["run_hrom"].GetString() if data.Has("run_hrom") else False
                self.train_hrom = data["train_hrom"].GetString() if data.Has("train_rom") else False
                if run_hrom and self.train_hrom:
                    # Check that train an run HROM are not set at the same time
                    err_msg = "\'run_hrom\' and \'train_hrom\' are both \'true\'. Select either training or running (if training has been already done)."
                    raise Exception(err_msg)
                elif self.train_hrom:
                    # Create the training utility to calculate the HROM weights
                    self.__hrom_training_utility = HRomTrainingUtility(
                        self._GetSolver(),
                        data)
                elif run_hrom:
                    # Set the HROM weights in elements and conditions
                    for key in data["elements_and_weights"]["Elements"].keys():
                        computing_model_part.GetElement(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, data["elements_and_weights"]["Elements"][key])
                    for key in data["elements_and_weights"]["Conditions"].keys():
                        computing_model_part.GetCondition(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, data["elements_and_weights"]["Conditions"][key])

        def FinalizeSolutionStep(self):
            # Call the HROM training utility to append the current step residuals
            # Note that this needs to be done prior to the other processes to avoid unfixing the BCs
            if self.train_hrom:
                self.__hrom_training_utility.AppendCurrentStepResiduals()

            # This calls the physics FinalizeSolutionStep (e.g. BCs)
            super().FinalizeSolutionStep()

        def Finalize(self):
            # This calls the physics Finalize
            super().Finalize()

            # Once simulation is completed, calculate and save the HROM weights
            if self.train_hrom:
                self.__hrom_training_utility.CalculateAndSaveHRomWeights()
                self.__hrom_training_utility.CreateHRomModelParts()

    return RomAnalysis(global_model, parameters)

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
