import importlib

import KratosMultiphysics
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics.RomApplication import new_python_solvers_wrapper_rom
from KratosMultiphysics.RomApplication.hrom_training_utility import HRomTrainingUtility
from KratosMultiphysics.RomApplication.petrov_galerkin_training_utility import PetrovGalerkinTrainingUtility
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess

from glob import glob 
from os import remove

def CreateRomAnalysisInstance(cls, global_model, parameters):
    class RomAnalysis(cls):

        def __init__(self,global_model, parameters):
            super().__init__(global_model, parameters)

        def _CreateSolver(self):
            """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """

            # Get the ROM settings from the RomParameters.json input file
            with open('RomParameters.json') as rom_parameters:
                self.rom_parameters = KratosMultiphysics.Parameters(rom_parameters.read())

            # Set the ROM settings in the "solver_settings" of the solver introducing the physics
            self.project_parameters["solver_settings"].AddValue("rom_settings", self.rom_parameters["rom_settings"])

            # HROM operations flags
            self.rom_basis_process_list_check = True
            self.rom_basis_output_process_check = True
            self.run_hrom = self.rom_parameters["run_hrom"].GetBool() if self.rom_parameters.Has("run_hrom") else False
            self.train_hrom = self.rom_parameters["train_hrom"].GetBool() if self.rom_parameters.Has("train_hrom") else False
            if self.run_hrom and self.train_hrom:
                # Check that train an run HROM are not set at the same time
                err_msg = "\'run_hrom\' and \'train_hrom\' are both \'true\'. Select either training or running (if training has been already done)."
                raise Exception(err_msg)

            # Petrov Galerking Training settings
            self.train_petrov_galerkin = self.rom_parameters["train_petrov_galerkin"]["train"].GetBool() if self.rom_parameters.Has("train_petrov_galerkin") else False
            if self.train_hrom and self.train_petrov_galerkin:
                err_msg = "\'train_petrov_galerkin\' and \'train_hrom\' are both \'true\'. Select only one training strategy."
                raise Exception(err_msg)
            if (self.train_hrom or self.train_petrov_galerkin) and (self.run_hrom):
                err_msg = "\'train_petrov_galerkin\' or \'train_hrom\' and \'run_hrom\' are both \'true\'. Select either training or running (if training has been already done)."
                raise Exception(err_msg)
            
            # ROM solving strategy
            self.solving_strategy = self.rom_parameters["projection_strategy"].GetString() if self.rom_parameters.Has("projection_strategy") else "galerkin"
            self.project_parameters["solver_settings"].AddString("projection_strategy",self.solving_strategy)

            # Add or remove parameters depending on the solving strategy
            ##LSPG
            if self.solving_strategy=="lspg":
                self.project_parameters["solver_settings"]["rom_settings"].AddBool("train_petrov_galerkin", self.train_petrov_galerkin)
            ##Petrov Galerkin
            if self.solving_strategy=="petrov_galerkin":
                self.petrov_galerkin_rom_dofs = self.project_parameters["solver_settings"]["rom_settings"]["petrov_galerkin_number_of_rom_dofs"].GetInt()
            else:    
                self.project_parameters["solver_settings"]["rom_settings"].RemoveValue("petrov_galerkin_number_of_rom_dofs")

            # Create the ROM solver
            return new_python_solvers_wrapper_rom.CreateSolver(
                self.model,
                self.project_parameters)

        def _GetListOfProcesses(self):
            # Get the already existent processes list
            list_of_processes = super()._GetListOfProcesses()

            # Check if there is any instance of ROM basis output
            if self.rom_basis_process_list_check:
                for process in list_of_processes:
                    if isinstance(process, KratosROM.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess):
                        warn_msg = "\'CalculateRomBasisOutputProcess\' instance found in ROM stage. Basis must be already stored in \'RomParameters.json\'. Removing instance from processes list."
                        KratosMultiphysics.Logger.PrintWarning("RomAnalysis", warn_msg)
                        list_of_processes.remove(process)
                self.rom_basis_process_list_check = False

            return list_of_processes

        def _GetListOfOutputProcesses(self):
            # Get the already existent output processes list
            list_of_output_processes = super()._GetListOfOutputProcesses()

            # Check if there is any instance of ROM basis output
            if self.rom_basis_output_process_check:
                for process in list_of_output_processes:
                    if isinstance(process, KratosROM.calculate_rom_basis_output_process.CalculateRomBasisOutputProcess):
                        warn_msg = "\'CalculateRomBasisOutputProcess\' instance found in ROM stage. Basis must be already stored in \'RomParameters.json\'. Removing instance from output processes list."
                        KratosMultiphysics.Logger.PrintWarning("RomAnalysis", warn_msg)
                        list_of_output_processes.remove(process)
                self.rom_basis_output_process_check = False

            return list_of_output_processes

        def _GetSimulationName(self):
            return "::[ROM Simulation]:: "

        def ModifyAfterSolverInitialize(self):
            """Here is where the ROM_BASIS is imposed to each node"""
            super().ModifyAfterSolverInitialize()

            # Get the model part where the ROM is to be applied
            computing_model_part = self._GetSolver().GetComputingModelPart().GetRootModelPart()
            # computing_model_part = self._GetSolver().GetComputingModelPart()

            # Set ROM basis
            nodal_modes = self.rom_parameters["nodal_modes"]
            nodal_dofs = len(self.project_parameters["solver_settings"]["rom_settings"]["nodal_unknowns"].GetStringArray())
            rom_dofs = self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()

            # Set the right nodal ROM basis 
            aux = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
            for node in computing_model_part.Nodes:
                node_id = str(node.Id)
                for j in range(nodal_dofs):
                    for i in range(rom_dofs):
                        aux[j,i] = nodal_modes[node_id][j][i].GetDouble()
                node.SetValue(KratosROM.ROM_BASIS, aux)
            
            # Set the left nodal ROM basis if it is different than the right nodal ROM basis (i.e. Petrov-Galerkin)
            if (self.solving_strategy == "petrov_galerkin"):
                petrov_galerkin_nodal_modes = self.rom_parameters["petrov_galerkin_nodal_modes"]
                petrov_galerkin_nodal_dofs = len(self.project_parameters["solver_settings"]["rom_settings"]["nodal_unknowns"].GetStringArray())
                aux = KratosMultiphysics.Matrix(petrov_galerkin_nodal_dofs, self.petrov_galerkin_rom_dofs)
                for node in computing_model_part.Nodes:
                    node_id = str(node.Id)
                    for j in range(petrov_galerkin_nodal_dofs):
                        for i in range(self.petrov_galerkin_rom_dofs):
                            aux[j,i] = petrov_galerkin_nodal_modes[node_id][j][i].GetDouble()
                    node.SetValue(KratosROM.ROM_LEFT_BASIS, aux)

            # Check for HROM stages
            if self.train_hrom:
                # Create the training utility to calculate the HROM weights
                self.__hrom_training_utility = HRomTrainingUtility(
                    self._GetSolver(),
                    self.rom_parameters)
            elif self.run_hrom:
                # Set the HROM weights in elements and conditions
                hrom_weights_elements = self.rom_parameters["elements_and_weights"]["Elements"]
                for key,value in zip(hrom_weights_elements.keys(), hrom_weights_elements.values()):
                    computing_model_part.GetElement(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, value.GetDouble()) #FIXME: FIX THE +1

                hrom_weights_condtions = self.rom_parameters["elements_and_weights"]["Conditions"]
                for key,value in zip(hrom_weights_condtions.keys(), hrom_weights_condtions.values()):
                    computing_model_part.GetCondition(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, value.GetDouble()) #FIXME: FIX THE +1

            # Check and Initialize Petrov Galerkin Training stage
            if self.train_petrov_galerkin:
                self.__petrov_galerkin_training_utility = PetrovGalerkinTrainingUtility(
                    self._GetSolver(),
                    self.rom_parameters)

        def FinalizeSolutionStep(self):
            if self.train_petrov_galerkin:
                self.__petrov_galerkin_training_utility.AppendCurrentStepProjectedSystem()
                ## Delete all .res.mm files when training Petrov-Galerkin with AssembledResiduals
                files_to_delete_list = glob('*.res.mm')
                for to_erase_file in files_to_delete_list:
                    remove(to_erase_file)

            # Call the HROM training utility to append the current step residuals
            # Note that this needs to be done prior to the other processes to avoid unfixing the BCs
            if self.train_hrom:
                self.__hrom_training_utility.AppendCurrentStepResiduals()

            # #FIXME: Make this optional. This must be a process
            # # Project the ROM solution onto the visualization modelparts
            # if self.run_hrom:
            #     model_part_name = self._GetSolver().settings["model_part_name"].GetString()
            #     visualization_model_part = self.model.GetModelPart("{}.{}Visualization".format(model_part_name, model_part_name))
            #     KratosROM.RomAuxiliaryUtilities.ProjectRomSolutionIncrementToNodes(
            #         self.rom_parameters["rom_settings"]["nodal_unknowns"].GetStringArray(),
            #         visualization_model_part)

            # This calls the physics FinalizeSolutionStep (e.g. BCs)
            super().FinalizeSolutionStep()

        def Finalize(self):
            # This calls the physics Finalize
            super().Finalize()

            # Once simulation is completed, calculate and save the Petrov Galerkin ROM basis
            if self.train_petrov_galerkin:
                self.__petrov_galerkin_training_utility.CalculateAndSaveBasis()

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
    simulation.Run()
