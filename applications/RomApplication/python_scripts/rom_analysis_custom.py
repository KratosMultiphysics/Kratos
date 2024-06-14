import importlib

import KratosMultiphysics
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom
from KratosMultiphysics.RomApplication.hrom_training_utility import HRomTrainingUtility
from KratosMultiphysics.RomApplication.petrov_galerkin_training_utility import PetrovGalerkinTrainingUtility
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess
import numpy as np

from glob import glob
from os import remove
from pathlib import Path

def CreateCustomRomAnalysisInstance(cls, global_model, parameters):
    class CustomRomAnalysis(cls):

        def __init__(self,global_model, parameters):
            super().__init__(global_model, parameters)
            print('STARTING CUSTOM ANNPROM SIMULATION')

        def _EncodeSnapshot(self, s):
            """Define method to encode a given snapshot
            e.g. for POD: return self.U.T@u
            """

        def _DecodeSnapshot(self, q):
            """Define method to encode a given snapshot
            e.g. for POD: return self.U@q
            """

        def _GetSVDPhiMatrices(self):
            """Define method to obtain the Phi_inf and Phi_sup matrices"""

        def _GetNNLayers(self):
            """Define method to obtain the NN layers' weights"""

        def _GetReferenceSnapshot(self):
            """Define method to obtain the reference snapshot"""

        def _CreateSolver(self):
            """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """

            # Assign rom basis output folder and file name
            self.rom_basis_output_name = 'RomParameters' #Default
            self.rom_basis_output_folder = 'rom_data' #Default
            if self.project_parameters.Has("output_processes"):
                for name in self.project_parameters["output_processes"].keys():
                    if name=="rom_output":
                        rom_output_paramaters = self.project_parameters["output_processes"]["rom_output"]
                        self.rom_basis_output_name = rom_output_paramaters[0]["Parameters"]["rom_basis_output_name"].GetString()
                        self.rom_basis_output_folder = rom_output_paramaters[0]["Parameters"]["rom_basis_output_folder"].GetString()
            self.rom_basis_output_name = Path(self.rom_basis_output_name)
            self.rom_basis_output_folder = Path(self.rom_basis_output_folder)

            # Get the ROM settings from the RomParameters.json input file
            with open(self.rom_basis_output_folder / self.rom_basis_output_name.with_suffix('.json')) as rom_parameters:
                self.rom_parameters = KratosMultiphysics.Parameters(rom_parameters.read())

            # Set the ROM settings in the "solver_settings" of the solver introducing the physics
            self.project_parameters["solver_settings"].AddValue("rom_settings", self.rom_parameters["rom_settings"])

            # Inner ROM settings
            self.rom_bns_settings = self.rom_parameters["rom_settings"]["rom_bns_settings"] if self.rom_parameters["rom_settings"].Has("rom_bns_settings") else KratosMultiphysics.Parameters("""{}""")

            # Retrieve the ROM and HROM formats from settings, defaulting to 'json' for backward compatibility with earlier file formats and configurations.
            # TODO: Consider switching the default to 'numpy' for ROM and another more efficient format for HROM
            # for improved efficiency with large datasets in future versions.
            self.rom_format = self.rom_parameters["rom_format"].GetString() if self.rom_parameters.Has("rom_format") else "json"
            self.hrom_format = (self.rom_parameters["hrom_settings"]["hrom_format"].GetString()
                                if self.rom_parameters["hrom_settings"].Has("hrom_format")
                                else "json")

            # HROM operations flags
            self.rom_basis_process_list_check = True
            self.rom_basis_output_process_check = True
            self.run_hrom = self.rom_parameters["run_hrom"].GetBool() if self.rom_parameters.Has("run_hrom") else False
            self.train_hrom = self.rom_parameters["train_hrom"].GetBool() if self.rom_parameters.Has("train_hrom") else False
            self.rom_manager = self.rom_parameters["rom_manager"].GetBool() if self.rom_parameters.Has("rom_manager") else False
            if self.run_hrom and self.train_hrom:
                # Check that train an run HROM are not set at the same time
                err_msg = "\'run_hrom\' and \'train_hrom\' are both \'true\'. Select either training or running (if training has been already done)."
                raise Exception(err_msg)

            # Petrov Galerking Training settings
            self.train_petrov_galerkin = self.rom_bns_settings["train_petrov_galerkin"].GetBool() if self.rom_bns_settings.Has("train_petrov_galerkin") else False
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
                solving_technique = self.rom_bns_settings["solving_technique"].GetString() if self.rom_bns_settings.Has("solving_technique") else "normal_equations"
                # Check if the solving technique is either "normal_equations" or "qr_decomposition"
                if solving_technique not in ["normal_equations", "qr_decomposition"]:
                    err_msg = f"'{solving_technique}' is not a valid solving technique. Choose either 'normal_equations' or 'qr_decomposition'."
                    raise Exception(err_msg)

                self.project_parameters["solver_settings"]["rom_settings"]["rom_bns_settings"].AddString("solving_technique", solving_technique)
                self.project_parameters["solver_settings"]["rom_settings"]["rom_bns_settings"].AddBool("train_petrov_galerkin", self.train_petrov_galerkin)
                #Adding the basis strategy for generating the left ROB for the Petrov-Galerkin ROM.
                petrov_galerkin_basis_strategy = self.rom_bns_settings["basis_strategy"].GetString() if self.rom_bns_settings.Has("basis_strategy") else "residuals"
                self.project_parameters["solver_settings"]["rom_settings"]["rom_bns_settings"].AddString("basis_strategy", petrov_galerkin_basis_strategy)

            ##Petrov Galerkin
            if self.solving_strategy=="petrov_galerkin":
                self.petrov_galerkin_rom_dofs = self.project_parameters["solver_settings"]["rom_settings"]["petrov_galerkin_number_of_rom_dofs"].GetInt()
            else:
                self.project_parameters["solver_settings"]["rom_settings"].RemoveValue("petrov_galerkin_number_of_rom_dofs")

            # ROM assembling strategy
            self.assembling_strategy = self.rom_parameters["assembling_strategy"].GetString() if self.rom_parameters.Has("assembling_strategy") else "global"
            self.project_parameters["solver_settings"].AddString("assembling_strategy",self.assembling_strategy)

            # Create the ROM solver
            return python_solvers_wrapper_rom.CreateSolver(
                self.model,
                self.project_parameters)

        def _GetListOfProcesses(self):
            # Get the already existent processes list
            list_of_processes = super()._GetListOfProcesses()

            # Check if there is any instance of ROM basis output
            if self.rom_basis_process_list_check and not self.rom_manager: #eliminate the RomOutputProcess if not run inside the RomManager
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
            if self.rom_basis_output_process_check and not self.rom_manager: #eliminate the RomOutputProcess if not run inside the RomManager
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
            # nodal_modes = self.rom_parameters["nodal_modes"]
            nodal_dofs = len(self.project_parameters["solver_settings"]["rom_settings"]["nodal_unknowns"].GetStringArray())
            # rom_dofs = self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()

            # Set the right nodal ROM basis
            if self.rom_format == "json":
                print('JSON format not supported for ANNPROM architecture')
                exit()
            elif self.rom_format == "numpy":
                # Set the nodal ROM basis
                # node_ids = np.load(self.rom_basis_output_folder / "NodeIds.npy")
                node_ids = np.arange(computing_model_part.NumberOfNodes())+1
                # right_modes = np.load(self.rom_basis_output_folder / "RightBasisMatrix.npy")
                # if right_modes.ndim ==1: #check if matrix contains a single mode (a 1D numpy array)
                #     right_modes.reshape(-1,1)
                # right_modes = right_modes[:,:rom_dofs]

                # Handle left modes for Petrov-Galerkin strategy
                if self.solving_strategy == "petrov_galerkin":
                    petrov_galerkin_rom_dofs = self.project_parameters["solver_settings"]["rom_settings"]["petrov_galerkin_number_of_rom_dofs"].GetInt()
                    left_modes = np.load(self.rom_basis_output_folder / "LeftBasisMatrix.npy")

                    # Ensure left_modes is a 2D array
                    if left_modes.ndim == 1:
                        left_modes = left_modes.reshape(-1, 1)

                    # Trim left_modes based on Petrov-Galerkin ROM DOFs
                    left_modes = left_modes[:, :petrov_galerkin_rom_dofs]

                # Preprocess to create a mapping from node ID to index
                node_id_to_index = {nid: idx for idx, nid in enumerate(node_ids)}

                # Loop over each node in computing_model_part
                for node in computing_model_part.Nodes:
                    index = node_id_to_index[node.Id]
                    offset = index * nodal_dofs

                    # Set ROM basis for the node
                    # node.SetValue(KratosROM.ROM_BASIS, KratosMultiphysics.Matrix(right_modes[offset:offset+nodal_dofs, :])) # ROM basis

                    # Set Petrov-Galerkin ROM basis if applicable
                    if self.solving_strategy == "petrov_galerkin":
                        node.SetValue(KratosROM.ROM_LEFT_BASIS, KratosMultiphysics.Matrix(left_modes[offset:offset + nodal_dofs, :]))

            # Check for HROM stages
            if self.train_hrom:
                # Pass the name of the Rom Parameters file and folder
                self.rom_parameters.AddString("rom_basis_output_name", str(self.rom_basis_output_name))
                self.rom_parameters.AddString("rom_basis_output_folder", str(self.rom_basis_output_folder))
                # Create the training utility to calculate the HROM weights
                self.__hrom_training_utility = HRomTrainingUtility(
                    self._GetSolver(),
                    self.rom_parameters)
            elif self.run_hrom:
                if self.hrom_format == "json":
                    # Set the HROM weights in elements and conditions
                    hrom_weights_elements = self.rom_parameters["elements_and_weights"]["Elements"]
                    for key,value in zip(hrom_weights_elements.keys(), hrom_weights_elements.values()):
                        computing_model_part.GetElement(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, value.GetDouble()) #FIXME: FIX THE +1
                    hrom_weights_condtions = self.rom_parameters["elements_and_weights"]["Conditions"]
                    for key,value in zip(hrom_weights_condtions.keys(), hrom_weights_condtions.values()):
                        computing_model_part.GetCondition(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, value.GetDouble()) #FIXME: FIX THE +1
                elif self.hrom_format == "numpy":
                    # Set the HROM weights in elements and conditions
                    element_indexes = np.load(f"{self.rom_basis_output_folder}/HROM_ElementIds.npy")
                    element_weights = np.load(f"{self.rom_basis_output_folder}/HROM_ElementWeights.npy")
                    condition_indexes = np.load(f"{self.rom_basis_output_folder}/HROM_ConditionIds.npy")
                    conditon_weights = np.load(f"{self.rom_basis_output_folder}/HROM_ConditionWeights.npy")
                    for i in range(np.size(element_indexes)):
                        computing_model_part.GetElement(int( element_indexes[i])+1).SetValue(KratosROM.HROM_WEIGHT, element_weights[i]  ) #FIXME: FIX THE +1
                    for i in range(np.size(condition_indexes)):
                        computing_model_part.GetCondition(int( condition_indexes[i])+1).SetValue(KratosROM.HROM_WEIGHT, conditon_weights[i]  ) #FIXME: FIX THE +1


            # Check and Initialize Petrov Galerkin Training stage
            if self.train_petrov_galerkin:
                # Pass the name of the Rom Parameters file and folder
                self.rom_parameters.AddString("rom_basis_output_name", str(self.rom_basis_output_name))
                self.rom_parameters.AddString("rom_basis_output_folder", str(self.rom_basis_output_folder))
                self.__petrov_galerkin_training_utility = PetrovGalerkinTrainingUtility(
                    self._GetSolver(),
                    self.rom_parameters)
                

        def Initialize(self):
            super().Initialize()

        def ModifyInitialGeometry(self):
            super().ModifyInitialGeometry()

            NNLayers=self._GetNNLayers()
            SVDPhiMatrices=self._GetSVDPhiMatrices()
            refSnapshot=self._GetReferenceSnapshot()
            numberOfROMModes = SVDPhiMatrices[0].Size2()
            self._GetSolver()._GetBuilderAndSolver().SetNumberOfROMModes(numberOfROMModes)
            self._GetSolver()._GetBuilderAndSolver().SetNumberOfNNLayers(len(NNLayers))
            for i, layer in enumerate(NNLayers):
                self._GetSolver()._GetBuilderAndSolver().SetNNLayer(i, layer)
            self._GetSolver()._GetBuilderAndSolver().SetPhiMatrices(SVDPhiMatrices[0], SVDPhiMatrices[1], SVDPhiMatrices[2])
            self._GetSolver()._GetBuilderAndSolver().SetReferenceSnapshot(refSnapshot)
            
            computing_model_part = self._GetSolver().GetComputingModelPart().GetRootModelPart()
            nodal_unknown_names= self.project_parameters["solver_settings"]["rom_settings"]["nodal_unknowns"].GetStringArray()
            nodal_dofs = len(nodal_unknown_names)
            nodal_unknowns=[]

            for node_var_name in nodal_unknown_names:
                nodal_unknowns.append(KratosMultiphysics.KratosGlobals.GetVariable(node_var_name))

            s = []
            for node in computing_model_part.Nodes:
                for nodal_var in nodal_unknowns:
                    s.append(node.GetSolutionStepValue(nodal_var))

            print('Nodal variables: ', nodal_unknown_names)

            s_default = np.array(s, copy=False)
            print(s_default.shape)

            q, aux_norm_data = self._EncodeSnapshot(s_default)
            print(q.shape)

            computing_model_part.SetValue(KratosROM.ROM_SOLUTION_BASE, KratosMultiphysics.Vector(np.squeeze(q, axis=0)))
            computing_model_part.SetValue(KratosROM.ROM_SOLUTION_TOTAL, KratosMultiphysics.Vector(np.squeeze(q, axis=0)))
            computing_model_part.SetValue(KratosROM.ROM_SOLUTION_INCREMENT, KratosMultiphysics.Vector(np.zeros(q.shape[1])))

            s_init = np.squeeze(self._DecodeSnapshot(q), axis=0)
            print(s_init.shape)

            i = 0
            for node in computing_model_part.Nodes:
                for nodal_var in nodal_unknowns:
                    node.SetSolutionStepValue(nodal_var, s_init[i])
                    i+=1
                    
            computing_model_part.SetValue(KratosROM.SOLUTION_BASE, KratosMultiphysics.Vector(s_init))

            # Initialize nodal ROM_BASIS to zeros
            for node in computing_model_part.Nodes:
                node.SetValue(KratosROM.ROM_BASIS, np.zeros((nodal_dofs, numberOfROMModes)))


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


        def GetHROM_utility(self):
            return self.__hrom_training_utility

        def GetPetrovGalerkinTrainUtility(self):
            return self.__petrov_galerkin_training_utility


        def Finalize(self):
            # This calls the physics Finalize
            super().Finalize()

            # Once simulation is completed, calculate and save the HROM weights
            if self.train_hrom and not self.rom_manager:
                self.__hrom_training_utility.CalculateAndSaveHRomWeights()
                self.__hrom_training_utility.CreateHRomModelParts()

            # Once simulation is completed, calculate and save the Petrov Galerkin ROM basis
            if self.train_petrov_galerkin and not self.rom_manager:
                self.__petrov_galerkin_training_utility.CalculateAndSaveBasis()

    return CustomRomAnalysis(global_model, parameters)

if __name__ == "__main__":

    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    analysis_stage_module_name = parameters["analysis_stage"].GetString()
    analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
    analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

    analysis_stage_module = importlib.import_module(analysis_stage_module_name)
    analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

    global_model = KratosMultiphysics.Model()
    simulation = CreateCustomRomAnalysisInstance(analysis_stage_class, global_model, parameters)
    simulation.Run()