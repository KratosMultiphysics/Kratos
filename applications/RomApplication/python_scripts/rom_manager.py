import KratosMultiphysics
from KratosMultiphysics.RomApplication.rom_testing_utilities import SetUpSimulationInstance
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
import numpy as np
import importlib
import json
from pathlib import Path



class RomManager(object):

    def __init__(self,project_parameters_name, general_rom_manager_parameters, CustomizeSimulation, UpdateProjectParameters,UpdateMaterialParametersFile):
        #FIXME:
        # - Use a method (upcoming) for smothly retrieving solutions. In here we are using the RomBasisOutput process in order to store the solutions
        # - There is some redundancy between the methods that launch the simulations. Can we create a single method?
        # - We must use the Hyper-Reduced model part for the HROM simulations. (The one in which we have eliminated the non selected elements and conditions)
        # - Not yet paralellised with COMPSs
        self.project_parameters_name = project_parameters_name
        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.rom_training_parameters = self._SetRomTrainintParameters()
        self.hrom_training_parameters = self.SetHromTrainingParameters()
        self.CustomizeSimulation = CustomizeSimulation
        self.UpdateProjectParameters = UpdateProjectParameters
        self.UpdateMaterialParametersFile = UpdateMaterialParametersFile



    def Fit(self, mu_train=[None], store_all_snapshots=False, store_fom_snapshots=False, store_rom_snapshots=False, store_hrom_snapshots=False, store_residuals_projected = False):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            if any(item == "ROM" for item in training_stages):
                fom_snapshots = self.__LaunchTrainROM(mu_train)
                if store_all_snapshots or store_fom_snapshots:
                    self._StoreSnapshotsMatrix('fom_snapshots', fom_snapshots)
                self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
                rom_snapshots = self.__LaunchROM(mu_train)
                if store_all_snapshots or store_rom_snapshots:
                    self._StoreSnapshotsMatrix('rom_snapshots', rom_snapshots)
                self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)

            if any(item == "HROM" for item in training_stages):
                #FIXME there will be an error if we only train HROM, but not ROM
                self._ChangeRomFlags(simulation_to_run = "trainHROMGalerkin")
                self.__LaunchTrainHROM(mu_train, store_residuals_projected)
                self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
                hrom_snapshots = self.__LaunchHROM(mu_train)
                if store_all_snapshots or store_hrom_snapshots:
                    self._StoreSnapshotsMatrix('hrom_snapshots', hrom_snapshots)
                self.ROMvsHROM_train = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)
        #######################

        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif chosen_projection_strategy == "lspg":
            if any(item == "ROM" for item in training_stages):
                fom_snapshots = self.__LaunchTrainROM(mu_train)
                if store_all_snapshots or store_fom_snapshots:
                    self._StoreSnapshotsMatrix('fom_snapshots', fom_snapshots)
                self._ChangeRomFlags(simulation_to_run = "lspg")
                rom_snapshots = self.__LaunchROM(mu_train)
                if store_all_snapshots or store_rom_snapshots:
                    self._StoreSnapshotsMatrix('rom_snapshots', rom_snapshots)
                self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)
            if any(item == "HROM" for item in training_stages):
                # Change the flags to train the HROM for LSPG
                self._ChangeRomFlags(simulation_to_run = "trainHROMLSPG")
                self.__LaunchTrainHROM(mu_train, store_residuals_projected)

                # Change the flags to run the HROM for LSPG
                self._ChangeRomFlags(simulation_to_run = "runHROMLSPG")
                hrom_snapshots = self.__LaunchHROM(mu_train)

                if store_all_snapshots or store_hrom_snapshots:
                    self._StoreSnapshotsMatrix('hrom_snapshots', hrom_snapshots)

                self.ROMvsHROM_train = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)
                #######################################

        ##########################
        ###  Petrov Galerkin   ###
        elif chosen_projection_strategy == "petrov_galerkin":
            if any(item == "ROM" for item in training_stages):
                fom_snapshots = self.__LaunchTrainROM(mu_train)
                if store_all_snapshots or store_fom_snapshots:
                    self._StoreSnapshotsMatrix('fom_snapshots', fom_snapshots)
                self._ChangeRomFlags(simulation_to_run = "TrainPG")
                self.__LaunchTrainPG(mu_train)
                self._ChangeRomFlags(simulation_to_run = "PG")
                rom_snapshots = self.__LaunchROM(mu_train)
                if store_all_snapshots or store_rom_snapshots:
                    self._StoreSnapshotsMatrix('rom_snapshots', rom_snapshots)
                self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)
            if any(item == "HROM" for item in training_stages):
                #FIXME there will be an error if we only train HROM, but not ROM
                self._ChangeRomFlags(simulation_to_run = "trainHROMPetrovGalerkin")
                self.__LaunchTrainHROM(mu_train, store_residuals_projected)
                self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
                hrom_snapshots = self.__LaunchHROM(mu_train)
                if store_all_snapshots or store_hrom_snapshots:
                    self._StoreSnapshotsMatrix('hrom_snapshots', hrom_snapshots)
                self.ROMvsHROM_train = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)
        ##########################
        else:
            err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
            raise Exception(err_msg)



    def Test(self, mu_test=[None]):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            if any(item == "ROM" for item in testing_stages):
                fom_snapshots = self.__LaunchTestFOM(mu_test)
                self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
                rom_snapshots = self.__LaunchTestROM(mu_test)
                self.ROMvsFOM_test = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)

            if any(item == "HROM" for item in testing_stages):
                #FIXME there will be an error if we only test HROM, but not ROM
                self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
                hrom_snapshots = self.__LaunchTestHROM(mu_test)
                self.ROMvsHROM_test = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)


        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif chosen_projection_strategy == "lspg":
            if any(item == "ROM" for item in testing_stages):
                fom_snapshots = self.__LaunchTestFOM(mu_test)
                self._ChangeRomFlags(simulation_to_run = "lspg")
                rom_snapshots = self.__LaunchTestROM(mu_test)
                self.ROMvsFOM_test = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)
            if any(item == "HROM" for item in testing_stages):
                self._ChangeRomFlags(simulation_to_run = "runHROMLSPG")
                hrom_snapshots = self.__LaunchTestHROM(mu_test)
                self.ROMvsHROM_test = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)
        #######################################


        ##########################
        ###  Petrov Galerkin   ###
        elif chosen_projection_strategy == "petrov_galerkin":
            if any(item == "ROM" for item in testing_stages):
                fom_snapshots = self.__LaunchTestFOM(mu_test)
                self._ChangeRomFlags(simulation_to_run = "PG")
                rom_snapshots = self.__LaunchTestROM(mu_test)
                self.ROMvsFOM_test = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)
            if any(item == "HROM" for item in testing_stages):
                #FIXME there will be an error if we only train HROM, but not ROM
                self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
                hrom_snapshots = self.__LaunchTestHROM(mu_test)
                self.ROMvsHROM_test = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)
        ##########################
        else:
            err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
            raise Exception(err_msg)




    def RunFOM(self, mu_run=[None]):
        self.__LaunchRunFOM(mu_run)

    def RunROM(self, mu_run=[None]):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif chosen_projection_strategy == "lspg":
            self._ChangeRomFlags(simulation_to_run = "lspg")
        ##########################
        ###  Petrov Galerkin   ###
        elif chosen_projection_strategy == "petrov_galerkin":
            self._ChangeRomFlags(simulation_to_run = "PG")
        else:
            err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
            raise Exception(err_msg)
        self.__LaunchRunROM(mu_run)

    def RunHROM(self, mu_run=[None], use_full_model_part = False):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif chosen_projection_strategy == "lspg":
            self._ChangeRomFlags(simulation_to_run = "runHROMLSPG")
        ##########################
        ###  Petrov Galerkin   ###
        elif chosen_projection_strategy == "petrov_galerkin":
            self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
        else:
            err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
            raise Exception(err_msg)
        self.__LaunchRunHROM(mu_run, use_full_model_part)


    def PrintErrors(self):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
        testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
        if any(item == "ROM" for item in training_stages):
            print("approximation error in train set FOM vs ROM: ", self.ROMvsFOM_train)
        if any(item == "HROM" for item in training_stages):
            print("approximation error in train set ROM vs HROM: ", self.ROMvsHROM_train)
        if any(item == "ROM" for item in testing_stages):
            print("approximation error in test set FOM vs ROM: ", self.ROMvsFOM_test)
        if any(item == "HROM" for item in testing_stages):
            print("approximation error in test set ROM vs HROM: ",  self.ROMvsHROM_test)


    def __LaunchTrainROM(self, mu_train):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_train):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy) #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters_copy = self._StoreResultsByName(parameters_copy,'FOM_Fit',mu,Id)
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters_copy)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)
        BasisOutputProcess._PrintRomBasis(SnapshotsMatrix) #Calling the RomOutput Process for creating the RomParameter.json

        return SnapshotsMatrix


    def __LaunchROM(self, mu_train):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_train):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)  #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters_copy = self._StoreResultsByName(parameters_copy,'ROM_Fit',mu,Id)
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix


    def __LaunchTrainPG(self, mu_train):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        PetrovGalerkinTrainMatrix = []
        for Id, mu in enumerate(mu_train):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)
            parameters_copy = self._StoreNoResults(parameters_copy)
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            PetrovGalerkinTrainMatrix.append(simulation.GetPetrovGalerkinTrainUtility()._GetSnapshotsMatrix()) #TODO is the best way of extracting the Projected Residuals calling the HROM residuals utility?
        simulation.GetPetrovGalerkinTrainUtility().CalculateAndSaveBasis(np.block(PetrovGalerkinTrainMatrix))


    def __LaunchTrainHROM(self, mu_train, store_residuals_projected=False):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        RedidualsSnapshotsMatrix = []
        for mu in mu_train:
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)
            parameters_copy = self._StoreNoResults(parameters_copy)
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            RedidualsSnapshotsMatrix.append(simulation.GetHROM_utility()._GetResidualsProjectedMatrix()) #TODO is the best way of extracting the Projected Residuals calling the HROM residuals utility?
        RedidualsSnapshotsMatrix = np.block(RedidualsSnapshotsMatrix)
        if store_residuals_projected:
            self._StoreSnapshotsMatrix('residuals_projected',RedidualsSnapshotsMatrix)
        u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(RedidualsSnapshotsMatrix,
        self.hrom_training_parameters["element_selection_svd_truncation_tolerance"].GetDouble())
        simulation.GetHROM_utility().hyper_reduction_element_selector.SetUp(u, InitialCandidatesSet = simulation.GetHROM_utility().candidate_ids)
        simulation.GetHROM_utility().hyper_reduction_element_selector.Run()
        if not simulation.GetHROM_utility().hyper_reduction_element_selector.success:
            KratosMultiphysics.Logger.PrintWarning("HRomTrainingUtility", "The Empirical Cubature Method did not converge using the initial set of candidates. Launching again without initial candidates.")
            #Imposing an initial candidate set can lead to no convergence. Restart without imposing the initial candidate set
            simulation.GetHROM_utility().hyper_reduction_element_selector.SetUp(u, InitialCandidatesSet = None)
            simulation.GetHROM_utility().hyper_reduction_element_selector.Run()
        simulation.GetHROM_utility().AppendHRomWeightsToRomParameters()
        simulation.GetHROM_utility().CreateHRomModelParts()


    def __LaunchHROM(self, mu_train):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_train):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)
            parameters_copy = self._StoreResultsByName(parameters_copy,'HROM_Fit',mu,Id)
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix


    def __LaunchTestFOM(self, mu_test):
        """
        This method should be parallel capable
        """
        # FIXME We must include a method to retrive solutions in the nodes and stop using the CalculateRomBasisOutputProcess to stote snapshots
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_test):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy) #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters_copy = self._StoreResultsByName(parameters_copy,'FOM_Test',mu,Id)
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters_copy)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)


        return SnapshotsMatrix


    def __LaunchTestROM(self, mu_test):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_test):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)  #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters_copy = self._StoreResultsByName(parameters_copy,'ROM_Test',mu,Id)
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix



    def __LaunchTestHROM(self, mu_test):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_test):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)
            parameters_copy = self._StoreResultsByName(parameters_copy,'HROM_Test',mu,Id)
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix


    def __LaunchRunFOM(self, mu_run):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        for Id, mu in enumerate(mu_run):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._StoreResultsByName(parameters_copy,'FOM_Run',mu,Id)
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters_copy)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()


    def __LaunchRunROM(self, mu_run):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        for Id, mu in enumerate(mu_run):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._StoreResultsByName(parameters_copy,'ROM_Run',mu,Id)
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()


    def __LaunchRunHROM(self, mu_run, use_full_model_part):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        if not use_full_model_part:
            model_part_name = parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
            parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(f"{model_part_name}HROM")

        for Id, mu in enumerate(mu_run):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._StoreResultsByName(parameters_copy,'HROM_Run',mu,Id)
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy)
            simulation.Run()


    def _AddHromParametersToRomParameters(self,f):
        f["rom_settings"]['rom_bns_settings']['monotonicity_preserving'] = self.general_rom_manager_parameters["ROM"]["galerkin_rom_bns_settings"]["monotonicity_preserving"].GetBool() if self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"].Has("monotonicity_preserving") else False
        f["hrom_settings"]["element_selection_type"] = self.general_rom_manager_parameters["HROM"]["element_selection_type"].GetString()
        f["hrom_settings"]["element_selection_svd_truncation_tolerance"] = self.general_rom_manager_parameters["HROM"]["element_selection_svd_truncation_tolerance"].GetDouble()
        f["hrom_settings"]["create_hrom_visualization_model_part"] = self.general_rom_manager_parameters["HROM"]["create_hrom_visualization_model_part"].GetBool()
        f["hrom_settings"]["echo_level"] = self.general_rom_manager_parameters["HROM"]["echo_level"].GetInt()
        f["hrom_settings"]["include_condition_parents"] = self.general_rom_manager_parameters["HROM"]["include_condition_parents"].GetBool() if self.general_rom_manager_parameters["HROM"].Has("include_condition_parents") else False
        f["hrom_settings"]["initial_candidate_elements_model_part_list"] = self.general_rom_manager_parameters["HROM"]["initial_candidate_elements_model_part_list"].GetStringArray() if self.general_rom_manager_parameters["HROM"].Has("initial_candidate_elements_model_part_list") else []
        f["hrom_settings"]["initial_candidate_conditions_model_part_list"] = self.general_rom_manager_parameters["HROM"]["initial_candidate_conditions_model_part_list"].GetStringArray() if self.general_rom_manager_parameters["HROM"].Has("initial_candidate_conditions_model_part_list") else []

    def _ChangeRomFlags(self, simulation_to_run = 'ROM'):
        """
        This method updates the Flags present in the RomParameters.json file
        for launching the correct part of the ROM workflow
        """
        #other options: "trainHROM", "runHROM"
        parameters_file_folder = self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].GetString() if self.general_rom_manager_parameters["ROM"].Has("rom_basis_output_folder") else "rom_data"
        parameters_file_name = self.general_rom_manager_parameters["ROM"]["rom_basis_output_name"].GetString() if self.general_rom_manager_parameters["ROM"].Has("rom_basis_output_name") else "RomParameters"

        # Convert to Path objects
        parameters_file_folder = Path(parameters_file_folder)
        parameters_file_name = Path(parameters_file_name)

        parameters_file_path = parameters_file_folder / parameters_file_name.with_suffix('.json')

        with parameters_file_path.open('r+') as parameter_file:
            f=json.load(parameter_file)
            f['assembling_strategy'] = self.general_rom_manager_parameters['assembling_strategy'].GetString() if self.general_rom_manager_parameters.Has('assembling_strategy') else 'global'
            self._AddHromParametersToRomParameters(f)
            if simulation_to_run=='GalerkinROM':
                f['projection_strategy']="galerkin"
                f['train_hrom']=False
                f['run_hrom']=False
                f["rom_settings"]['rom_bns_settings']['monotonicity_preserving'] = self.general_rom_manager_parameters["ROM"]["galerkin_rom_bns_settings"]["monotonicity_preserving"].GetBool() if self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"].Has("monotonicity_preserving") else False
            elif simulation_to_run=='trainHROMGalerkin':
                f['train_hrom']=True
                f['run_hrom']=False
                f["rom_settings"]['rom_bns_settings']['monotonicity_preserving'] = self.general_rom_manager_parameters["ROM"]["galerkin_rom_bns_settings"]["monotonicity_preserving"].GetBool() if self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"].Has("monotonicity_preserving") else False
            elif simulation_to_run=='runHROMGalerkin':
                f['projection_strategy']="galerkin"
                f['train_hrom']=False
                f['run_hrom']=True
                f["rom_settings"]['rom_bns_settings']['monotonicity_preserving'] = self.general_rom_manager_parameters["ROM"]["galerkin_rom_bns_settings"]["monotonicity_preserving"].GetBool() if self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"].Has("monotonicity_preserving") else False
            elif simulation_to_run == 'lspg':
                f['train_hrom'] = False
                f['run_hrom'] = False
                f['projection_strategy'] = "lspg"
                self.set_lspg_rom_bns_settings(f)
            elif simulation_to_run == 'trainHROMLSPG':
                f['train_hrom'] = True
                f['run_hrom'] = False
                f['projection_strategy'] = "lspg"
                self.set_lspg_rom_bns_settings(f)
            elif simulation_to_run == 'runHROMLSPG':
                f['train_hrom'] = False
                f['run_hrom'] = True
                f['projection_strategy'] = "lspg"
                self.set_lspg_rom_bns_settings(f)
            elif simulation_to_run == 'TrainPG':
                f['train_hrom'] = False
                f['run_hrom'] = False
                f['projection_strategy'] = "lspg"
                self.set_lspg_rom_bns_settings(f)
                f["rom_settings"]['rom_bns_settings']['train_petrov_galerkin'] = True  # Override the default
            elif simulation_to_run=='PG':
                f['train_hrom']=False
                f['run_hrom']=False
                f['projection_strategy']="petrov_galerkin"
                f["rom_settings"]['rom_bns_settings']['train_petrov_galerkin'] = False
            elif simulation_to_run=='trainHROMPetrovGalerkin':
                f['train_hrom']=True
                f['run_hrom']=False
                f['projection_strategy']="petrov_galerkin"
                f["rom_settings"]['rom_bns_settings']['train_petrov_galerkin'] = False
            elif simulation_to_run=='runHROMPetrovGalerkin':
                f['train_hrom']=False
                f['run_hrom']=True
                f['projection_strategy']="petrov_galerkin"
                f["rom_settings"]['rom_bns_settings']['train_petrov_galerkin'] = False
            else:
                raise Exception(f'Unknown flag "{simulation_to_run}" change for RomParameters.json')
            parameter_file.seek(0)
            json.dump(f,parameter_file,indent=4)
            parameter_file.truncate()


    def set_lspg_rom_bns_settings(self, f):
        settings = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]
        keys_with_defaults = {
            'train_petrov_galerkin': False,
            'basis_strategy': "residuals",
            'include_phi': False,
            'svd_truncation_tolerance': 1e-4,
            'solving_technique': "normal_equations",
            'monotonicity_preserving': False
        }

        for key, default_value in keys_with_defaults.items():
            if settings.Has(key):
                if isinstance(default_value, bool):
                    f["rom_settings"]['rom_bns_settings'][key] = settings[key].GetBool()
                elif isinstance(default_value, str):
                    f["rom_settings"]['rom_bns_settings'][key] = settings[key].GetString()
                elif isinstance(default_value, float):
                    f["rom_settings"]['rom_bns_settings'][key] = settings[key].GetDouble()
                # Add more types as needed
            else:
                f["rom_settings"]['rom_bns_settings'][key] = default_value



    def _AddBasisCreationToProjectParameters(self, parameters):
        #FIXME make sure no other rom_output already existed. If so, erase the prior and keep only the one in self.rom_training_parameters
        parameters["output_processes"].AddEmptyArray("rom_output")
        parameters["output_processes"]["rom_output"].Append(self.rom_training_parameters)

        return parameters




    def _StoreResultsByName(self,parameters,results_name,mu, Id):

        if  self.general_rom_manager_parameters["output_name"].GetString() == "mu":
            case_name = (", ".join(map(str, mu)))
        elif self.general_rom_manager_parameters["output_name"].GetString() == "id":
            case_name = str(Id)
        if self.general_rom_manager_parameters["save_gid_output"].GetBool():
            parameters["output_processes"]["gid_output"][0]["Parameters"]["output_name"].SetString('Results/'+ results_name +  case_name)
        else:
            parameters["output_processes"].RemoveValue("gid_output")
        if self.general_rom_manager_parameters["save_vtk_output"].GetBool():
            parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString('Results/vtk_output_'+ results_name + case_name)
        else:
            parameters["output_processes"].RemoveValue("vtk_output")

        return parameters


    def _StoreNoResults(self, parameters):
        parameters["output_processes"].RemoveValue("gid_output")
        parameters["output_processes"].RemoveValue("vtk_output")

        return parameters



    def _SetRomTrainintParameters(self):
        defaults = self._GetDefaulRomBasisOutputParameters()
        defaults["Parameters"]["rom_manager"].SetBool(True) #we set the flag to true when inside the RomManager to trigger particular behaviour for multiple parameters

        if self.general_rom_manager_parameters["ROM"].Has("svd_truncation_tolerance"):
            defaults["Parameters"]["svd_truncation_tolerance"].SetDouble(self.general_rom_manager_parameters["ROM"]["svd_truncation_tolerance"].GetDouble())
        if self.general_rom_manager_parameters["ROM"].Has("model_part_name"):
            defaults["Parameters"]["model_part_name"].SetString(self.general_rom_manager_parameters["ROM"]["model_part_name"].GetString())
        if self.general_rom_manager_parameters["ROM"].Has("rom_basis_output_format"):
            defaults["Parameters"]["rom_basis_output_format"].SetString(self.general_rom_manager_parameters["ROM"]["rom_basis_output_format"].GetString())
        if self.general_rom_manager_parameters["ROM"].Has("rom_basis_output_name"):
            defaults["Parameters"]["rom_basis_output_name"].SetString(self.general_rom_manager_parameters["ROM"]["rom_basis_output_name"].GetString())
        if self.general_rom_manager_parameters["ROM"].Has("rom_basis_output_folder"):
            defaults["Parameters"]["rom_basis_output_folder"].SetString(self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].GetString())
        if self.general_rom_manager_parameters["ROM"].Has("nodal_unknowns"):
            defaults["Parameters"]["nodal_unknowns"].SetStringArray(self.general_rom_manager_parameters["ROM"]["nodal_unknowns"].GetStringArray())
        if self.general_rom_manager_parameters["ROM"].Has("snapshots_interval"):
            defaults["Parameters"]["snapshots_interval"].SetDouble(self.general_rom_manager_parameters["ROM"]["snapshots_interval"].GetDouble())

        return defaults

    def SetHromTrainingParameters(self):
        defaults = self._GetDefaulHromTrainingParameters()
        defaults["element_selection_type"].SetString(self.general_rom_manager_parameters["HROM"]["element_selection_type"].GetString())
        defaults["element_selection_svd_truncation_tolerance"].SetDouble(self.general_rom_manager_parameters["HROM"]["element_selection_svd_truncation_tolerance"].GetDouble())
        defaults["create_hrom_visualization_model_part"].SetBool(self.general_rom_manager_parameters["HROM"]["create_hrom_visualization_model_part"].GetBool())
        defaults["echo_level"].SetInt(self.general_rom_manager_parameters["HROM"]["echo_level"].GetInt())

        return defaults


    def _GetAnalysisStageClass(self, parameters):

        analysis_stage_module_name = parameters["analysis_stage"].GetString()
        analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
        analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

        analysis_stage_module = importlib.import_module(analysis_stage_module_name)
        analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

        return analysis_stage_class


    def _GetDefaulRomBasisOutputParameters(self):
        rom_training_parameters = KratosMultiphysics.Parameters("""{
                "python_module" : "calculate_rom_basis_output_process",
                "kratos_module" : "KratosMultiphysics.RomApplication",
                "process_name"  : "CalculateRomBasisOutputProcess",
                "help"          : "This process should write the Rom basis",
                "Parameters"    :
                {
                    "model_part_name": "",
                    "rom_manager" : false,      // set to false for manual manipulation of ROM via flags in the RomParameters
                    "snapshots_control_type": "step",
                    "snapshots_interval": 1.0,
                    "nodal_unknowns":  [],
                    "rom_basis_output_format": "json",
                    "rom_basis_output_name": "RomParameters",
                    "rom_basis_output_folder": "rom_data",
                    "svd_truncation_tolerance": 1e-3
                }
            }""")
        return rom_training_parameters



    def _GetDefaulHromTrainingParameters(self):
        hrom_training_parameters = KratosMultiphysics.Parameters("""{
                "hrom_format": "numpy",
                "element_selection_type": "empirical_cubature",
                "element_selection_svd_truncation_tolerance": 1.0e-6,
                "echo_level" : 0,
                "create_hrom_visualization_model_part" : true,
                "projection_strategy": "galerkin",
                "include_conditions_model_parts_list": [],
                "include_elements_model_parts_list": [],
                "initial_candidate_elements_model_part_list" : [],
                "initial_candidate_conditions_model_part_list" : [],
                "include_nodal_neighbouring_elements_model_parts_list":[],
                "include_minimum_condition": false,
                "include_condition_parents": false
            }""")
        return hrom_training_parameters



    def _StoreSnapshotsMatrix(self, string_numpy_array_name, numpy_array):

        # Define the directory and file path
        rom_output_folder_name = self.rom_training_parameters["Parameters"]["rom_basis_output_folder"].GetString()
        directory = Path(rom_output_folder_name) / 'SnapshotsMatrices'
        file_path = directory / f'{string_numpy_array_name}.npy'

        # Create the directory if it doesn't exist
        directory.mkdir(parents=True, exist_ok=True)

        #save the array inside the chosen directory
        np.save(file_path, numpy_array)
