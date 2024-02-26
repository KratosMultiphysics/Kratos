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
        self.rom_training_parameters = self._SetRomTrainingParameters()
        self.hrom_training_parameters = self.SetHromTrainingParameters()
        self.CustomizeSimulation = CustomizeSimulation
        self.UpdateProjectParameters = UpdateProjectParameters
        self.UpdateMaterialParametersFile = UpdateMaterialParametersFile
        self.SetUpQuantityOfInterestContainers()


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
            self.QoI_Fit_FOM.append(simulation.GetFinalData())
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
            self.QoI_Fit_ROM.append(simulation.GetFinalData())
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
            self.QoI_Fit_HROM.append(simulation.GetFinalData())
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
            self.QoI_Test_FOM.append(simulation.GetFinalData())
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
            self.QoI_Test_ROM.append(simulation.GetFinalData())
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
            self.QoI_Test_HROM.append(simulation.GetFinalData())
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
            self.QoI_Run_FOM.append(simulation.GetFinalData())


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
            self.QoI_Run_ROM.append(simulation.GetFinalData())


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
            self.QoI_Run_HROM.append(simulation.GetFinalData())


    def _AddHromParametersToRomParameters(self,f):
        f["hrom_settings"]["element_selection_type"] = self.hrom_training_parameters["element_selection_type"].GetString()
        f["hrom_settings"]["element_selection_svd_truncation_tolerance"] = self.hrom_training_parameters["element_selection_svd_truncation_tolerance"].GetDouble()
        f["hrom_settings"]["create_hrom_visualization_model_part"] = self.hrom_training_parameters["create_hrom_visualization_model_part"].GetBool()
        f["hrom_settings"]["echo_level"] = self.hrom_training_parameters["echo_level"].GetInt()
        f["hrom_settings"]["include_condition_parents"] = self.hrom_training_parameters["include_condition_parents"].GetBool()
        f["hrom_settings"]["initial_candidate_elements_model_part_list"] = self.hrom_training_parameters["initial_candidate_elements_model_part_list"].GetStringArray()
        f["hrom_settings"]["initial_candidate_conditions_model_part_list"] = self.hrom_training_parameters["initial_candidate_conditions_model_part_list"].GetStringArray()
        f["hrom_settings"]["constraint_sum_weights"] = self.hrom_training_parameters["constraint_sum_weights"].GetBool()

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
                f["rom_settings"]['rom_bns_settings'] = self._SetGalerkinBnSParameters()
            elif simulation_to_run=='trainHROMGalerkin':
                f['train_hrom']=True
                f['run_hrom']=False
                f["rom_settings"]['rom_bns_settings'] = self._SetGalerkinBnSParameters()
            elif simulation_to_run=='runHROMGalerkin':
                f['projection_strategy']="galerkin"
                f['train_hrom']=False
                f['run_hrom']=True
                f["rom_settings"]['rom_bns_settings'] = self._SetGalerkinBnSParameters()
            elif simulation_to_run == 'lspg':
                f['train_hrom'] = False
                f['run_hrom'] = False
                f['projection_strategy'] = "lspg"
                f["rom_settings"]['rom_bns_settings'] = self._SetLSPGBnSParameters()
            elif simulation_to_run == 'trainHROMLSPG':
                f['train_hrom'] = True
                f['run_hrom'] = False
                f['projection_strategy'] = "lspg"
                f["rom_settings"]['rom_bns_settings'] = self._SetLSPGBnSParameters()
            elif simulation_to_run == 'runHROMLSPG':
                f['train_hrom'] = False
                f['run_hrom'] = True
                f['projection_strategy'] = "lspg"
                f["rom_settings"]['rom_bns_settings'] = self._SetLSPGBnSParameters()
            elif simulation_to_run == 'TrainPG':
                f['train_hrom'] = False
                f['run_hrom'] = False
                f['projection_strategy'] = "lspg"
                f["rom_settings"]['rom_bns_settings'] = self._SetLSPGBnSParameters()
                f["rom_settings"]['rom_bns_settings']['train_petrov_galerkin'] = True  # Override the default
            elif simulation_to_run=='PG':
                f['train_hrom']=False
                f['run_hrom']=False
                f['projection_strategy']="petrov_galerkin"
                f["rom_settings"]['rom_bns_settings'] = self._SetPetrovGalerkinBnSParameters()
            elif simulation_to_run=='trainHROMPetrovGalerkin':
                f['train_hrom']=True
                f['run_hrom']=False
                f['projection_strategy']="petrov_galerkin"
                f["rom_settings"]['rom_bns_settings'] = self._SetPetrovGalerkinBnSParameters()
            elif simulation_to_run=='runHROMPetrovGalerkin':
                f['train_hrom']=False
                f['run_hrom']=True
                f['projection_strategy']="petrov_galerkin"
                f["rom_settings"]['rom_bns_settings'] = self._SetPetrovGalerkinBnSParameters()
            else:
                raise Exception(f'Unknown flag "{simulation_to_run}" change for RomParameters.json')
            parameter_file.seek(0)
            json.dump(f,parameter_file,indent=4)
            parameter_file.truncate()

    def _SetGalerkinBnSParameters(self):
        # Retrieve the default parameters as a JSON string and parse it into a dictionary
        defaults_json = self._GetGalerkinBnSParameters()
        defaults = json.loads(defaults_json)

        # Ensure 'galerkin_rom_bns_settings' exists in ROM parameters
        if not self.general_rom_manager_parameters["ROM"].Has("galerkin_rom_bns_settings"):
            self.general_rom_manager_parameters["ROM"].AddEmptyValue("galerkin_rom_bns_settings")

        # Get the ROM parameters for Galerkin
        rom_params = self.general_rom_manager_parameters["ROM"]["galerkin_rom_bns_settings"]

        # Update defaults with any existing ROM parameters
        self._UpdateDefaultsWithRomParams(defaults, rom_params)

        return defaults

    def _SetLSPGBnSParameters(self):
        # Retrieve the default parameters as a JSON string and parse it into a dictionary
        defaults_json = self._GetLSPGBnSParameters()
        defaults = json.loads(defaults_json)

        # Ensure 'lspg_rom_bns_settings' exists in ROM parameters
        if not self.general_rom_manager_parameters["ROM"].Has("lspg_rom_bns_settings"):
            self.general_rom_manager_parameters["ROM"].AddEmptyValue("lspg_rom_bns_settings")

        # Get the ROM parameters for LSPG
        rom_params = self.general_rom_manager_parameters["ROM"]["lspg_rom_bns_settings"]

        # Update defaults with any existing ROM parameters
        self._UpdateDefaultsWithRomParams(defaults, rom_params)

        return defaults

    def _SetPetrovGalerkinBnSParameters(self):
        # Retrieve the default parameters as a JSON string and parse it into a dictionary
        defaults_json = self._GetPetrovGalerkinBnSParameters()
        defaults = json.loads(defaults_json)

        # Ensure 'petrov_galerkin_rom_bns_settings' exists in ROM parameters
        if not self.general_rom_manager_parameters["ROM"].Has("petrov_galerkin_rom_bns_settings"):
            self.general_rom_manager_parameters["ROM"].AddEmptyValue("petrov_galerkin_rom_bns_settings")

        # Get the ROM parameters for Petrov-Galerkin
        rom_params = self.general_rom_manager_parameters["ROM"]["petrov_galerkin_rom_bns_settings"]

        # Update defaults with any existing ROM parameters
        self._UpdateDefaultsWithRomParams(defaults, rom_params)

        return defaults

    def _UpdateDefaultsWithRomParams(self, defaults, rom_params):
        for key, default_value in defaults.items():
            if rom_params.Has(key):
                if isinstance(default_value, bool):
                    defaults[key] = rom_params[key].GetBool()
                elif isinstance(default_value, str):
                    defaults[key] = rom_params[key].GetString()
                elif isinstance(default_value, float):
                    defaults[key] = rom_params[key].GetDouble()
        return defaults

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



    def _SetRomTrainingParameters(self):
        defaults = self._GetDefaulRomBasisOutputParameters()
        defaults["Parameters"]["rom_manager"].SetBool(True)  # Set the flag to true when inside the RomManager to trigger particular behavior for multiple parameters

        rom_params = self.general_rom_manager_parameters["ROM"]

        keys_to_copy = [
            "svd_truncation_tolerance",
            "model_part_name",
            "rom_basis_output_format",
            "rom_basis_output_name",
            "rom_basis_output_folder",
            "nodal_unknowns",
            "snapshots_interval",
        ]

        for key in keys_to_copy:
            if key in rom_params.keys():
                defaults["Parameters"][key] = rom_params[key]

        return defaults




    def SetHromTrainingParameters(self):
        defaults = self._GetDefaulHromTrainingParameters()
        hrom_params = self.general_rom_manager_parameters["HROM"]

        keys_to_copy = [
            "element_selection_type",
            "element_selection_svd_truncation_tolerance",
            "create_hrom_visualization_model_part",
            "echo_level",
            "constraint_sum_weights",
        ]

        for key in keys_to_copy:
            if key in hrom_params.keys():
                defaults[key] = hrom_params[key]

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
                "include_condition_parents": false,
                "constraint_sum_weights": true
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


    def SetUpQuantityOfInterestContainers(self):
        #TODO implement more options if the QoI is too large to keep in RAM
        self.QoI_Fit_FOM = []
        self.QoI_Fit_ROM = []
        self.QoI_Fit_HROM = []
        self.QoI_Test_FOM = []
        self.QoI_Test_ROM = []
        self.QoI_Test_HROM = []
        self.QoI_Run_FOM = []
        self.QoI_Run_ROM = []
        self.QoI_Run_HROM = []


    def _GetGalerkinBnSParameters(self):
        # Define the default settings in JSON format for Galerkin BnS
        rom_bns_settings = """{
            "monotonicity_preserving": false
        }"""
        return rom_bns_settings

    def _GetPetrovGalerkinBnSParameters(self):
        # Define the default settings in JSON format for Petrov-Galerkin BnS
        rom_bns_settings = """{
            "monotonicity_preserving": false
        }"""
        return rom_bns_settings

    def _GetLSPGBnSParameters(self):
        # Define the default settings in JSON format
        # Comments:
        # - basis_strategy: Options include 'residuals', 'jacobian', 'reactions'
        # - solving_technique: Options include 'normal_equations', 'qr_decomposition'
        rom_bns_settings = """{
            "train_petrov_galerkin": false,
            "basis_strategy": "residuals",
            "include_phi": false,
            "svd_truncation_tolerance": 1e-8,
            "solving_technique": "normal_equations",
            "monotonicity_preserving": false
        }"""
        return rom_bns_settings