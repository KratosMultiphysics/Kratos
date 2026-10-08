import numpy as np
import importlib
import shutil
import json
import shutil
from pathlib import Path
import types
import KratosMultiphysics
from KratosMultiphysics.RomApplication.rom_database import RomDatabase
from KratosMultiphysics.RomApplication.rom_testing_utilities import SetUpSimulationInstance
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
from KratosMultiphysics.RomApplication.rom_nn_interface import NN_ROM_Interface
import re


class RomManager(object):

    def __init__(self,project_parameters_name="ProjectParameters.json", general_rom_manager_parameters=None, CustomizeSimulation=None, UpdateProjectParameters=None,UpdateMaterialParametersFile=None, mu_names=None):
        #FIXME:
        # - There is some redundancy between the methods that launch the simulations. Can we create a single method?
        # - Not yet paralellised with COMPSs
        self.project_parameters_name = project_parameters_name
        self._SetUpRomManagerParameters(general_rom_manager_parameters)
        if CustomizeSimulation is None:
            self.CustomizeSimulation = self.DefaultCustomizeSimulation
        else:
            self.CustomizeSimulation = CustomizeSimulation
        if UpdateProjectParameters is None:
            self.UpdateProjectParameters = self.DefaultUpdateProjectParameters
        else:
            self.UpdateProjectParameters = UpdateProjectParameters
        if UpdateMaterialParametersFile is None:
            self.UpdateMaterialParametersFile = self.DefaultUpdateMaterialParametersFile
        else:
            self.UpdateMaterialParametersFile = UpdateMaterialParametersFile
        self.SetUpQuantityOfInterestContainers()
        self.data_base = RomDatabase(self.general_rom_manager_parameters, mu_names)
        self.SetupErrorsDictionaries()

    def Fit(self, mu_train=[None],mu_validation=[None]):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
        type_of_decoder = self.general_rom_manager_parameters["type_of_decoder"].GetString()
        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            if type_of_decoder =="ann_enhanced":
                if any(item == "ROM" for item in training_stages):
                    self._LaunchTrainROM(mu_train)
                    self._LaunchFOM(mu_validation) #What to do here with the gid and vtk results?
                    self.TrainAnnEnhancedROM(mu_train,mu_validation)
                    self._ChangeRomFlags(simulation_to_run = "GalerkinROM_ANN")
                    nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base)
                    self._LaunchROM(mu_train, nn_rom_interface=nn_rom_interface)
                if any(item == "HROM" for item in training_stages):
                    self._ChangeRomFlags(simulation_to_run = "trainHROMGalerkin_ANN")
                    self._LaunchTrainHROM(mu_train, nn_rom_interface=nn_rom_interface)
                    self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin_ANN")
                    self._LaunchHROM(mu_train,nn_rom_interface=nn_rom_interface)
            elif type_of_decoder =="linear":
                if any(item == "ROM" for item in training_stages):
                    self._LaunchTrainROM(mu_train)
                    self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
                    self._LaunchROM(mu_train)
                if any(item == "HROM" for item in training_stages):
                    #FIXME there will be an error if we only train HROM, but not ROM
                    self._ChangeRomFlags(simulation_to_run = "trainHROMGalerkin")
                    self._LaunchTrainHROM(mu_train)
                    self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
                    self._LaunchHROM(mu_train)
        #######################

        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif chosen_projection_strategy == "lspg":
            if type_of_decoder =="ann_enhanced":
                if any(item == "ROM" for item in training_stages):
                    self._LaunchTrainROM(mu_train)
                    self._LaunchFOM(mu_validation) #What to do here with the gid and vtk results?
                    self.TrainAnnEnhancedROM(mu_train,mu_validation)
                    self._ChangeRomFlags(simulation_to_run = "lspg_ANN")
                    nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base)
                    self._LaunchROM(mu_train, nn_rom_interface=nn_rom_interface)
                if any(item == "HROM" for item in training_stages):
                    self._ChangeRomFlags(simulation_to_run = "trainHROMlspg_ANN")
                    self._LaunchTrainHROM(mu_train, nn_rom_interface=nn_rom_interface)
                    self._ChangeRomFlags(simulation_to_run = "runHROMlspg_ANN")
                    self._LaunchHROM(mu_train,nn_rom_interface=nn_rom_interface)
            elif type_of_decoder =="linear":
                if any(item == "ROM" for item in training_stages):
                    self._LaunchTrainROM(mu_train)
                    self._ChangeRomFlags(simulation_to_run = "lspg")
                    self._LaunchROM(mu_train)
                if any(item == "HROM" for item in training_stages):
                    # Change the flags to train the HROM for LSPG
                    self._ChangeRomFlags(simulation_to_run = "trainHROMLSPG")
                    self._LaunchTrainHROM(mu_train)
                    # Change the flags to run the HROM for LSPG
                    self._ChangeRomFlags(simulation_to_run = "runHROMLSPG")
                    self._LaunchHROM(mu_train)
        #######################################

        ##########################
        ###  Petrov Galerkin   ###
        elif chosen_projection_strategy == "petrov_galerkin":
            if type_of_decoder =="ann_enhanced":
                err_msg = f'ann_enhanced rom only available for Galerkin Rom and LSPG ROM for the moment'
                raise Exception(err_msg)
            elif type_of_decoder =="linear":
            ##########################
                if any(item == "ROM" for item in training_stages):
                    self._LaunchTrainROM(mu_train)
                    self._ChangeRomFlags(simulation_to_run = "TrainPG")
                    self._LaunchTrainPG(mu_train)
                    self._ChangeRomFlags(simulation_to_run = "PG")
                    self._LaunchROM(mu_train)

                if any(item == "HROM" for item in training_stages):
                    #FIXME there will be an error if we only train HROM, but not ROM
                    self._ChangeRomFlags(simulation_to_run = "trainHROMPetrovGalerkin")
                    self._LaunchTrainHROM(mu_train)
                    self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
                    self._LaunchHROM(mu_train)
            ##########################
        else:
            err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
            raise Exception(err_msg)
        self.ComputeErrors(mu_train)

    def TrainAnnEnhancedROM(self, mu_train, mu_validation):
        counter = 0
        self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["online"]["model_number"].SetInt(counter)
        in_database, _ = self.data_base.check_if_in_database("Neural_Network", mu_train)
        if not in_database:
            self._LaunchTrainNeuralNetwork(mu_train,mu_validation)
        elif in_database and self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["training"]["retrain_if_exists"].GetBool():
            while in_database:
                counter+=1
                self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["online"]["model_number"].SetInt(counter)
                #using Fit(), the model launched will be the one trained last.
                #For Test() or Run() methods, it is the one privided in "model_number"
                in_database, _ = self.data_base.check_if_in_database("Neural_Network", mu_train)
            self._LaunchTrainNeuralNetwork(mu_train,mu_validation)

    def TestNeuralNetworkReconstruction(self, mu_train, mu_validation):
        self._LaunchTestNeuralNetworkReconstruction( mu_train, mu_validation)


    def Test(self, mu_test=[None], mu_train=[None]):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
        type_of_decoder = self.general_rom_manager_parameters["type_of_decoder"].GetString()

        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            if type_of_decoder =="ann_enhanced":
                nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base)
                if any(item == "ROM" for item in testing_stages):
                    self._LoadSolutionBasis(mu_train)
                    self._LaunchFOM(mu_test, gid_and_vtk_name='FOM_Test')
                    self._ChangeRomFlags(simulation_to_run = "GalerkinROM_ANN")
                    self._LaunchROM(mu_test, gid_and_vtk_name='ROM_Test', nn_rom_interface=nn_rom_interface)
                if any(item == "HROM" for item in testing_stages):
                    self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin_ANN")
                    self._LaunchHROM(mu_test, nn_rom_interface=nn_rom_interface, gid_and_vtk_name='HROM_Test')
            elif type_of_decoder =="linear":
                if any(item == "ROM" for item in testing_stages):
                    self._LoadSolutionBasis(mu_train)
                    self._LaunchFOM(mu_test, gid_and_vtk_name='FOM_Test')
                    self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
                    self._LaunchROM(mu_test, gid_and_vtk_name='ROM_Test')
                if any(item == "HROM" for item in testing_stages):
                    #FIXME there will be an error if we only test HROM, but not ROM
                    self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
                    self._LaunchHROM(mu_test,gid_and_vtk_name='HROM_Test')

        #######################

        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif chosen_projection_strategy == "lspg":
            if type_of_decoder =="ann_enhanced":
                nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base)
                if any(item == "ROM" for item in testing_stages):
                    self._LoadSolutionBasis(mu_train)
                    self._LaunchFOM(mu_test, gid_and_vtk_name='FOM_Test')
                    self._ChangeRomFlags(simulation_to_run = "lspg_ANN")
                    self._LaunchROM(mu_test, gid_and_vtk_name='ROM_Test', nn_rom_interface=nn_rom_interface)
                if any(item == "HROM" for item in testing_stages):
                    self._ChangeRomFlags(simulation_to_run = "runHROMlspg_ANN")
                    self._LaunchHROM(mu_test, nn_rom_interface=nn_rom_interface, gid_and_vtk_name='HROM_Test')
            elif type_of_decoder =="linear":
                if any(item == "ROM" for item in testing_stages):
                    self._LoadSolutionBasis(mu_train)
                    self._LaunchFOM(mu_test,gid_and_vtk_name='FOM_Test')
                    self._ChangeRomFlags(simulation_to_run = "lspg")
                    self._LaunchROM(mu_test,gid_and_vtk_name='ROM_Test')
                if any(item == "HROM" for item in testing_stages):
                    self._ChangeRomFlags(simulation_to_run = "runHROMLSPG")
                    self._LaunchHROM(mu_test,gid_and_vtk_name='HROM_Test')
        #######################################


        ##########################
        ###  Petrov Galerkin   ###
        elif chosen_projection_strategy == "petrov_galerkin":
            if type_of_decoder =="ann_enhanced":
                err_msg = f'ann_enhanced rom only available for Galerkin Rom and LSPG ROM for the moment'
                raise Exception(err_msg)
            elif type_of_decoder =="linear":
                if any(item == "ROM" for item in testing_stages):
                    self._LoadSolutionBasis(mu_train)
                    self._LaunchFOM(mu_test,gid_and_vtk_name='FOM_Test')
                    self._ChangeRomFlags(simulation_to_run = "PG")
                    self._LaunchROM(mu_test,gid_and_vtk_name='ROM_Test')
                if any(item == "HROM" for item in testing_stages):
                    #FIXME there will be an error if we only train HROM, but not ROM
                    self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
                    self._LaunchHROM(mu_test,gid_and_vtk_name='HROM_Test')
        ##########################
        else:
            err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
            raise Exception(err_msg)
        self.ComputeErrors(mu_test, 'Test')


    def RunFOM(self, mu_run=[None]):
        self._LaunchRunFOM(mu_run)


    def RunROM(self, mu_run=[None], mu_train=[None]):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        type_of_decoder = self.general_rom_manager_parameters["type_of_decoder"].GetString()
        nn_rom_interface = None
        self._LoadSolutionBasis(mu_train)
        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            if type_of_decoder =="ann_enhanced":
                self._ChangeRomFlags(simulation_to_run = "GalerkinROM_ANN")
                nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base)
            elif type_of_decoder =="linear":
                self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif chosen_projection_strategy == "lspg":
            if type_of_decoder =="ann_enhanced":
                self._ChangeRomFlags(simulation_to_run = "lspg_ANN")
                nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base)
            elif type_of_decoder =="linear":
                self._ChangeRomFlags(simulation_to_run = "lspg")
        ##########################
        ###  Petrov Galerkin   ###
        elif chosen_projection_strategy == "petrov_galerkin":
            if type_of_decoder =="ann_enhanced":
                err_msg = f'ann_enhanced rom only available for Galerkin Rom and LSPG ROM for the moment'
                raise Exception(err_msg)
            elif type_of_decoder =="linear":
                self._ChangeRomFlags(simulation_to_run = "PG")
        #########################################
        else:
            err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
            raise Exception(err_msg)
        self._LaunchRunROM(mu_run, nn_rom_interface=nn_rom_interface)




    def RunHROM(self, mu_run=[None], mu_train=[None], use_full_model_part = False):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        type_of_decoder = self.general_rom_manager_parameters["type_of_decoder"].GetString()
        nn_rom_interface = None
        self._LoadSolutionBasis(mu_train)
        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            if type_of_decoder =="ann_enhanced":
                self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin_ANN")
                nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base)
            elif type_of_decoder =="linear":
                self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif chosen_projection_strategy == "lspg":
            if type_of_decoder =="ann_enhanced":
                self._ChangeRomFlags(simulation_to_run = "runHROMlspg_ANN")
                nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base)
            elif type_of_decoder =="linear":
                self._ChangeRomFlags(simulation_to_run = "runHROMLSPG")
        ##########################
        ###  Petrov Galerkin   ###
        elif chosen_projection_strategy == "petrov_galerkin":
            if type_of_decoder =="ann_enhanced":
                err_msg = f'HROM only supports linear projection strategy for the moment'
                raise Exception(err_msg)
            elif type_of_decoder =="linear":
                self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
        else:
            err_msg = f'Provided projection strategy {chosen_projection_strategy} is not supported. Available options are \'galerkin\', \'lspg\' and \'petrov_galerkin\'.'
            raise Exception(err_msg)
        self._LaunchRunHROM(mu_run, use_full_model_part,nn_rom_interface)


    def ExportAnnEnhancedRom(self, mu_train=[None], export_folder="rom_data_standalone"):
        """
        Exports the ANN-enhanced ROM trained with 'mu_train' to a folder that can be run without the RomManager and its database.
        The exported folder replaces the ROM folder (default 'rom_data') of the case to be run standalone.
        The HROM files are copied if they exist. Set 'run_hrom' to true in the exported RomParameters to use them.
        """
        if self.general_rom_manager_parameters["type_of_decoder"].GetString() != "ann_enhanced":
            raise Exception("'ExportAnnEnhancedRom' requires 'type_of_decoder' to be 'ann_enhanced'.")
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        if chosen_projection_strategy == "galerkin":
            simulation_to_run = "GalerkinROM_ANN"
        elif chosen_projection_strategy == "lspg":
            simulation_to_run = "lspg_ANN"
        else:
            raise Exception(f'ann_enhanced rom only available for Galerkin Rom and LSPG ROM, not for \'{chosen_projection_strategy}\'.')

        # Write the RomParameters of the ANN-enhanced ROM in the ROM folder
        self._LoadSolutionBasis(mu_train)
        self._ChangeRomFlags(simulation_to_run = simulation_to_run)
        nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base)

        rom_folder = Path(self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].GetString())
        rom_parameters_file_name = Path(self.general_rom_manager_parameters["ROM"]["rom_basis_output_name"].GetString()).with_suffix('.json')
        export_folder = Path(export_folder)
        export_folder.mkdir(parents=True, exist_ok=True)

        # Files read by NN_ROM_Interface.FromNumpyFiles
        np.save(export_folder / "RightBasisMatrix.npy", nn_rom_interface.phi[:, :nn_rom_interface.n_sup])
        np.save(export_folder / "SingularValues.npy", nn_rom_interface.sigma[:nn_rom_interface.n_sup])
        shutil.copy(nn_rom_interface.network_weights_path, export_folder / "model_weights.npy")
        for file_name in ["NodeIds.npy", "HROM_ElementIds.npy", "HROM_ElementWeights.npy", "HROM_ConditionIds.npy", "HROM_ConditionWeights.npy"]:
            if (rom_folder / file_name).exists():
                shutil.copy(rom_folder / file_name, export_folder / file_name)

        with open(rom_folder / rom_parameters_file_name, 'r') as parameter_file:
            rom_parameters = json.load(parameter_file)
        rom_parameters["rom_manager"] = False
        rom_parameters["ann_enhanced_settings"] = {"modes" : [nn_rom_interface.n_inf, nn_rom_interface.n_sup]}
        with open(export_folder / rom_parameters_file_name, 'w') as parameter_file:
            json.dump(rom_parameters, parameter_file, indent=4)



    def ComputeErrors(self, mu_list, case="Fit"):
        fom_snapshots = self.data_base.get_snapshots_matrix_from_database(mu_list, table_name=f'FOM')
        if case=="Fit":
            stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
        elif case=="Test":
            stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
        stages = {"ROM", "HROM"} & set(stages)  # Ensures only "ROM" or "HROM" if present
        rom_snapshots = None
        hrom_snapshots = None
        if "ROM" in stages:
            rom_snapshots = self.data_base.get_snapshots_matrix_from_database(mu_list, table_name=f'ROM')
            error_rom_fom = np.linalg.norm(fom_snapshots - rom_snapshots) / np.linalg.norm(fom_snapshots)
            self.ROMvsFOM[case] = error_rom_fom
        if "HROM" in stages:
            if rom_snapshots is None:  # Only fetch if not already fetched
                rom_snapshots = self.data_base.get_snapshots_matrix_from_database(mu_list, table_name=f'ROM')
            hrom_snapshots = self.data_base.get_snapshots_matrix_from_database(mu_list, table_name=f'HROM')
            error_rom_hrom = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)
            error_fom_hrom = np.linalg.norm(fom_snapshots - hrom_snapshots) / np.linalg.norm(fom_snapshots)
            self.ROMvsHROM[case] = error_rom_hrom
            self.FOMvsHROM[case] = error_fom_hrom



    def PrintErrors(self):
        training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
        testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()

        training_set = set(training_stages)
        testing_set = set(testing_stages)

        # Check in Fit
        if "ROM" in training_set:
            self.aux_print_errors(self.ROMvsFOM['Fit'], 'train', 'FOM vs ROM')
        if "HROM" in training_set:
            self.aux_print_errors(self.ROMvsHROM['Fit'], 'train', 'ROM vs HROM')
            self.aux_print_errors(self.FOMvsHROM['Fit'], 'train', 'FOM vs HROM')

        # Check in Test
        if "ROM" in testing_set:
            self.aux_print_errors(self.ROMvsFOM['Test'], 'test', 'FOM vs ROM')
        if "HROM" in testing_set:
            self.aux_print_errors(self.ROMvsHROM['Test'], 'test', 'ROM vs HROM')
            self.aux_print_errors(self.FOMvsHROM['Test'], 'test', 'FOM vs HROM')

    def aux_print_errors(self, error, train_or_test, comparison_in_string):
        message = f"approximation error in {train_or_test} set {comparison_in_string}"
        if error is None:
            print(f"{message} not computed")
        else:
            print(f"{message}: {error}")

    def SetupErrorsDictionaries(self):
        self.ROMvsFOM = {'Fit': None, 'Test': None}
        self.ROMvsHROM = {'Fit': None, 'Test': None}
        self.FOMvsHROM = {'Fit': None, 'Test': None}


    def _LaunchTrainROM(self, mu_train):
        """
        This method should be parallel capable
        """
        self._LaunchFOM(mu_train)
        self._LaunchComputeSolutionBasis(mu_train)



    def _LaunchFOM(self, mu_train, gid_and_vtk_name='FOM_Fit'):
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        NonConvergedSolutionsGathering = self.general_rom_manager_parameters["store_nonconverged_fom_solutions"].GetBool()
        for Id, mu in enumerate(mu_train):
            fom_in_database, _ = self.data_base.check_if_in_database("FOM", mu)
            nonconverged_fom_in_database, _ = self.data_base.check_if_in_database("NonconvergedFOM", mu)
            if not fom_in_database or (NonConvergedSolutionsGathering and not nonconverged_fom_in_database):
                parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
                parameters_copy = self._AddNumpyOutputToProjectParameters(parameters_copy)
                parameters_copy = self._StoreResultsByName(parameters_copy,gid_and_vtk_name,mu,Id)
                materials_file_name = self._GetMaterialsFileName(parameters_copy)
                self.UpdateMaterialParametersFile(materials_file_name, mu)
                model = KratosMultiphysics.Model()
                analysis_stage_class = self._GetAnalysisStageClass(parameters_copy)
                simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)
                if NonConvergedSolutionsGathering:
                    simulation = self.ActivateNonconvergedSolutionsGathering(simulation)
                simulation.Run()
                if not fom_in_database:
                    self.data_base.add_to_database("QoI_FOM", mu, simulation.GetFinalData())
                    SnapshotsMatrix = self._GetSnapshotsMatrixFromNumpyOutput()
                    self.data_base.add_to_database("FOM", mu, SnapshotsMatrix)
                if NonConvergedSolutionsGathering:
                    self.data_base.add_to_database("NonconvergedFOM", mu, simulation.GetNonconvergedSolutions())



    def _LaunchComputeSolutionBasis(self, mu_train):
        in_database, hash_basis = self.data_base.check_if_in_database("RightBasis", mu_train)
        if not in_database:
            if self.general_rom_manager_parameters["ROM"]["use_non_converged_sols"].GetBool():
                u,sigma = self._ComputeSVD(self.data_base.get_snapshots_matrix_from_database(mu_train, table_name='NonconvergedFOM')) #TODO this might be too large for single opeartion, add partitioned svd
            else:
                u,sigma = self._ComputeSVD(self.data_base.get_snapshots_matrix_from_database(mu_train, table_name='FOM'))
            self._PrintRomBasis(u, sigma)
            self.data_base.add_to_database("RightBasis", mu_train, u )
            self.data_base.add_to_database("SingularValues_Solution", mu_train, sigma )
        else:
            _ , hash_sigma = self.data_base.check_if_in_database("SingularValues_Solution", mu_train)
            self._PrintRomBasis(self.data_base.get_single_numpy_from_database(hash_basis), self.data_base.get_single_numpy_from_database(hash_sigma) ) #this updates the RomParameters.json
        self.GenerateDatabaseSummary()

    def _LoadSolutionBasis(self, mu_train):
        in_database, hash_basis = self.data_base.check_if_in_database("RightBasis", mu_train)
        basis_directory = self.data_base.database_root_directory.parent
        basis_path = basis_directory / 'RightBasisMatrix.npy'
        basis_exists_in_parent_directory = basis_path.exists()

        if not in_database and basis_exists_in_parent_directory:
            err_msg = f'ROM basis not in RomDatabase, using the one in the working directory {basis_path}'
            KratosMultiphysics.Logger.PrintWarning(err_msg)
        elif not in_database and not basis_exists_in_parent_directory:
            err_msg = f'ROM basis not found for indicated training snapshots. Please run the Train method() to create it.'
            raise Exception(err_msg)
        else:
            _ , hash_sigma = self.data_base.check_if_in_database("SingularValues_Solution", mu_train)
            self._PrintRomBasis(self.data_base.get_single_numpy_from_database(hash_basis), self.data_base.get_single_numpy_from_database(hash_sigma) ) #this updates the RomParameters.json
        self.GenerateDatabaseSummary()



    def _LaunchROM(self, mu_train, gid_and_vtk_name='ROM_Fit', nn_rom_interface = None):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        for Id, mu in enumerate(mu_train):
            in_database, _ = self.data_base.check_if_in_database("ROM", mu)
            if not in_database:
                parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
                parameters_copy = self._AddNumpyOutputToProjectParameters(parameters_copy)
                parameters_copy = self._StoreResultsByName(parameters_copy,gid_and_vtk_name,mu,Id)
                materials_file_name = self._GetMaterialsFileName(parameters_copy)
                self.UpdateMaterialParametersFile(materials_file_name, mu)
                model = KratosMultiphysics.Model()
                analysis_stage_class = type(self._SetUpRomSimulationInstance(model, parameters_copy, nn_rom_interface=nn_rom_interface))
                simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)

                simulation.Run()
                self.data_base.add_to_database("QoI_ROM", mu, simulation.GetFinalData())
                SnapshotsMatrix = self._GetSnapshotsMatrixFromNumpyOutput()
                self.data_base.add_to_database("ROM", mu, SnapshotsMatrix )

        self.GenerateDatabaseSummary()



    def _LaunchTrainPG(self, mu_train):
        """
        This method should be parallel capable
        """
        in_database, hash_basis =  self.data_base.check_if_in_database("LeftBasis", mu_train)
        if not in_database:
            with open(self.project_parameters_name,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            PetrovGalerkinTrainingUtility = None
            for Id, mu in enumerate(mu_train):
                in_database, _ = self.data_base.check_if_in_database("PetrovGalerkinSnapshots", mu)
                if not in_database:
                    parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
                    parameters_copy = self._StoreNoResults(parameters_copy)
                    materials_file_name = self._GetMaterialsFileName(parameters_copy)
                    self.UpdateMaterialParametersFile(materials_file_name, mu)
                    model = KratosMultiphysics.Model()
                    analysis_stage_class = type(self._SetUpRomSimulationInstance(model, parameters_copy))
                    simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)
                    simulation.Run()
                    PetrovGalerkinTrainingUtility = simulation.GetPetrovGalerkinTrainUtility()
                    pretrov_galerkin_matrix = PetrovGalerkinTrainingUtility._GetSnapshotsMatrix() #TODO is the best way of extracting the Projected Residuals calling the HROM residuals utility?
                    self.data_base.add_to_database("PetrovGalerkinSnapshots", mu, pretrov_galerkin_matrix)
            if PetrovGalerkinTrainingUtility is None:
                PetrovGalerkinTrainingUtility = self.InitializeDummySimulationForPetrovGalerkinTrainingUtility()
            snapshots_matrix = self.data_base.get_snapshots_matrix_from_database(mu_train, table_name="PetrovGalerkinSnapshots")
            u = PetrovGalerkinTrainingUtility._CalculateResidualBasis(snapshots_matrix)
            PetrovGalerkinTrainingUtility._AppendNewBasisToRomParameters(u)
            self.data_base.add_to_database("LeftBasis", mu_train, u )
        else:
            PetrovGalerkinTrainingUtility = self.InitializeDummySimulationForPetrovGalerkinTrainingUtility()
            PetrovGalerkinTrainingUtility._AppendNewBasisToRomParameters(self.data_base.get_single_numpy_from_database(hash_basis)) #this updates the RomParameters.json

        self.GenerateDatabaseSummary()



    def _LaunchTrainHROM(self, mu_train, nn_rom_interface = None):
        """
        This method should be parallel capable
        The projected residuals are written by a ProjectedResidualsOutputProcess (one per sub-solver in coupled problems),
        and a set of HROM weights is computed for each of them (i.e. one column of HROM weights per sub-solver).
        """
        in_database_elems, hash_z =  self.data_base.check_if_in_database("HROM_Elements", mu_train)
        in_database_weights, hash_w =  self.data_base.check_if_in_database("HROM_Weights", mu_train)
        if not in_database_elems and not in_database_weights:
            with open(self.project_parameters_name,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            residuals_output_settings = self._GetResidualsProjectedOutputSettings()
            # The residuals are collected by the output processes, so the HRomTrainingUtility does not collect them
            self._SetTrainHROMFlag(False)
            for mu in mu_train:
                in_database, _ = self.data_base.check_if_in_database("ResidualsProjected", mu)
                if not in_database:
                    parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
                    parameters_copy = self._AddResidualsProjectedOutputProcessToProjectParameters(parameters_copy)  #This deals with the creation of residuals
                    parameters_copy = self._StoreNoResults(parameters_copy)
                    # Remove the residuals of previous runs, so that only the ones of this mu are collected
                    for _, _, folder in residuals_output_settings:
                        for residual_file in folder.glob("Residual_*.npy"):
                            residual_file.unlink()
                    materials_file_name = self._GetMaterialsFileName(parameters_copy)
                    self.UpdateMaterialParametersFile(materials_file_name, mu)
                    model = KratosMultiphysics.Model()
                    analysis_stage_class = type(self._SetUpRomSimulationInstance(model, parameters_copy, nn_rom_interface=nn_rom_interface))
                    simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)
                    simulation.Run()
                    ResidualProjected = self._GetResidualProjected() # this method fetches the residuals projected from the corresponding folder.
                    self.data_base.add_to_database("ResidualsProjected", mu, ResidualProjected )
                    # TODO for later PR: Implement path-based database storage to avoid RAM bottlenecks and erase the temporary .npy files from disk here.
            self._SetTrainHROMFlag(True)
            RedidualsSnapshotsMatrix = self.data_base.get_snapshots_matrix_from_database(mu_train, table_name="ResidualsProjected")
            HROM_utility = self.InitializeDummySimulationForHromTrainingUtility(nn_rom_interface=nn_rom_interface)
            z, w = self._LocalEmpiricalCubature(np.split(RedidualsSnapshotsMatrix, len(residuals_output_settings)), HROM_utility.candidate_ids)
            HROM_utility.hyper_reduction_element_selector.z = z
            HROM_utility.hyper_reduction_element_selector.w = w
            self.data_base.add_to_database("HROM_Elements", mu_train, np.squeeze(z))
            self.data_base.add_to_database("HROM_Weights", mu_train, np.squeeze(w))
        # elif (in_database_elems and not in_database_weights) or (not in_database_elems and in_database_weights):
        #     error #if one of the two is there but not the other, the database is corrupted
        else:
            HROM_utility = self.InitializeDummySimulationForHromTrainingUtility(nn_rom_interface=nn_rom_interface)
            #doing this ensures the elements and weights npy files contained in the rom_data folder are the ones to use.i.e. they are not from old runs of the rom manager
            HROM_utility.hyper_reduction_element_selector.w = self.data_base.get_single_numpy_from_database(hash_w)
            HROM_utility.hyper_reduction_element_selector.z = self.data_base.get_single_numpy_from_database(hash_z)
        HROM_utility.AppendHRomWeightsToRomParameters()
        #HROM_utility.CreateHRomModelParts() TODO Fix HROM model part creation. This is broken for generic model parts and MED files
        self.GenerateDatabaseSummary()

    def _LocalEmpiricalCubature(self, residuals_list, initial_candidates=None):
        """
        Runs the ECM for each matrix of projected residuals (one per solver), with the HROM settings of that solver. The entities
        selected so far are the initial candidates of the next one. Returns the selected indexes and their weights (one column per coupled solver).
        """
        hrom_settings_list = [coupled_solver["HROM"] for coupled_solver in self._GetCoupledSolvers()] or [self.general_rom_manager_parameters["HROM"]]
        weights_matrix = np.zeros((residuals_list[0].shape[0], len(residuals_list)))
        selected = np.empty(0, dtype=int)
        candidates = initial_candidates
        for i, residuals in enumerate(residuals_list):
            svd_truncation_tolerance = hrom_settings_list[i]["element_selection_svd_truncation_tolerance"].GetDouble()
            constrain_sum_of_weights = hrom_settings_list[i]["constraint_sum_weights"].GetBool()
            u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(residuals, svd_truncation_tolerance)
            element_selector = EmpiricalCubatureMethod()
            element_selector.SetUp(u, InitialCandidatesSet = candidates, constrain_sum_of_weights = constrain_sum_of_weights)
            element_selector.Run()
            if not element_selector.success:
                KratosMultiphysics.Logger.PrintWarning("RomManager", f"The Empirical Cubature Method did not converge for solver {i} using the initial set of candidates. Launching again without initial candidates.")
                element_selector.SetUp(u, InitialCandidatesSet = None, constrain_sum_of_weights = constrain_sum_of_weights)
                element_selector.Run()
            z_i = np.atleast_1d(np.squeeze(element_selector.z)).astype(int)
            weights_matrix[z_i, i] = np.squeeze(element_selector.w)
            selected = np.union1d(selected, z_i)
            candidates = selected
        return selected, weights_matrix[selected, :]

    def _GetResidualProjected(self):
        """
        Fetches the projected residuals from the corresponding output folders.
        In the case of coupled physics, it creates a consolidated vertical matrix.
        """

        def load_and_block_residuals(folder_path):
            """Robustly loads all .npy files in a folder in the correct chronological order."""
            path_obj = Path(folder_path)
            if not path_obj.exists():
                raise FileNotFoundError(f"Expected residual folder not found: {folder_path}")

            # Natural sort ensures 'Residual_2' comes before 'Residual_10'
            def natural_sort_key(s):
                return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', str(s))]

            files = sorted(path_obj.glob("*.npy"), key=lambda x: natural_sort_key(x.name))

            if not files:
                raise ValueError(f"No .npy files found in {folder_path}")

            # Load and consolidate into a single block
            residuals_list = [np.load(f) for f in files]
            return np.block(residuals_list)

        # Vertical stack of the residuals of each solver (one folder per coupled sub-solver, or a single one)
        return np.vstack([load_and_block_residuals(folder) for _, _, folder in self._GetResidualsProjectedOutputSettings()])


    def _LaunchHROM(self, mu_train, nn_rom_interface=None,gid_and_vtk_name ='HROM_Fit'):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        for Id, mu in enumerate(mu_train):
            in_database, _ = self.data_base.check_if_in_database("HROM", mu)
            if not in_database:
                parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
                parameters_copy = self._AddNumpyOutputToProjectParameters(parameters_copy)
                parameters_copy = self._StoreResultsByName(parameters_copy,gid_and_vtk_name,mu,Id)
                materials_file_name = self._GetMaterialsFileName(parameters_copy)
                self.UpdateMaterialParametersFile(materials_file_name, mu)
                model = KratosMultiphysics.Model()
                analysis_stage_class = type(self._SetUpRomSimulationInstance(model, parameters_copy, nn_rom_interface=nn_rom_interface))
                simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)
                simulation.Run()
                self.data_base.add_to_database("QoI_HROM", mu, simulation.GetFinalData())
                SnapshotsMatrix = self._GetSnapshotsMatrixFromNumpyOutput()
                self.data_base.add_to_database("HROM", mu, SnapshotsMatrix)

        self.GenerateDatabaseSummary()


    def _LaunchRunFOM(self, mu_run):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        for Id, mu in enumerate(mu_run):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._StoreResultsByName(parameters_copy,'FOM_Run',mu,Id)
            materials_file_name = self._GetMaterialsFileName(parameters_copy)
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters_copy)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)
            simulation.Run()
            self.QoI_Run_FOM.append(simulation.GetFinalData())

    def _LaunchRunROM(self, mu_run, nn_rom_interface=None):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        for Id, mu in enumerate(mu_run):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._StoreResultsByName(parameters_copy,'ROM_Run',mu,Id)
            materials_file_name = self._GetMaterialsFileName(parameters_copy)
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(self._SetUpRomSimulationInstance(model, parameters_copy, nn_rom_interface=nn_rom_interface))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)
            simulation.Run()
            self.QoI_Run_ROM.append(simulation.GetFinalData())


    def _LaunchRunHROM(self, mu_run, use_full_model_part, nn_rom_interface):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        if not use_full_model_part:
            model_import_settings = parameters["solver_settings"]["model_import_settings"]
            if model_import_settings.Has("input_filename"):
                model_part_name = model_import_settings["input_filename"].GetString()
                model_import_settings["input_filename"].SetString(f"{model_part_name}HROM")
            else:
                KratosMultiphysics.Logger.PrintWarning("RomManager", "No 'input_filename' in 'model_import_settings'. Running the HROM on the full model part.")

        for Id, mu in enumerate(mu_run):
            parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
            parameters_copy = self._StoreResultsByName(parameters_copy,'HROM_Run',mu,Id)
            materials_file_name = self._GetMaterialsFileName(parameters_copy)
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(self._SetUpRomSimulationInstance(model, parameters_copy, nn_rom_interface=nn_rom_interface))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)
            simulation.Run()
            self.QoI_Run_HROM.append(simulation.GetFinalData())

    def _LaunchTrainNeuralNetwork(self, mu_train, mu_validation):
        RomNeuralNetworkTrainer = self._TryImportNNTrainer()
        rom_nn_trainer = RomNeuralNetworkTrainer(self.general_rom_manager_parameters, mu_train, mu_validation, self.data_base)
        rom_nn_trainer.TrainNetwork()
        self.data_base.add_to_database("Neural_Network", mu_train , None)
        rom_nn_trainer.EvaluateNetwork()


    def _LaunchTestNeuralNetworkReconstruction(self,mu_train, mu_validation):
        RomNeuralNetworkTrainer = self._TryImportNNTrainer()
        rom_nn_trainer = RomNeuralNetworkTrainer(self.general_rom_manager_parameters, mu_train, mu_validation, self.data_base)
        rom_nn_trainer.EvaluateNetwork()

    def InitializeDummySimulationForSnapshotsModelPart(self):
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters = self._StoreNoResults(parameters)
        model = KratosMultiphysics.Model()
        analysis_stage_class = self._GetAnalysisStageClass(parameters)
        simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters,[])
        simulation.Initialize()
        return model[self.general_rom_manager_parameters["ROM"]["model_part_name"].GetString()]


    def InitializeDummySimulationForHromTrainingUtility(self, nn_rom_interface=None):
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters = self._StoreNoResults(parameters)
        model = KratosMultiphysics.Model()
        analysis_stage_class = type(self._SetUpRomSimulationInstance(model, parameters, nn_rom_interface=nn_rom_interface))
        simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters,[])
        simulation.Initialize()
        return simulation.GetHROM_utility()


    def InitializeDummySimulationForPetrovGalerkinTrainingUtility(self):
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters = self._StoreNoResults(parameters)
        model = KratosMultiphysics.Model()
        analysis_stage_class = type(self._SetUpRomSimulationInstance(model, parameters))
        simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters,[])
        simulation.Initialize()
        return simulation.GetPetrovGalerkinTrainUtility()


    def _AddHromParametersToRomParameters(self,f):
        f["hrom_settings"]["hrom_format"] = self.general_rom_manager_parameters["HROM"]["hrom_format"].GetString()
        f["hrom_settings"]["element_selection_type"] = self.general_rom_manager_parameters["HROM"]["element_selection_type"].GetString()
        f["hrom_settings"]["element_selection_svd_truncation_tolerance"] = self.general_rom_manager_parameters["HROM"]["element_selection_svd_truncation_tolerance"].GetDouble()
        f["hrom_settings"]["constraint_sum_weights"] = self.general_rom_manager_parameters["HROM"]["constraint_sum_weights"].GetBool()
        f["hrom_settings"]["svd_type"] = self.general_rom_manager_parameters["HROM"]["svd_type"].GetString()
        f["hrom_settings"]["create_hrom_visualization_model_part"] = self.general_rom_manager_parameters["HROM"]["create_hrom_visualization_model_part"].GetBool()
        f["hrom_settings"]["include_elements_model_parts_list"] = self.general_rom_manager_parameters["HROM"]["include_elements_model_parts_list"].GetStringArray()
        f["hrom_settings"]["include_conditions_model_parts_list"] = self.general_rom_manager_parameters["HROM"]["include_conditions_model_parts_list"].GetStringArray()
        f["hrom_settings"]["initial_candidate_elements_model_part_list"] = self.general_rom_manager_parameters["HROM"]["initial_candidate_elements_model_part_list"].GetStringArray()
        f["hrom_settings"]["initial_candidate_conditions_model_part_list"] = self.general_rom_manager_parameters["HROM"]["initial_candidate_conditions_model_part_list"].GetStringArray()
        f["hrom_settings"]["include_nodal_neighbouring_elements_model_parts_list"] = self.general_rom_manager_parameters["HROM"]["include_nodal_neighbouring_elements_model_parts_list"].GetStringArray()
        f["hrom_settings"]["include_minimum_condition"] = self.general_rom_manager_parameters["HROM"]["include_minimum_condition"].GetBool()
        f["hrom_settings"]["include_condition_parents"] = self.general_rom_manager_parameters["HROM"]["include_condition_parents"].GetBool()
        f["hrom_settings"]["echo_level"] = self.general_rom_manager_parameters["HROM"]["echo_level"].GetInt()

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
            # The order of the coupled solvers defines the HROM weights column (weight_vector_index) used by each of them
            coupled_solvers = self._GetCoupledSolverNames()
            if coupled_solvers:
                f['coupled_solvers'] = coupled_solvers
            else:
                f.pop('coupled_solvers', None)
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
            elif simulation_to_run=='GalerkinROM_ANN':
                f['train_hrom']=False
                f['run_hrom']=False
                f['projection_strategy']="galerkin_ann"
                f["rom_settings"]['rom_bns_settings'] = self._SetGalerkinBnSParameters()
            elif simulation_to_run=="trainHROMGalerkin_ANN":
                f['projection_strategy']="galerkin_ann"
                f['train_hrom']=True
                f['run_hrom']=False
                f["rom_settings"]['rom_bns_settings'] = self._SetGalerkinBnSParameters()
            elif simulation_to_run=="runHROMGalerkin_ANN":
                f['projection_strategy']="galerkin_ann"
                f['train_hrom']=False
                f['run_hrom']=True
                f["rom_settings"]['rom_bns_settings'] = self._SetGalerkinBnSParameters()
            elif simulation_to_run=='lspg_ANN':
                f['train_hrom']=False
                f['run_hrom']=False
                f['projection_strategy']="lspg_ann"
                f["rom_settings"]['rom_bns_settings'] = self._SetLSPGBnSParameters()
            elif simulation_to_run=="trainHROMlspg_ANN":
                f['projection_strategy']="lspg_ann"
                f['train_hrom']=True
                f['run_hrom']=False
                f["rom_settings"]['rom_bns_settings'] = self._SetLSPGBnSParameters()
            elif simulation_to_run=="runHROMlspg_ANN":
                f['projection_strategy']="lspg_ann"
                f['train_hrom']=False
                f['run_hrom']=True
                f["rom_settings"]['rom_bns_settings'] = self._SetLSPGBnSParameters()
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

    def _SetUpRomSimulationInstance(self, model, parameters, nn_rom_interface=None):
        rom_params = self.general_rom_manager_parameters["ROM"]
        return SetUpSimulationInstance(model, parameters, nn_rom_interface=nn_rom_interface,
                                       rom_basis_output_folder=rom_params["rom_basis_output_folder"].GetString(),
                                       rom_basis_output_name=rom_params["rom_basis_output_name"].GetString())


    def _GetSnapshotsOutputPath(self):
        return Path(self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].GetString()) / "numpy_snapshots"


    def _GetMaterialsFileName(self, parameters):
        """Returns the materials file name, or None if it is not defined at the solver_settings level (e.g. coupled solvers)."""
        solver_settings = parameters["solver_settings"]
        if solver_settings.Has("material_import_settings"):
            return solver_settings["material_import_settings"]["materials_filename"].GetString()
        return None

    # Settings that can be given per coupled solver. The rest of the 'ROM' and 'HROM' settings are shared by all of them
    _COUPLED_SOLVER_SETTINGS = {
        "ROM" : ["model_part_name", "nodal_unknowns", "svd_truncation_tolerance"],
        "HROM" : ["element_selection_svd_truncation_tolerance", "constraint_sum_weights"]
    }

    def _SetUpCoupledSolvers(self):
        """
        Validates the 'coupled_solvers' and completes the 'ROM' and 'HROM' settings of each of them with the general ones.
        A single (monolithic) ROM basis shared by all the coupled solvers is currently supported. Hence, its snapshots
        contain the 'nodal_unknowns' of all of them and the smallest 'svd_truncation_tolerance' among them is used.
        """
        coupled_solvers = self.general_rom_manager_parameters["coupled_solvers"]
        for i in range(coupled_solvers.size()):
            self._CompleteCoupledSolverSettings(coupled_solvers[i])
        if coupled_solvers.size() == 0:
            return

        rom_settings = self.general_rom_manager_parameters["ROM"]
        tolerances = {coupled_solvers[i]["sub_solver_name"].GetString() : coupled_solvers[i]["ROM"]["svd_truncation_tolerance"].GetDouble() for i in range(coupled_solvers.size())}
        if len(set(tolerances.values())) > 1:
            KratosMultiphysics.Logger.PrintWarning("RomManager", f"Different 'svd_truncation_tolerance' set for the coupled solvers {tolerances}. Only a monolithic ROM basis (shared by all the coupled solvers) is currently supported, so the smallest one ({min(tolerances.values())}) is used.")
        rom_settings["svd_truncation_tolerance"].SetDouble(min(tolerances.values()))
        nodal_unknowns = set()
        for i in range(coupled_solvers.size()):
            nodal_unknowns.update(coupled_solvers[i]["ROM"]["nodal_unknowns"].GetStringArray())
        rom_settings["nodal_unknowns"].SetStringArray(sorted(nodal_unknowns))

    def _CompleteCoupledSolverSettings(self, coupled_solver):
        coupled_solver.ValidateAndAssignDefaults(KratosMultiphysics.Parameters("""{
            "sub_solver_name": "",
            "ROM": {},
            "HROM": {}
        }"""))
        sub_solver_name = coupled_solver["sub_solver_name"].GetString()
        if not sub_solver_name:
            raise Exception("Each of the 'coupled_solvers' must provide its 'sub_solver_name' (e.g. 'fluid_solver').")
        for block, supported_keys in self._COUPLED_SOLVER_SETTINGS.items():
            general_settings = self.general_rom_manager_parameters[block]
            for key in coupled_solver[block].keys():
                if key not in supported_keys:
                    raise Exception(f"'{key}' cannot be set in the '{block}' settings of the coupled solver '{sub_solver_name}'. The ones supported per coupled solver are {supported_keys}; the rest are shared by all the coupled solvers and are set in the general '{block}' settings.")
            # Check the types against the general settings
            supported_defaults = KratosMultiphysics.Parameters()
            for key in supported_keys:
                supported_defaults.AddValue(key, general_settings[key])
            coupled_solver[block].ValidateDefaults(supported_defaults)
            for key in supported_keys:
                if not coupled_solver[block].Has(key):
                    coupled_solver[block].AddValue(key, general_settings[key])

    def _GetCoupledSolvers(self):
        """
        Returns the settings of the coupled sub-solvers ('sub_solver_name' and complete 'ROM' and 'HROM' settings).
        These are taken from 'coupled_solvers' or, if not provided, from the '<sub_solver_name>_settings' of a coupled solver
        in the project parameters. Empty for a single solver.
        """
        coupled_solvers = self.general_rom_manager_parameters["coupled_solvers"]
        if coupled_solvers.size() > 0:
            return [coupled_solvers[i] for i in range(coupled_solvers.size())]
        with open(self.project_parameters_name,'r') as parameter_file:
            solver_settings = KratosMultiphysics.Parameters(parameter_file.read())["solver_settings"]
        if solver_settings["solver_type"].GetString() not in ["ThermallyCoupled", "ThermoMechanicallyCoupled"]: # Coupled solvers supported by the RomAnalysis
            return []
        detected_solvers = []
        for key, value in solver_settings.items():
            if key.endswith("_solver_settings") and value.Has("model_part_name"):
                coupled_solver = KratosMultiphysics.Parameters("""{"ROM": {}, "HROM": {}}""")
                coupled_solver.AddString("sub_solver_name", key[:-len("_settings")])
                coupled_solver["ROM"].AddString("model_part_name", value["model_part_name"].GetString())
                self._CompleteCoupledSolverSettings(coupled_solver)
                detected_solvers.append(coupled_solver)
        return detected_solvers

    def _GetCoupledSolverNames(self):
        """Returns the names of the coupled sub-solvers (e.g. ['fluid_solver', 'thermal_solver']). Empty for a single solver."""
        return [coupled_solver["sub_solver_name"].GetString() for coupled_solver in self._GetCoupledSolvers()]

    def _GetResidualsProjectedOutputSettings(self):
        """Returns the (sub_solver_name, model_part_name, output folder) of each ProjectedResidualsOutputProcess: one per coupled sub-solver, or one for the solver."""
        rom_basis_output_folder = Path(self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].GetString())
        coupled_solvers = self._GetCoupledSolvers()
        if coupled_solvers:
            return [(coupled_solver["sub_solver_name"].GetString(), coupled_solver["ROM"]["model_part_name"].GetString(), rom_basis_output_folder / f"Residuals_{coupled_solver['sub_solver_name'].GetString()}") for coupled_solver in coupled_solvers]
        return [("", self.general_rom_manager_parameters["ROM"]["model_part_name"].GetString(), rom_basis_output_folder / "Residuals")]

    def _SetTrainHROMFlag(self, train_hrom):
        rom_parameters_file_path = (Path(self.general_rom_manager_parameters["ROM"]["rom_basis_output_folder"].GetString()) / self.general_rom_manager_parameters["ROM"]["rom_basis_output_name"].GetString()).with_suffix('.json')
        with rom_parameters_file_path.open('r+') as parameter_file:
            f = json.load(parameter_file)
            f['train_hrom'] = train_hrom
            parameter_file.seek(0)
            json.dump(f, parameter_file, indent=4)
            parameter_file.truncate()

    def _AddNumpyOutputToProjectParameters(self, parameters):
        rom_params = self.general_rom_manager_parameters["ROM"]
        numpy_output_parameters = KratosMultiphysics.Parameters("""{
            "python_module" : "numpy_output_process",
            "kratos_module" : "KratosMultiphysics.RomApplication",
            "process_name"  : "NumpyOutputProcess",
            "Parameters"    : {}
        }""")
        numpy_output_parameters["Parameters"].AddString("model_part_name", rom_params["model_part_name"].GetString())
        numpy_output_parameters["Parameters"].AddString("output_control_type", rom_params["snapshots_control_type"].GetString())
        numpy_output_parameters["Parameters"].AddDouble("output_interval", rom_params["snapshots_interval"].GetDouble())
        numpy_output_parameters["Parameters"].AddValue("nodal_results", rom_params["nodal_unknowns"])
        numpy_output_parameters["Parameters"].AddString("output_path", str(self._GetSnapshotsOutputPath()))
        parameters["output_processes"].AddEmptyArray("rom_manager_numpy_output")
        parameters["output_processes"]["rom_manager_numpy_output"].Append(numpy_output_parameters)

        # Remove leftovers from previous runs so that only the upcoming simulation snapshots are fetched
        shutil.rmtree(self._GetSnapshotsOutputPath(), ignore_errors=True)

        return parameters


    def _GetSnapshotsMatrixFromNumpyOutput(self):
        snapshots_output_path = self._GetSnapshotsOutputPath()
        # Files are named 'solution_<step>' or 'solution_<time>', sort them numerically
        snapshot_files = sorted(snapshots_output_path.glob("solution_*.npy"), key=lambda f: float(f.stem[len("solution_"):]))
        snapshots_matrix = np.block([np.load(f) for f in snapshot_files])
        shutil.rmtree(snapshots_output_path)

        return snapshots_matrix


    def _ComputeSVD(self, snapshots_matrix):
        # Calculate the randomized SVD of the snapshots matrix
        svd_truncation_tolerance = self.general_rom_manager_parameters["ROM"]["svd_truncation_tolerance"].GetDouble()
        u,sigma,_,_= RandomizedSingularValueDecomposition().Calculate(snapshots_matrix, svd_truncation_tolerance)
        return u, sigma


    def _PrintRomBasis(self, u, sigma):
        rom_params = self.general_rom_manager_parameters["ROM"]
        rom_basis_output_format = rom_params["rom_basis_output_format"].GetString()
        rom_basis_output_folder = Path(rom_params["rom_basis_output_folder"].GetString())
        rom_basis_output_name = rom_params["rom_basis_output_name"].GetString()
        print_singular_values = rom_params["print_singular_values"].GetBool()
        if print_singular_values and rom_basis_output_format == "json":
            err_msg = 'Cannot print singular values if using the "json" output format. Please use "numpy" instead.'
            raise Exception(err_msg)

        # Note that the nodal unknowns are sorted alphabetically, consistently with the NumpyOutputProcess
        nodal_unknowns = sorted(rom_params["nodal_unknowns"].GetStringArray())
        n_nodal_unknowns = len(nodal_unknowns)
        model_part = self.InitializeDummySimulationForSnapshotsModelPart()

        # Initialize the Python dictionary with the default settings
        rom_basis_dict = {
            "rom_manager" : True,
            "train_hrom": False,
            "run_hrom": False,
            "projection_strategy": "galerkin",
            "assembling_strategy": "global",
            "rom_format": rom_basis_output_format,
            "rom_settings": {
                "rom_bns_settings": {},
                "nodal_unknowns": nodal_unknowns,
                "number_of_rom_dofs": np.shape(u)[1],
                "petrov_galerkin_number_of_rom_dofs": 0
            },
            "hrom_settings": {
                "hrom_format": rom_basis_output_format
            },
            "nodal_modes": {},
            "elements_and_weights" : {}
        }

        # Create the folder if it doesn't already exist
        if not rom_basis_output_folder.exists():
            rom_basis_output_folder.mkdir(parents=True)

        if rom_basis_output_format == "json":
            # Storing modes in JSON format
            i = 0
            for node in model_part.Nodes:
                rom_basis_dict["nodal_modes"][node.Id] = u[i:i+n_nodal_unknowns].tolist()
                i += n_nodal_unknowns
        elif rom_basis_output_format == "numpy":
            # Storing modes in Numpy format
            node_ids = np.array([node.Id for node in model_part.Nodes])
            np.save(rom_basis_output_folder / "RightBasisMatrix.npy", u)
            np.save(rom_basis_output_folder / "NodeIds.npy", node_ids)
            if print_singular_values:
                np.save(rom_basis_output_folder / "SingularValuesVector.npy", sigma)
        else:
            err_msg = f"Unsupported output format {rom_basis_output_format}. Available options are 'json' and 'numpy'."
            raise Exception(err_msg)

        # Creating the ROM JSON file containing or not the modes depending on the output format
        output_filename = rom_basis_output_folder / f"{rom_basis_output_name}.json"
        with output_filename.open('w') as f:
            json.dump(rom_basis_dict, f, indent = 4)


    def _AddResidualsProjectedOutputProcessToProjectParameters(self, parameters):
        if not parameters.Has("output_processes"):
            parameters.AddEmptyValue("output_processes")

        # Remove existing to avoid duplicate process
        if parameters["output_processes"].Has("rom_residuals_output"):
            parameters["output_processes"].RemoveValue("rom_residuals_output")

        parameters["output_processes"].AddEmptyArray("rom_residuals_output")
        # One process per coupled sub-solver, or a single one
        for sub_solver_name, model_part_name, _ in self._GetResidualsProjectedOutputSettings():
            process_parameters = self._SetUpProjectedResidualsOutputProcessParameters(model_part_name, sub_solver_name)
            parameters["output_processes"]["rom_residuals_output"].Append(process_parameters)

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



    def _SetUpRomManagerParameters(self, input):

        if input is None:
            input = KratosMultiphysics.Parameters()
        self.general_rom_manager_parameters = input


        default_settings = KratosMultiphysics.Parameters("""{
            "rom_stages_to_train" : ["ROM","HROM"],             // ["ROM","HROM"]
            "rom_stages_to_test" : [],              // ["ROM","HROM"]
            "paralellism" : null,                        // null, TODO: add "compss"
            "projection_strategy": "galerkin",            // "lspg", "galerkin", "petrov_galerkin"
            "type_of_decoder" : "linear",               // "linear" "ann_enhanced",  TODO: add "quadratic"
            "assembling_strategy": "global",            // "global", "elemental"
            "save_gid_output": false,                    // false, true #if true, it must exits previously in the ProjectParameters.json
            "save_vtk_output": false,                    // false, true #if true, it must exits previously in the ProjectParameters.json
            "output_name": "id",                         // "id" , "mu"
            "store_nonconverged_fom_solutions": false,
            "coupled_solvers": [],                       // [] for a single solver. e.g. [{"sub_solver_name": "fluid_solver", "ROM": {"model_part_name": "FluidModelPart"}, "HROM": {}}, {"sub_solver_name": "thermal_solver", "ROM": {"model_part_name": "ThermalModelPart"}, "HROM": {}}]
            "ROM":{
                "svd_truncation_tolerance": 1e-5,
                "model_part_name": "Structure",                            // This changes depending on the simulation: Structure, FluidModelPart, ThermalPart #TODO: Idenfity it automatically
                "nodal_unknowns": ["DISPLACEMENT_X","DISPLACEMENT_Y"],     // Main unknowns. Snapshots are taken from these
                "rom_basis_output_format": "numpy",
                "rom_basis_output_name": "RomParameters",
                "rom_basis_output_folder": "rom_data",
                "snapshots_control_type": "step",                          // "step", "time"
                "snapshots_interval": 1,
                "print_singular_values": false,
                "use_non_converged_sols" : false,
                "galerkin_rom_bns_settings": {
                    "monotonicity_preserving": false
                },
                "lspg_rom_bns_settings": {
                    "train_petrov_galerkin": false,
                    "basis_strategy": "residuals",                        // 'residuals', 'jacobian', 'reactions'
                    "include_phi": false,
                    "svd_truncation_tolerance": 1e-12,
                    "solving_technique": "normal_equations",              // 'normal_equations', 'qr_decomposition'
                    "monotonicity_preserving": false
                },
                "petrov_galerkin_rom_bns_settings": {
                    "monotonicity_preserving": false
                },
                "ann_enhanced_settings": {
                    "modes":[5,50],
                    "layers_size":[200,200],
                    "batch_size":2,
                    "epochs":800,
                    "NN_gradient_regularisation_weight": 0.0,
                    "lr_strategy":{
                        "scheduler": "sgdr",
                        "base_lr": 0.001,
                        "additional_params": [1e-4, 10, 400]
                    },
                    "training":{
                        "retrain_if_exists" : false  // If false only one model will be trained for each the mu_train and NN hyperparameters combination
                    },
                    "online":{
                        "model_number": 0   // out of the models existing for the same parameters, this is the model that will be lauched
                    }
                }
            },
            "HROM":{
                "hrom_format": "numpy",                            //  "json", "numpy"
                "element_selection_type": "empirical_cubature",
                "element_selection_svd_truncation_tolerance": 1.0e-6,
                "constraint_sum_weights": true,                   // if true, then sum(w) = num_elems (this avoids trivial solutions sum(w)=0)
                "svd_type": "numpy_rsvd",                         //  "numpy_svd", "numpy_rsvd"
                "create_hrom_visualization_model_part" : true,
                "include_elements_model_parts_list": [],          //The elements of the submodel parts included in this list will be considered in the HROM model part
                "include_conditions_model_parts_list": [],         //The conditions of the submodel parts included in this list will be considered in the HROM model part
                "initial_candidate_elements_model_part_list" : [],        //These elements will be given priority when creating the HROM model part
                "initial_candidate_conditions_model_part_list" : [],       //These conditions will be given priority when creating the HROM model part
                "include_nodal_neighbouring_elements_model_parts_list":[],
                "include_minimum_condition": false,       // if true, keep at least one condition per submodelpart
                "include_condition_parents": false,       // if true, when a condition is chosen by the ECM algorithm, the parent element is also included
                "echo_level" : 0                          // if >0, get hrom training status promts
            }
        }""")

        self.general_rom_manager_parameters.RecursivelyValidateAndAssignDefaults(default_settings)
        self._SetUpCoupledSolvers()




    def DefaultCustomizeSimulation(self, cls, global_model, parameters, mu=None):
        # Default function that does nothing special
        class DefaultCustomSimulation(cls):
            def __init__(self, model, project_parameters):
                super().__init__(model, project_parameters)

            def Initialize(self):
                super().Initialize()

            def FinalizeSolutionStep(self):
                super().FinalizeSolutionStep()

            def CustomMethod(self):
                pass  # Do nothing special

        return DefaultCustomSimulation(global_model, parameters)


    def ActivateNonconvergedSolutionsGathering(self, simulation):

        # Patch the RomAnalysis class to save the selected time steps results
        def Initialize(cls):
            super(type(simulation), cls).Initialize()
            cls._GetSolver()._GetSolutionStrategy().SetUpNonconvergedSolutionsFlag(True) # this assumes the strategy used contains this method

        def GetNonconvergedSolutions(cls):
            a,_ = cls._GetSolver()._GetSolutionStrategy().GetNonconvergedSolutions()
            return np.asarray(a)

        simulation.Initialize  = types.MethodType(Initialize, simulation)
        simulation.GetNonconvergedSolutions  = types.MethodType(GetNonconvergedSolutions, simulation)

        return simulation




    def DefaultUpdateProjectParameters(self, parameters, mu=None):
        return parameters

    def DefaultUpdateMaterialParametersFile(self, material_parametrs_file_name=None, mu=None):
        pass
        # with open(material_parametrs_file_name, mode="r+") as f:
        #     data = json.load(f)
        #     #change the angles of 1st and 2nd layer
        #     data["properties"][0]["Material"]["Variables"]["EULER_ANGLES"][0] = mu[0]
        #     data["properties"][1]["Material"]["Variables"]["EULER_ANGLES"][0] = mu[1]
        #     #write to file and save file
        #     f.seek(0)
        #     json.dump(data, f, indent=4)
        #     f.truncate()

    def _SetUpProjectedResidualsOutputProcessParameters(self, model_part_name, sub_solver_name=""):
        process_settings = self._GetDefaulProjectedResidualsOutputProcessParameters()
        rom_params = self.general_rom_manager_parameters["ROM"]

        # Set the specific solver and model part names
        process_settings["Parameters"]["model_part_name"].SetString(model_part_name)
        if sub_solver_name:
            process_settings["Parameters"]["sub_solver_name"].SetString(sub_solver_name)

        # Map global ROM settings
        if rom_params.Has("snapshots_interval"):
            process_settings["Parameters"]["output_interval"] = rom_params["snapshots_interval"].Clone()

        if rom_params.Has("snapshots_control_type"):
            process_settings["Parameters"]["output_control_type"] = rom_params["snapshots_control_type"].Clone()

        # Dynamically set output path to prevent conflicts (e.g., rom_data/Residuals_fluid_solver)
        base_folder = "rom_data"
        if rom_params.Has("rom_basis_output_folder"):
            base_folder = rom_params["rom_basis_output_folder"].GetString()

        folder_suffix = f"_{sub_solver_name}" if sub_solver_name else ""
        process_settings["Parameters"]["output_path"].SetString(f"{base_folder}/Residuals{folder_suffix}")

        return process_settings





    def _GetAnalysisStageClass(self, parameters):

        analysis_stage_module_name = parameters["analysis_stage"].GetString()
        analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
        analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

        analysis_stage_module = importlib.import_module(analysis_stage_module_name)
        analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

        return analysis_stage_class


    def _GetDefaulProjectedResidualsOutputProcessParameters(self):
        return KratosMultiphysics.Parameters("""{
            "python_module": "projected_residuals_output_process",
            "kratos_module": "KratosMultiphysics.RomApplication",
            "process_name": "ProjectedResidualsOutputProcess",
            "Parameters": {
                "model_part_name": "",
                "output_control_type": "step",
                "output_interval": 1,
                "output_path": "rom_data/Residuals",
                "sub_solver_name": "",
                "range_of_entity_ids_to_fetch_residuals_projected": ["1","end"]
            }
        }""")



    def SetUpQuantityOfInterestContainers(self):
        #TODO implement more options if the QoI is too large to keep in RAM
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


    def GenerateDatabaseSummary(self):
        self.data_base.generate_database_summary()

    def GenerateDatabaseCompleteDump(self):
        self.data_base.dump_database_as_excel()

    def _TryImportNNTrainer(self):
        try:
            from KratosMultiphysics.RomApplication.rom_nn_trainer import RomNeuralNetworkTrainer
            return RomNeuralNetworkTrainer
        except ImportError:
            err_msg = f'Failed to import the RomNeuralNetworkTrainer class. Make sure TensorFlow is properly installed.'
            raise Exception(err_msg)
