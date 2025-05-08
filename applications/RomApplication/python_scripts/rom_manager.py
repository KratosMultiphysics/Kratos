import numpy as np
import importlib
import json
from pathlib import Path
import types
import KratosMultiphysics
from KratosMultiphysics.RomApplication.rom_database import RomDatabase
from KratosMultiphysics.RomApplication.rom_testing_utilities import SetUpSimulationInstance
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
from KratosMultiphysics.RomApplication.rom_nn_interface import NN_ROM_Interface
from KratosMultiphysics.RomApplication.initialize_from_snapshot_process import InitializeFromSnapshotProcess


class RomManager(object):

    def __init__(self,project_parameters_name="ProjectParameters.json", general_rom_manager_parameters=None, CustomizeSimulation=None, UpdateProjectParameters=None,UpdateMaterialParametersFile=None, mu_names=None):
        #FIXME:
        # - Use a method (upcoming) for smothly retrieving solutions. In here we are using the RomBasisOutput process in order to store the solutions
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

    # def Fit(self, mu_train=[None],mu_validation=[None], mu_init_train=[None], mu_init_validation=[None]):
    def Fit(self, mu_train=[None],mu_validation=[None], start_from_closest_mu=False, filter_nan = False):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
        type_of_decoder = self.general_rom_manager_parameters["type_of_decoder"].GetString()
        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            if type_of_decoder =="ann_enhanced":
                if any(item == "ROM" for item in training_stages):
                    self._LaunchTrainROM(mu_train, start_from_closest_mu=start_from_closest_mu)
                    self._LaunchFOM(mu_validation, start_from_closest_mu=start_from_closest_mu) #What to do here with the gid and vtk results?
                    self.TrainAnnEnhancedROM(mu_train,mu_validation)
                    self._ChangeRomFlags(simulation_to_run = "GalerkinROM_ANN")
                    nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base)
                    self._LaunchROM(mu_train, mu_train = mu_train, nn_rom_interface=nn_rom_interface)
                if any(item == "HROM" for item in training_stages):
                    err_msg = f'HROM is not available yet for ann_enhanced decoders.'
                    raise Exception(err_msg)
            elif type_of_decoder =="linear":
                if any(item == "ROM" for item in training_stages):
                    self._LaunchTrainROM(mu_train, start_from_closest_mu=start_from_closest_mu)
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
                    self._LaunchTrainROM(mu_train, start_from_closest_mu=start_from_closest_mu)
                    self._LaunchFOM(mu_validation, start_from_closest_mu=start_from_closest_mu) #What to do here with the gid and vtk results?
                    self.TrainAnnEnhancedROM(mu_train,mu_validation)
                    self._ChangeRomFlags(simulation_to_run = "lspg_ANN")
                    nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base)
                    self._LaunchROM(mu_train, mu_train = mu_train, nn_rom_interface=nn_rom_interface)
                if any(item == "HROM" for item in training_stages):
                    err_msg = f'HROM is not available yet for ann_enhanced decoders.'
                    raise Exception(err_msg)
            elif type_of_decoder =="linear":
                if any(item == "ROM" for item in training_stages):
                    self._LaunchTrainROM(mu_train, start_from_closest_mu=start_from_closest_mu)
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
                    self._LaunchTrainROM(mu_train, start_from_closest_mu=start_from_closest_mu)
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
        if type_of_decoder =="ann_enhanced":
            if nn_rom_interface.get_using_manual_model_path():
                self.ComputeErrors(mu_train, case='Test', filter_nan=filter_nan, mu_aux=None, nn_model_name=nn_rom_interface.get_nn_model_name())
            else:
                self.ComputeErrors(mu_train, case='Test', filter_nan=filter_nan, mu_aux=mu_train)
        else:
            self.ComputeErrors(mu_train, case='Test', filter_nan=filter_nan)

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

    def TestNeuralNetworkReconstruction(self, mu_train, mu_validation, mu_test):
        custom_model_path = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["online"]["custom_model_path"].GetString()
        if custom_model_path == '':
            custom_model_path = None
        self._LaunchTestNeuralNetworkReconstruction( mu_train, mu_validation, mu_test, custom_model_path)


    def Test(self, mu_test=[None], mu_train=[None], start_from_closest_mu=False, filter_nan=False):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
        type_of_decoder = self.general_rom_manager_parameters["type_of_decoder"].GetString()
        custom_model_path = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["online"]["custom_model_path"].GetString()
        if custom_model_path == '':
            custom_model_path = None

        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            if type_of_decoder =="ann_enhanced":
                if any(item == "ROM" for item in testing_stages):
                    self._LoadSolutionBasis(mu_train)
                    self._LaunchFOM(mu_test, start_from_closest_mu=start_from_closest_mu, gid_and_vtk_name='FOM_Test')
                    self._ChangeRomFlags(simulation_to_run = "GalerkinROM_ANN")
                    nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base, custom_model_path=custom_model_path)
                    self._LaunchROM(mu_test, mu_train = mu_train, gid_and_vtk_name='ROM_Test', nn_rom_interface=nn_rom_interface)
                if any(item == "HROM" for item in testing_stages):
                    err_msg = f'HROM is not available yet for ann_enhanced decoders.'
                    raise Exception(err_msg)
            elif type_of_decoder =="linear":
                if any(item == "ROM" for item in testing_stages):
                    self._LoadSolutionBasis(mu_train)
                    self._LaunchFOM(mu_test, start_from_closest_mu=start_from_closest_mu, gid_and_vtk_name='FOM_Test')
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
                if any(item == "ROM" for item in testing_stages):
                    self._LoadSolutionBasis(mu_train)
                    self._LaunchFOM(mu_test, start_from_closest_mu=start_from_closest_mu, gid_and_vtk_name='FOM_Test')
                    self._ChangeRomFlags(simulation_to_run = "lspg_ANN")
                    nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base, custom_model_path=custom_model_path)
                    self._LaunchROM(mu_test, mu_train = mu_train, gid_and_vtk_name='ROM_Test', nn_rom_interface=nn_rom_interface)
                if any(item == "HROM" for item in testing_stages):
                    err_msg = f'HROM is not available yet for ann_enhanced decoders.'
            elif type_of_decoder =="linear":
                if any(item == "ROM" for item in testing_stages):
                    self._LoadSolutionBasis(mu_train)
                    self._LaunchFOM(mu_test, start_from_closest_mu=start_from_closest_mu, gid_and_vtk_name='FOM_Test')
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
                    self._LaunchFOM(mu_test, start_from_closest_mu=start_from_closest_mu, gid_and_vtk_name='FOM_Test')
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
        if type_of_decoder =="ann_enhanced":
            if nn_rom_interface.get_using_manual_model_path():
                self.ComputeErrors(mu_test, case='Test', filter_nan=filter_nan, mu_aux=None, nn_model_name=nn_rom_interface.get_nn_model_name())
            else:
                self.ComputeErrors(mu_test, case='Test', filter_nan=filter_nan, mu_aux=mu_train)
        else:
            self.ComputeErrors(mu_test, case='Test', filter_nan=filter_nan)


    def RunFOM(self, mu_run=[None]):
        self._LaunchRunFOM(mu_run)


    def RunROM(self, mu_run=[None], mu_train=[None]):
        chosen_projection_strategy = self.general_rom_manager_parameters["projection_strategy"].GetString()
        type_of_decoder = self.general_rom_manager_parameters["type_of_decoder"].GetString()
        custom_model_path = self.general_rom_manager_parameters["ROM"]["ann_enhanced_settings"]["online"]["custom_model_path"].GetString()
        if custom_model_path == "":
            custom_model_path = None

        nn_rom_interface = None
        self._LoadSolutionBasis(mu_train)
        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            if type_of_decoder =="ann_enhanced":
                self._ChangeRomFlags(simulation_to_run = "GalerkinROM_ANN")
                nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base, custom_model_path=custom_model_path)
            elif type_of_decoder =="linear":
                self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif chosen_projection_strategy == "lspg":
            if type_of_decoder =="ann_enhanced":
                self._ChangeRomFlags(simulation_to_run = "lspg_ANN")
                nn_rom_interface = NN_ROM_Interface(mu_train, self.data_base, custom_model_path=custom_model_path)
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
        self._LoadSolutionBasis(mu_train)
        #######################
        ######  Galerkin ######
        if chosen_projection_strategy == "galerkin":
            if type_of_decoder =="ann_enhanced":
                err_msg = f'HROM only supports linear projection strategy for the moment'
                raise Exception(err_msg)
            elif type_of_decoder =="linear":
                self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif chosen_projection_strategy == "lspg":
            if type_of_decoder =="ann_enhanced":
                err_msg = f'HROM only supports linear projection strategy for the moment'
                raise Exception(err_msg)
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
        self._LaunchRunHROM(mu_run, use_full_model_part)



    def ComputeErrors(self, mu_list, case="Fit", filter_nan = False, mu_aux = None, nn_model_name = None):
        fom_snapshots = self.data_base.get_snapshots_matrix_from_database(mu_list, mu_aux = mu_aux, table_name=f'FOM', nn_model_name=nn_model_name)
        if case=="Fit":
            stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
        elif case=="Test":
            stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
        stages = {"ROM", "HROM"} & set(stages)  # Ensures only "ROM" or "HROM" if present
        rom_snapshots = None
        hrom_snapshots = None
        if "ROM" in stages:
            rom_snapshots = self.data_base.get_snapshots_matrix_from_database(mu_list, mu_aux = mu_aux, table_name=f'ROM', nn_model_name=nn_model_name)
            nan_cols_ROM = np.unique(np.argwhere(np.isnan(rom_snapshots))[:,1])
            self.ROM_non_converged_samples[case]=nan_cols_ROM.shape[0]
            filter_ROM = nan_cols_ROM if filter_nan else None
            self.ROMvsFOM[case] = self._GetGeometricMeanError(rom_snapshots, fom_snapshots, filter = filter_ROM)
            # error_rom_fom = np.linalg.norm(fom_snapshots - rom_snapshots) / np.linalg.norm(fom_snapshots)
            # self.ROMvsFOM[case] = error_rom_fom
        if "HROM" in stages:
            if rom_snapshots is None:  # Only fetch if not already fetched
                rom_snapshots = self.data_base.get_snapshots_matrix_from_database(mu_list, mu_aux = mu_aux, table_name=f'ROM', nn_model_name=nn_model_name)
                nan_cols_ROM = np.unique(np.argwhere(np.isnan(rom_snapshots))[:,1])
                self.ROM_non_converged_samples[case]=nan_cols_ROM.shape[0]
                filter_ROM = nan_cols_ROM if filter_nan else None
            hrom_snapshots = self.data_base.get_snapshots_matrix_from_database(mu_list, mu_aux = mu_aux, table_name=f'HROM', nn_model_name=nn_model_name)
            nan_cols_HROM = np.unique(np.argwhere(np.isnan(hrom_snapshots))[:,1])
            self.HROM_non_converged_samples[case]=nan_cols_HROM.shape[0]
            filter_HROM = nan_cols_HROM if filter_nan else None
            self.ROMvsHROM[case] = self._GetGeometricMeanError(hrom_snapshots, rom_snapshots, filter = filter_ROM) # This should be the union between filter_ROM nad filter_HROM
            self.FOMvsHROM[case] = self._GetGeometricMeanError(hrom_snapshots, fom_snapshots, filter = filter_HROM)
            # error_rom_hrom = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)
            # error_fom_hrom = np.linalg.norm(fom_snapshots - hrom_snapshots) / np.linalg.norm(fom_snapshots)
            # self.ROMvsHROM[case] = error_rom_hrom
            # self.FOMvsHROM[case] = error_fom_hrom

    def _GetGeometricMeanError(self, approx_snapshots, true_snapshots, filter = None):
        if filter is not None:
            approx_snapshots = np.delete(approx_snapshots, filter, axis=1)
            true_snapshots = np.delete(true_snapshots, filter, axis=1)

        N=true_snapshots.shape[1]
        err_numer=np.linalg.norm(true_snapshots-approx_snapshots, ord=2, axis=0)
        err_denom=np.linalg.norm(true_snapshots, ord=2, axis=0)
        return np.exp(np.sum(np.log(err_numer/err_denom))/N)
    

    def PrintErrors(self):
        training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
        testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()

        training_set = set(training_stages)
        testing_set = set(testing_stages)

        # Check in Fit
        if "ROM" in training_set:
            self.aux_print_errors(self.ROMvsFOM['Fit'], 'train', 'FOM vs ROM')
            self.aux_print_non_converged(self.ROM_non_converged_samples['Fit'], 'train', 'ROM')
        if "HROM" in training_set:
            self.aux_print_errors(self.ROMvsHROM['Fit'], 'train', 'ROM vs HROM')
            self.aux_print_errors(self.FOMvsHROM['Fit'], 'train', 'FOM vs HROM')
            self.aux_print_non_converged(self.ROM_non_converged_samples['Fit'], 'train', 'ROM')
            self.aux_print_non_converged(self.HROM_non_converged_samples['Fit'], 'train', 'HROM')

        # Check in Test
        if "ROM" in testing_set:
            self.aux_print_errors(self.ROMvsFOM['Test'], 'test', 'FOM vs ROM')
            self.aux_print_non_converged(self.ROM_non_converged_samples['Test'], 'test', 'ROM')
        if "HROM" in testing_set:
            self.aux_print_errors(self.ROMvsHROM['Test'], 'test', 'ROM vs HROM')
            self.aux_print_errors(self.FOMvsHROM['Test'], 'test', 'FOM vs HROM')
            self.aux_print_non_converged(self.ROM_non_converged_samples['Test'], 'test', 'ROM')
            self.aux_print_non_converged(self.HROM_non_converged_samples['Test'], 'test', 'HROM')

    def aux_print_errors(self, error, train_or_test, comparison_in_string):
        message = f"approximation error in {train_or_test} set {comparison_in_string}"
        if error is None:
            print(f"{message} not computed")
        else:
            print(f"{message}: {error}")

    def aux_print_non_converged(self, number_of_non_converged, train_or_test, snapshot_type):
        message = f"Number of non-converged samples in {train_or_test} set {snapshot_type}"
        if number_of_non_converged is None:
            print(f"{message}: not computed")
        else:
            print(f"{message}: {number_of_non_converged}")


    def SetupErrorsDictionaries(self):
        self.ROMvsFOM = {'Fit': None, 'Test': None}
        self.ROMvsHROM = {'Fit': None, 'Test': None}
        self.FOMvsHROM = {'Fit': None, 'Test': None}
        self.ROM_non_converged_samples = {'Fit': None, 'Test': None}
        self.HROM_non_converged_samples = {'Fit': None, 'Test': None}


    # def _LaunchTrainROM(self, mu_train, mu_init_train):
    def _LaunchTrainROM(self, mu_train, start_from_closest_mu=False):
        """
        This method should be parallel capable
        """
        self._LaunchFOM(mu_train, start_from_closest_mu = start_from_closest_mu)
        self._LaunchComputeSolutionBasis(mu_train)



    # def _LaunchFOM(self, mu_train, mu_train_init, gid_and_vtk_name='FOM_Fit'):
    def _LaunchFOM(self, mu_train, start_from_closest_mu=False, gid_and_vtk_name='FOM_Fit'):
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        NonConvergedSolutionsGathering = self.general_rom_manager_parameters["store_nonconverged_fom_solutions"].GetBool()
        RemoveFixedDoFs = self.general_rom_manager_parameters["ROM"]["remove_fixed_dofs"].GetBool()
        failed_mu_list = []
        for Id, mu in enumerate(mu_train):
            try:
                print(f'Starting FOM simulation for mu: {mu}')
                fom_in_database, _ = self.data_base.check_if_in_database("FOM", mu)
                nonconverged_fom_in_database, _ = self.data_base.check_if_in_database("NonconvergedFOM", mu)
                if not fom_in_database or (NonConvergedSolutionsGathering and not nonconverged_fom_in_database):
                    parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
                    parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy) #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
                    parameters_copy = self._StoreResultsByName(parameters_copy,gid_and_vtk_name,mu,Id)
                    materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
                    self.UpdateMaterialParametersFile(materials_file_name, mu)
                    # if mu_train_init[0] is not None:
                    if start_from_closest_mu is True:
                        # print('Init mu: ', mu_train_init[Id])
                        mu_init = self.data_base.get_closest_mu_in_FOM(mu)
                        if mu_init is not None:
                            print(f'Starting from initial mu: {mu_init}')
                            parameters_copy = self._AddInitializeFromSnapshotProcessToProjectParameters(parameters_copy, mu_init)
                    model = KratosMultiphysics.Model()
                    analysis_stage_class = self._GetAnalysisStageClass(parameters_copy)
                    simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)
                    if NonConvergedSolutionsGathering:
                        simulation = self.ActivateNonconvergedSolutionsGathering(simulation)
                    simulation.Run()
                    if not fom_in_database:
                        self.data_base.add_to_database("QoI_FOM", mu, simulation.GetFinalData())
                        for process in simulation._GetListOfOutputProcesses():
                            if isinstance(process, CalculateRomBasisOutputProcess):
                                BasisOutputProcess = process
                        SnapshotsMatrix = BasisOutputProcess._GetSnapshotsMatrix() #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
                        self.data_base.add_to_database("FOM", mu, SnapshotsMatrix)
                        if RemoveFixedDoFs:
                            FixityConfig = BasisOutputProcess._GetFixityConfiguration()
                            self.data_base.add_to_database("Fixity", mu, FixityConfig)
                if NonConvergedSolutionsGathering:
                        self.data_base.add_to_database("NonconvergedFOM", mu, simulation.GetNonconvergedSolutions())
            except:
                print(f'FOM simulation for mu: {mu} failed. Moving on to next one')
                failed_mu_list.append(mu)

        self.GenerateDatabaseSummary()
        
        if len(failed_mu_list) > 0:
            err_msg = f'Failed to simulate mu values: {failed_mu_list}'
            raise Exception(err_msg)


    def _LaunchComputeSolutionBasis(self, mu_train):
        RemoveFixedDoFs = self.general_rom_manager_parameters["ROM"]["remove_fixed_dofs"].GetBool()
        if RemoveFixedDoFs:
            # Check that all fixity configs are equal and retrieve it
            fixity_config_list = self.data_base.check_constant_fixity_and_get_list(mu_train)
        else:
            fixity_config_list = []
        in_database, hash_basis = self.data_base.check_if_in_database("RightBasis", mu_train)
        if not in_database:
            BasisOutputProcess = self.InitializeDummySimulationForBasisOutputProcess()
            if self.general_rom_manager_parameters["ROM"]["use_non_converged_sols"].GetBool():
                snapshots_matrix = self.data_base.get_snapshots_matrix_from_database(mu_train, table_name='NonconvergedFOM')
            else:
                snapshots_matrix = self.data_base.get_snapshots_matrix_from_database(mu_train, table_name='FOM')
            if RemoveFixedDoFs:
                # Apply remove fixed components from snapshots_matrix before computing its SVD
                BasisOutputProcess._ApplyFixityOnSnapshotsMatrix(fixity_config_list,snapshots_matrix)
            u,sigma = BasisOutputProcess._ComputeSVD(snapshots_matrix) #TODO the case for non-converged might be too large for single opeartion, add partitioned svd
            BasisOutputProcess._PrintRomBasis(u, sigma, fixity_config_list) #Calling the RomOutput Process for creating the RomParameter.json
            self.data_base.add_to_database("RightBasis", mu_train, u )
            self.data_base.add_to_database("SingularValues_Solution", mu_train, sigma )
        else:
            BasisOutputProcess = self.InitializeDummySimulationForBasisOutputProcess()
            _ , hash_sigma = self.data_base.check_if_in_database("SingularValues_Solution", mu_train)
            BasisOutputProcess._PrintRomBasis(self.data_base.get_single_numpy_from_database(hash_basis),
                                              self.data_base.get_single_numpy_from_database(hash_sigma),
                                              fixity_config_list) #this updates the RomParameters.json
        self.GenerateDatabaseSummary()

    def _LoadSolutionBasis(self, mu_train):
        if self.general_rom_manager_parameters["ROM"]["remove_fixed_dofs"].GetBool():
            fixity_config_list = self.data_base.check_constant_fixity_and_get_list(mu_train)
        else:
            fixity_config_list = []
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
            BasisOutputProcess = self.InitializeDummySimulationForBasisOutputProcess()
            _ , hash_sigma = self.data_base.check_if_in_database("SingularValues_Solution", mu_train)
            BasisOutputProcess._PrintRomBasis(self.data_base.get_single_numpy_from_database(hash_basis),
                                              self.data_base.get_single_numpy_from_database(hash_sigma),
                                              fixity_config_list) #this updates the RomParameters.json
        self.GenerateDatabaseSummary()



    def _LaunchROM(self, mu_run, mu_train = None, gid_and_vtk_name='ROM_Fit', nn_rom_interface = None):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        BasisOutputProcess = None
        nn_model_name = None
        if nn_rom_interface is not None:
            if nn_rom_interface.get_using_manual_model_path():
                mu_train = None
                nn_model_name = nn_rom_interface.get_nn_model_name()
        for Id, mu in enumerate(mu_run):
            print('Launching ROM for ', mu)
            in_database, _ = self.data_base.check_if_in_database("ROM", mu, mu_aux = mu_train, nn_model_name = nn_model_name)
            if not in_database:
                parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
                parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)  #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
                parameters_copy = self._StoreResultsByName(parameters_copy,gid_and_vtk_name,mu,Id)
                materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
                self.UpdateMaterialParametersFile(materials_file_name, mu)
                model = KratosMultiphysics.Model()
                analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy, nn_rom_interface=nn_rom_interface))
                simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)

                simulation.Run()
                self.data_base.add_to_database("QoI_ROM", mu, simulation.GetFinalData(), mu_aux = mu_train, nn_model_name = nn_model_name)
                for process in simulation._GetListOfOutputProcesses():
                    if isinstance(process, CalculateRomBasisOutputProcess):
                        BasisOutputProcess = process
                SnapshotsMatrix = BasisOutputProcess._GetSnapshotsMatrix() #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
                self.data_base.add_to_database("ROM", mu, SnapshotsMatrix, mu_aux = mu_train, nn_model_name = nn_model_name)

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
                    parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)
                    parameters_copy = self._StoreNoResults(parameters_copy)
                    materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
                    self.UpdateMaterialParametersFile(materials_file_name, mu)
                    model = KratosMultiphysics.Model()
                    analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
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





    def _LaunchTrainHROM(self, mu_train):
        """
        This method should be parallel capable
        """
        in_database_elems, hash_z =  self.data_base.check_if_in_database("HROM_Elements", mu_train)
        in_database_weights, hash_w =  self.data_base.check_if_in_database("HROM_Weights", mu_train)
        if not in_database_elems and not in_database_weights:
            with open(self.project_parameters_name,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            simulation = None
            for mu in mu_train:
                in_database, _ = self.data_base.check_if_in_database("ResidualsProjected", mu)
                if not in_database:
                    parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
                    parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)
                    parameters_copy = self._StoreNoResults(parameters_copy)
                    materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
                    self.UpdateMaterialParametersFile(materials_file_name, mu)
                    model = KratosMultiphysics.Model()
                    analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
                    simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)
                    simulation.Run()
                    ResidualProjected = simulation.GetHROM_utility()._GetResidualsProjectedMatrix() #TODO flush intermediately the residuals projected to cope with large models.
                    self.data_base.add_to_database("ResidualsProjected", mu, ResidualProjected )
            RedidualsSnapshotsMatrix = self.data_base.get_snapshots_matrix_from_database(mu_train, table_name="ResidualsProjected")
            u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(RedidualsSnapshotsMatrix,
            self.general_rom_manager_parameters["HROM"]["element_selection_svd_truncation_tolerance"].GetDouble()) #TODO load basis for residuals projected. Can we truncate it only, not compute the whole SVD but only return the respective number of singular vectors?
            if simulation is None:
                HROM_utility = self.InitializeDummySimulationForHromTrainingUtility()
            else:
                HROM_utility = simulation.GetHROM_utility()
            HROM_utility.hyper_reduction_element_selector.SetUp(u, InitialCandidatesSet = HROM_utility.candidate_ids)
            HROM_utility.hyper_reduction_element_selector.Run()
            if not HROM_utility.hyper_reduction_element_selector.success:
                KratosMultiphysics.Logger.PrintWarning("HRomTrainingUtility", "The Empirical Cubature Method did not converge using the initial set of candidates. Launching again without initial candidates.")
                #Imposing an initial candidate set can lead to no convergence. Restart without imposing the initial candidate set
                HROM_utility.hyper_reduction_element_selector.SetUp(u, InitialCandidatesSet = None)
                HROM_utility.hyper_reduction_element_selector.Run()
            z = HROM_utility.hyper_reduction_element_selector.z
            w = HROM_utility.hyper_reduction_element_selector.w
            self.data_base.add_to_database("HROM_Elements", mu_train, np.squeeze(z))
            self.data_base.add_to_database("HROM_Weights", mu_train, np.squeeze(w))
            HROM_utility.AppendHRomWeightsToRomParameters()
            HROM_utility.CreateHRomModelParts()
        # elif (in_database_elems and not in_database_weights) or (not in_database_elems and in_database_weights):
        #     error #if one of the two is there but not the other, the database is corrupted
        else:
            HROM_utility = self.InitializeDummySimulationForHromTrainingUtility()
            #doing this ensures the elements and weights npy files contained in the rom_data folder are the ones to use.i.e. they are not from old runs of the rom manager
            HROM_utility.hyper_reduction_element_selector.w = self.data_base.get_single_numpy_from_database(hash_w)
            HROM_utility.hyper_reduction_element_selector.z = self.data_base.get_single_numpy_from_database(hash_z)
            HROM_utility.AppendHRomWeightsToRomParameters()
            HROM_utility.CreateHRomModelParts()
        self.GenerateDatabaseSummary()

    def _LaunchHROM(self, mu_train, gid_and_vtk_name ='HROM_Fit'):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        BasisOutputProcess = None
        for Id, mu in enumerate(mu_train):
            in_database, _ = self.data_base.check_if_in_database("HROM", mu)
            if not in_database:
                parameters_copy = self.UpdateProjectParameters(parameters.Clone(), mu)
                parameters_copy = self._AddBasisCreationToProjectParameters(parameters_copy)
                parameters_copy = self._StoreResultsByName(parameters_copy,gid_and_vtk_name,mu,Id)
                materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
                self.UpdateMaterialParametersFile(materials_file_name, mu)
                model = KratosMultiphysics.Model()
                analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy))
                simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)
                simulation.Run()
                self.data_base.add_to_database("QoI_HROM", mu, simulation.GetFinalData())
                for process in simulation._GetListOfOutputProcesses():
                    if isinstance(process, CalculateRomBasisOutputProcess):
                        BasisOutputProcess = process
                SnapshotsMatrix = BasisOutputProcess._GetSnapshotsMatrix() #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
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
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
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
            materials_file_name = parameters_copy["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
            self.UpdateMaterialParametersFile(materials_file_name, mu)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters_copy, nn_rom_interface=nn_rom_interface))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)
            simulation.Run()
            self.QoI_Run_ROM.append(simulation.GetFinalData())


    def _LaunchRunHROM(self, mu_run, use_full_model_part):
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
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters_copy, mu)
            simulation.Run()
            self.QoI_Run_HROM.append(simulation.GetFinalData())

    def _LaunchTrainNeuralNetwork(self, mu_train, mu_validation):
        RomNeuralNetworkTrainer = self._TryImportNNTrainer()
        rom_nn_trainer = RomNeuralNetworkTrainer(self.general_rom_manager_parameters, mu_train, mu_validation, self.data_base)
        rom_nn_trainer.TrainNetwork()
        self.data_base.add_to_database("Neural_Network", mu_train , None)
        rom_nn_trainer.EvaluateNetwork()


    def _LaunchTestNeuralNetworkReconstruction(self,mu_train, mu_validation, mu_test, custom_model_path = None):
        RomNeuralNetworkTrainer = self._TryImportNNTrainer()
        rom_nn_trainer = RomNeuralNetworkTrainer(self.general_rom_manager_parameters, mu_train, mu_validation, self.data_base)
        rom_nn_trainer.EvaluateNetwork(mu_eval = mu_test, custom_model_path = custom_model_path)

    def InitializeDummySimulationForBasisOutputProcess(self):
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters = self._AddBasisCreationToProjectParameters(parameters)
        parameters = self._StoreNoResults(parameters)
        model = KratosMultiphysics.Model()
        analysis_stage_class = self._GetAnalysisStageClass(parameters)
        simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
        simulation.Initialize()
        for process in simulation._GetListOfOutputProcesses():
            if isinstance(process, CalculateRomBasisOutputProcess):
                BasisOutputProcess = process
        return BasisOutputProcess


    def InitializeDummySimulationForHromTrainingUtility(self):
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters = self._AddBasisCreationToProjectParameters(parameters)
        parameters = self._StoreNoResults(parameters)
        model = KratosMultiphysics.Model()
        analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
        simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
        simulation.Initialize()
        return simulation.GetHROM_utility()


    def InitializeDummySimulationForPetrovGalerkinTrainingUtility(self):
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters = self._AddBasisCreationToProjectParameters(parameters)
        parameters = self._StoreNoResults(parameters)
        model = KratosMultiphysics.Model()
        analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
        simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
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
            elif simulation_to_run=='lspg_ANN':
                f['train_hrom']=False
                f['run_hrom']=False
                f['projection_strategy']="lspg_ann"
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

    def _AddBasisCreationToProjectParameters(self, parameters):
        #FIXME make sure no other rom_output already existed. If so, erase the prior and keep only the one in self.general_rom_manager_parameters["ROM"]
        parameters["output_processes"].AddEmptyArray("rom_output")
        rom_basis_parameters = self._SetUpRomBasisParameters()
        parameters["output_processes"]["rom_output"].Append(rom_basis_parameters)

        return parameters
    

    def _AddInitializeFromSnapshotProcessToProjectParameters(self, parameters, mu_init):
        if len(mu_init) > 0:
            init_snapshot = self.data_base.get_snapshots_matrix_from_database([mu_init], table_name='FOM').reshape(-1)
            if init_snapshot is None:
                raise Exception(f"No snapshot found for initial mu: {mu_init}.")
        
            parameters["processes"].AddEmptyArray("rom_process")

            process_parameters = KratosMultiphysics.Parameters("""{
                "python_module" : "initialize_from_snapshot_process",
                "kratos_module" : "KratosMultiphysics.RomApplication",
                "process_name"  : "InitializeFromSnapshotProcess",
                "help"          : "This process should initialize the simulation from a custom snapshot",
                "Parameters"    : {}
                }""")

            process_parameters["Parameters"] = InitializeFromSnapshotProcess.GetDefaultParameters()
            process_parameters["Parameters"]["nodal_unknowns"] = self.general_rom_manager_parameters["ROM"]["nodal_unknowns"]
            process_parameters["Parameters"]["model_part_name"] = self.general_rom_manager_parameters["ROM"]["model_part_name"]
            # process_parameters["Parameters"]["move_mesh"] = parameters["solver_settings"]["move_mesh"]
            process_parameters["Parameters"]["initial_snapshot"].SetVector(init_snapshot.tolist())

            parameters["processes"]["rom_process"].Append(process_parameters)

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
            "ROM":{
                "svd_truncation_tolerance": 1e-5,
                "model_part_name": "Structure",                            // This changes depending on the simulation: Structure, FluidModelPart, ThermalPart #TODO: Idenfity it automatically
                "nodal_unknowns": ["DISPLACEMENT_X","DISPLACEMENT_Y"],     // Main unknowns. Snapshots are taken from these
                "rom_basis_output_format": "numpy",
                "rom_basis_output_name": "RomParameters",
                "rom_basis_output_folder": "rom_data",
                "snapshots_control_type": "step",                          // "step", "time"
                "snapshots_interval": 1,
                "snapshots_control_is_periodic": false,
                "print_singular_values": false,
                "remove_fixed_dofs": false,
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
                        "model_number": 0,   // out of the models existing for the same parameters, this is the model that will be lauched
                        "custom_model_path": ""
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

    def _SetUpRomBasisParameters(self):
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
            "print_singular_values",
            "remove_fixed_dofs",
            "snapshots_control_is_periodic"
        ]

        for key in keys_to_copy:
            if key in rom_params.keys():
                defaults["Parameters"][key] = rom_params[key]

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
                    "snapshots_control_is_periodic": false,
                    "nodal_unknowns":  [],
                    "rom_basis_output_format": "json",
                    "rom_basis_output_name": "RomParameters",
                    "rom_basis_output_folder": "rom_data",
                    "svd_truncation_tolerance": 1e-3,
                    "print_singular_values": false,
                    "remove_fixed_dofs": false
                }
            }""")
        return rom_training_parameters


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
        

    def GenerateOrderedFOMSnapshotsMatrix(self, mu_list):
        for Id, mu in enumerate(mu_list):
            in_database, _ = self.data_base.check_if_in_database("FOM", mu)
            if not in_database:
                err_msg = f'Missing FOM simulation for mu: {mu}'
                raise Exception(err_msg)
            
        snapshots_matrix = self.data_base.get_ordered_snapshots_matrix_from_database(mu_list, table_name='FOM')
        # snapshots_matrix = self.data_base.get_snapshots_matrix_from_database(mu_list, table_name='FOM')
        return snapshots_matrix
        