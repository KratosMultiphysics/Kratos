import KratosMultiphysics
from KratosMultiphysics.RomApplication.rom_testing_utilities import SetUpSimulationInstance
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
import numpy as np
import importlib
import json




class RomManager(object):

    def __init__(self,project_parameters_name, general_rom_manager_parameters, CustomizeSimulation, UpdateProjectParameters):
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



    def Fit(self, mu_train=None):
        if mu_train is None:
            mu_train = ['single case with parameters already contained in the ProjectParameters.json and CustomSimulation']
        #######################
        ######  Galerkin ######
        if self.general_rom_manager_parameters["projection_strategy"].GetString() == "galerkin":
            training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
            if any(item == "ROM" for item in training_stages):
                fom_snapshots = self.LaunchTrainROM(mu_train)
                self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
                rom_snapshots = self.LaunchROM(mu_train)
                self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)

            if any(item == "HROM" for item in training_stages):
                #FIXME there will be an error if we only train HROM, but not ROM
                self._ChangeRomFlags(simulation_to_run = "trainHROMGalerkin")
                self.LaunchTrainHROM(mu_train)
                self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
                hrom_snapshots = self.LaunchHROM(mu_train)
                self.ROMvsHROM_train = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)
        #######################


        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif self.general_rom_manager_parameters["projection_strategy"].GetString() == "lspg":
            training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
            if any(item == "ROM" for item in training_stages):
                fom_snapshots = self.LaunchTrainROM(mu_train)
                self._ChangeRomFlags(simulation_to_run = "lspg")
                rom_snapshots = self.LaunchROM(mu_train)
                self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)
            if any(item == "HROM" for item in training_stages):
                raise Exception('Sorry, Hyper Reduction not yet implemented for lspg')
        #######################################


        ##########################
        ###  Petrov Galerkin   ###
        elif self.general_rom_manager_parameters["projection_strategy"].GetString() == "petrov_galerkin":
            training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
            if any(item == "ROM" for item in training_stages):
                fom_snapshots = self.LaunchTrainROM(mu_train)
                self._ChangeRomFlags(simulation_to_run = "TrainPG")
                self.TrainPG(mu_train)
                self._ChangeRomFlags(simulation_to_run = "PG")
                rom_snapshots = self.LaunchROM(mu_train)
                self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)
            if any(item == "HROM" for item in training_stages):
                #FIXME there will be an error if we only train HROM, but not ROM
                self._ChangeRomFlags(simulation_to_run = "trainHROMPetrovGalerkin")
                self.LaunchTrainHROM(mu_train)
                self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
                hrom_snapshots = self.LaunchHROM(mu_train)
                self.ROMvsHROM_train = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)
        ##########################



    def Test(self, mu_test=None):
        if mu_test is None:
            mu_test = ['single case with parameters already contained in the ProjectParameters.json and CustomSimulation']
        #######################
        ######  Galerkin ######
        if self.general_rom_manager_parameters["projection_strategy"].GetString() == "galerkin":
            testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
            if any(item == "ROM" for item in testing_stages):
                fom_snapshots = self.LaunchTestFOM(mu_test)
                self._ChangeRomFlags(simulation_to_run = "GalerkinROM")
                rom_snapshots = self.LaunchTestROM(mu_test)
                self.ROMvsFOM_test = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)

            if any(item == "HROM" for item in testing_stages):
                #FIXME there will be an error if we only test HROM, but not ROM
                self._ChangeRomFlags(simulation_to_run = "runHROMGalerkin")
                hrom_snapshots = self.LaunchTestHROM(mu_test)
                self.ROMvsHROM_test = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)


        #######################################
        ##  Least-Squares Petrov Galerkin   ###
        elif self.general_rom_manager_parameters["projection_strategy"].GetString() == "lspg":
            testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
            if any(item == "ROM" for item in testing_stages):
                fom_snapshots = self.LaunchTestFOM(mu_test)
                self._ChangeRomFlags(simulation_to_run = "lspg")
                rom_snapshots = self.LaunchTestROM(mu_test)
                self.ROMvsFOM_test = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)
            if any(item == "HROM" for item in testing_stages):
                raise Exception('Sorry, Hyper Reduction not yet implemented for lspg')
        #######################################


        ##########################
        ###  Petrov Galerkin   ###
        elif self.general_rom_manager_parameters["projection_strategy"].GetString() == "petrov_galerkin":
            testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
            if any(item == "ROM" for item in testing_stages):
                fom_snapshots = self.LaunchTestFOM(mu_test)
                self._ChangeRomFlags(simulation_to_run = "PG")
                rom_snapshots = self.LaunchTestROM(mu_test)
                self.ROMvsFOM_test = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)
            if any(item == "HROM" for item in testing_stages):
                #FIXME there will be an error if we only train HROM, but not ROM
                self._ChangeRomFlags(simulation_to_run = "runHROMPetrovGalerkin")
                hrom_snapshots = self.LaunchTestHROM(mu_test)
                self.ROMvsHROM_test = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)
        ##########################




    def PrintErrors(self):
        training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
        testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
        if any(item == "ROM" for item in training_stages):
            print("approximation error in train set FOM vs ROM: ", self.ROMvsFOM_train)
        if any(item == "HROM" for item in training_stages) and not(self.general_rom_manager_parameters["projection_strategy"].GetString() == "lspg"):
            print("approximation error in train set ROM vs HROM: ", self.ROMvsHROM_train)
        if any(item == "ROM" for item in testing_stages):
            print("approximation error in test set FOM vs ROM: ", self.ROMvsFOM_test)
        if any(item == "HROM" for item in testing_stages) and not(self.general_rom_manager_parameters["projection_strategy"].GetString() == "lspg"):
            print("approximation error in test set ROM vs HROM: ",  self.ROMvsHROM_test)


    def LaunchTrainROM(self, mu_train):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_train):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._AddBasisCreationToProjectParameters(parameters) #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters = self._StoreResultsByName(parameters,'FOM',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)
        BasisOutputProcess._PrintRomBasis(SnapshotsMatrix) #Calling the RomOutput Process for creating the RomParameter.json

        return SnapshotsMatrix


    def LaunchROM(self, mu_train):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_train):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._AddBasisCreationToProjectParameters(parameters)  #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters = self._StoreResultsByName(parameters,'ROM',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix


    def TrainPG(self, mu_train):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        PetrovGalerkinTrainMatrix = []
        for Id, mu in enumerate(mu_train):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._StoreResultsByName(parameters,'ROM',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            PetrovGalerkinTrainMatrix.append(simulation.GetPetrovGalerkinTrainUtility()._GetSnapshotsMatrix()) #TODO is the best way of extracting the Projected Residuals calling the HROM residuals utility?
        simulation.GetPetrovGalerkinTrainUtility().CalculateAndSaveBasis(np.block(PetrovGalerkinTrainMatrix))


    def LaunchTrainHROM(self, mu_train):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        RedidualsSnapshotsMatrix = []
        for mu in mu_train:
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._StoreNoResults(parameters)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            RedidualsSnapshotsMatrix.append(simulation.GetHROM_utility()._GetResidualsProjectedMatrix()) #TODO is the best way of extracting the Projected Residuals calling the HROM residuals utility?
        RedidualsSnapshotsMatrix = np.block(RedidualsSnapshotsMatrix)
        u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(RedidualsSnapshotsMatrix,
        self.hrom_training_parameters["element_selection_svd_truncation_tolerance"].GetDouble())
        simulation.GetHROM_utility().hyper_reduction_element_selector.SetUp(u)
        simulation.GetHROM_utility().hyper_reduction_element_selector.Run()
        simulation.GetHROM_utility().AppendHRomWeightsToRomParameters()


    def LaunchHROM(self, mu_train):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_train):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._AddBasisCreationToProjectParameters(parameters)
            parameters = self._StoreResultsByName(parameters,'HROM',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix


    def LaunchTestFOM(self, mu_test):
        """
        This method should be parallel capable
        """
        # FIXME We must include a method to retrive solutions in the nodes and stop using the CalculateRomBasisOutputProcess to stote snapshots
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_test):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._AddBasisCreationToProjectParameters(parameters) #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters = self._StoreResultsByName(parameters,'FOM_Test',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)


        return SnapshotsMatrix


    def LaunchTestROM(self, mu_test):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_test):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._AddBasisCreationToProjectParameters(parameters)  #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters = self._StoreResultsByName(parameters,'ROM_Test',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix



    def LaunchTestHROM(self, mu_test):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for Id, mu in enumerate(mu_test):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._AddBasisCreationToProjectParameters(parameters)
            parameters = self._StoreResultsByName(parameters,'HROM_Test',mu,Id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if isinstance(process, CalculateRomBasisOutputProcess):
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix



    def _ChangeRomFlags(self, simulation_to_run = 'ROM'):
        """
        This method updates the Flags present in the RomParameters.json file
        for launching the correct part of the ROM workflow
        """
        #other options: "trainHROM", "runHROM"
        #taken from code by Philipa & Catharina
        parameters_file_name = './RomParameters.json'
        with open(parameters_file_name, 'r+') as parameter_file:
            f=json.load(parameter_file)
            if simulation_to_run=='GalerkinROM':
                f['projection_strategy']="galerkin"
                f['train_hrom']=False
                f['run_hrom']=False
            elif simulation_to_run=='trainHROMGalerkin':
                f['train_hrom']=True
                f['run_hrom']=False
            elif simulation_to_run=='runHROMGalerkin':
                f['projection_strategy']="galerkin"
                f['train_hrom']=False
                f['run_hrom']=True
            elif simulation_to_run=='lspg':
                f['train_hrom']=False
                f['run_hrom']=False
                f['projection_strategy']="lspg"
            elif simulation_to_run=='TrainPG':
                f['train_hrom']=False
                f['run_hrom']=False
                f['projection_strategy']="lspg"
                f['train_petrov_galerkin']['train'] = True
                f['train_petrov_galerkin']['basis_strategy'] = self.general_rom_manager_parameters["ROM"]["petrov_galerkin_training_parameters"]["basis_strategy"].GetString()
                f['train_petrov_galerkin']['include_phi'] = self.general_rom_manager_parameters["ROM"]["petrov_galerkin_training_parameters"]["include_phi"].GetBool()
                f['train_petrov_galerkin']['svd_truncation_tolerance'] = self.general_rom_manager_parameters["ROM"]["petrov_galerkin_training_parameters"]["svd_truncation_tolerance"].GetDouble()
            elif simulation_to_run=='PG':
                f['train_hrom']=False
                f['run_hrom']=False
                f['projection_strategy']="petrov_galerkin"
                f['train_petrov_galerkin']['train'] = False
            elif simulation_to_run=='trainHROMPetrovGalerkin':
                f['train_hrom']=True
                f['run_hrom']=False
                f['projection_strategy']="petrov_galerkin"
                f['train_petrov_galerkin']['train'] = False
            elif simulation_to_run=='runHROMPetrovGalerkin':
                f['train_hrom']=False
                f['run_hrom']=True
                f['projection_strategy']="petrov_galerkin"
                f['train_petrov_galerkin']['train'] = False
            else:
                raise Exception(f'Unknown flag "{simulation_to_run}" change for RomParameters.json')
            parameter_file.seek(0)
            json.dump(f,parameter_file,indent=4)
            parameter_file.truncate()





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
                    "svd_truncation_tolerance": 1e-3
                }
            }""")
        return rom_training_parameters



    def _GetDefaulHromTrainingParameters(self):
        hrom_training_parameters = KratosMultiphysics.Parameters("""{
                "element_selection_type": "empirical_cubature",
                "element_selection_svd_truncation_tolerance": 0,
                "echo_level" : 0,
                "create_hrom_visualization_model_part" : true
            }""")
        return hrom_training_parameters












