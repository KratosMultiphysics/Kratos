import KratosMultiphysics
from KratosMultiphysics.RomApplication.rom_testing_utilities import SetUpSimulationInstance
from KratosMultiphysics.RomApplication.calculate_rom_basis_output_process import CalculateRomBasisOutputProcess
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
import numpy as np
import importlib
import json




class RomManager(object):

    def __init__(self,mu, project_parameters_name, general_rom_manager_parameters, rom_training_parameters, hrom_training_parameters, CustomizeSimulation, UpdateProjectParameters):
        #FIXME There is room for improvement of this prototype:
        # - Use a method (upcoming) for smothly retrieving solutions. In here we are using the RomBasisOutput process in order to store the solutions
        # - There is some redundancy between the methods that launch the simulations. Can we create a single method?
        # - We must use the Hyper-Reduced model part for the HROM simulations. (The one in which we have eliminated the non selected elements and conditions)
        # - Not yet paralellised with COMPSs
        self.mu = mu
        self.project_parameters_name = project_parameters_name
        self.rom_training_parameters = rom_training_parameters
        self.hrom_training_parameters = hrom_training_parameters
        self.general_rom_manager_parameters = general_rom_manager_parameters
        self.CustomizeSimulation = CustomizeSimulation
        self.UpdateProjectParameters = UpdateProjectParameters
        if len(mu)>1:
            self.rom_training_parameters["Parameters"]["multiple_simulations"].SetBool(True)


    def Fit(self):
        training_stages = self.general_rom_manager_parameters["rom_stages_to_train"].GetStringArray()
        if any(item == "ROM" for item in training_stages):
            fom_snapshots = self.LauchTraiROM()
            rom_snapshots = self.LauchROM()
            self.ROMvsFOM_train = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)

        if any(item == "HROM" for item in training_stages):
            #FIXME there will be an error if we only train HROM, but not ROM
            self.LauchTraiHROM()
            hrom_snapshots = self.LauchHROM()
            self.ROMvsHROM_train = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)



    def Test(self, mu_test):
        testing_stages = self.general_rom_manager_parameters["rom_stages_to_test"].GetStringArray()
        if any(item == "ROM" for item in testing_stages):
            fom_snapshots = self.LauchTestFOM(mu_test)
            rom_snapshots = self.LauchTestROM(mu_test)
            self.ROMvsFOM_test = np.linalg.norm(fom_snapshots - rom_snapshots)/ np.linalg.norm(fom_snapshots)

        if any(item == "HROM" for item in testing_stages):
            #FIXME there will be an error if we only test HROM, but not ROM
            hrom_snapshots = self.LauchTestHROM(mu_test)
            self.ROMvsHROM_test = np.linalg.norm(rom_snapshots - hrom_snapshots) / np.linalg.norm(rom_snapshots)


    def PrintErrors(self):
        print("approximation error in train set FOM vs ROM: ", self.ROMvsFOM_train)
        print("approximation error in train set ROM vs HROM: ", self.ROMvsHROM_train)
        print("approximation error in test set FOM vs ROM: ", self.ROMvsFOM_test)
        print("approximation error in test set ROM vs HROM: ",  self.ROMvsHROM_test)






    def LauchTraiROM(self):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        SnapshotsMatrix = []
        for id, mu in enumerate(self.mu):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._AddBasisCreationToProjectParameters(parameters) #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters = self._StoreResultsByName(parameters,'FOM',mu,id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if type(process) == CalculateRomBasisOutputProcess:
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)
        BasisOutputProcess._PrintRomBasis(SnapshotsMatrix) #Calling the RomOutput Process for creating the RomParameter.json

        return SnapshotsMatrix


    def LauchROM(self):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for id, mu in enumerate(self.mu):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._AddBasisCreationToProjectParameters(parameters)  #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters = self._StoreResultsByName(parameters,'ROM',mu,id)
            self._ChangeRomFlags(simulation_to_run = "ROM")
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if type(process) == CalculateRomBasisOutputProcess:
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix


    def LauchTraiHROM(self):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        RedidualsSnapshotsMatrix = []
        for id, mu in enumerate(self.mu):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._StoreNoResults(parameters)
            self._ChangeRomFlags(simulation_to_run = "trainHROM")
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


    def LauchHROM(self):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for id, mu in enumerate(self.mu):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._AddBasisCreationToProjectParameters(parameters)
            parameters = self._StoreResultsByName(parameters,'HROM',mu,id)
            self._ChangeRomFlags(simulation_to_run = "runHROM")
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if type(process) == CalculateRomBasisOutputProcess:
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix


    def LauchTestFOM(self, mu_test):
        """
        This method should be parallel capable
        """
        # FIXME Since we are using the RomBasisOutputProcess for storing the snapshots in this prototype, calling it with a single training parameter will automatically generate
        # the RomBasis. We must include a method to retrive solutions in the nodes
        # In this example it will not cause problems, because since there are more than 1 parameter "mu", the Finalize() method wont re-create the RomParameters.json
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        SnapshotsMatrix = []
        for id, mu in enumerate(mu_test):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._AddBasisCreationToProjectParameters(parameters) #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters = self._StoreResultsByName(parameters,'FOM_Test',mu,id)
            model = KratosMultiphysics.Model()
            analysis_stage_class = self._GetAnalysisStageClass(parameters)
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if type(process) == CalculateRomBasisOutputProcess:
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)


        return SnapshotsMatrix


    def LauchTestROM(self, mu_test):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for id, mu in enumerate(mu_test):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._AddBasisCreationToProjectParameters(parameters)  #TODO stop using the RomBasisOutputProcess to store the snapshots. Use instead the upcoming build-in function
            parameters = self._StoreResultsByName(parameters,'ROM_Test',mu,id)
            self._ChangeRomFlags(simulation_to_run = "ROM")
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if type(process) == CalculateRomBasisOutputProcess:
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix



    def LauchTestHROM(self, mu_test):
        """
        This method should be parallel capable
        """
        with open(self.project_parameters_name,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        SnapshotsMatrix = []
        for id, mu in enumerate(mu_test):
            parameters = self.UpdateProjectParameters(parameters, mu)
            parameters = self._AddBasisCreationToProjectParameters(parameters)
            parameters = self._StoreResultsByName(parameters,'HROM_Test',mu,id)
            self._ChangeRomFlags(simulation_to_run = "runHROM")
            model = KratosMultiphysics.Model()
            analysis_stage_class = type(SetUpSimulationInstance(model, parameters))
            simulation = self.CustomizeSimulation(analysis_stage_class,model,parameters)
            simulation.Run()
            for process in simulation._GetListOfOutputProcesses():
                if type(process) == CalculateRomBasisOutputProcess:
                    BasisOutputProcess = process
            SnapshotsMatrix.append(BasisOutputProcess._GetSnapshotsMatrix()) #TODO add a CustomMethod() as a standard method in the Analysis Stage to retrive some solution
        SnapshotsMatrix = np.block(SnapshotsMatrix)

        return SnapshotsMatrix



    def _ChangeRomFlags(self, simulation_to_run = 'ROM'):
        """
        This method updates the Flags present in the RomParameters.json file
        for launching the correct part of the ROM workflow: ROM, TrainHROM, RunHROM
        """
        #other options: "trainHROM", "runHROM"
        #taken from code by Philipa & Catharina
        parameters_file_name = './RomParameters.json'
        with open(parameters_file_name, 'r+') as parameter_file:
            f=json.load(parameter_file)
            if simulation_to_run=='ROM':
                f['train_hrom']=False
                f['run_hrom']=False
            elif simulation_to_run=='trainHROM':
                f['train_hrom']=True
                f['run_hrom']=False
            elif simulation_to_run=='runHROM':
                f['train_hrom']=False
                f['run_hrom']=True
            else:
                print('Unknown operation. Add new rule!')
            parameter_file.seek(0)
            json.dump(f,parameter_file,indent=4)
            parameter_file.truncate()





    def _AddBasisCreationToProjectParameters(self, parameters):
        parameters["output_processes"].AddEmptyArray("rom_output")
        parameters["output_processes"]["rom_output"].Append(self.rom_training_parameters)

        return parameters




    def _StoreResultsByName(self,parameters,results_name,mu, id):

        if  self.general_rom_manager_parameters["output_name"].GetString() == "mu":
            case_name = (", ".join(map(str, mu)))
        elif self.general_rom_manager_parameters["output_name"].GetString() == "id":
            case_name = str(id)
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





    @staticmethod
    def _GetAnalysisStageClass(parameters):

        analysis_stage_module_name = parameters["analysis_stage"].GetString()
        analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
        analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

        analysis_stage_module = importlib.import_module(analysis_stage_module_name)
        analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

        return analysis_stage_class

















