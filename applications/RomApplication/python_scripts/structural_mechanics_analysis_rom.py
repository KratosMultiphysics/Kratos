import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom as solver_wrapper
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

import json
import numpy as np

class StructuralMechanicsAnalysisROM(StructuralMechanicsAnalysis):

    def __init__(self,model,project_parameters, hyper_reduction_element_selector = None, ECM_SVD_tolerance = 1e-6):
        super().__init__(model,project_parameters)
        if hyper_reduction_element_selector != None :
            if hyper_reduction_element_selector == "EmpiricalCubature":
                self.hyper_reduction_element_selector = EmpiricalCubatureMethod()
                self.ECM_SVD_tolerance = ECM_SVD_tolerance
                self.time_step_residual_matrix_container = []
            else:
                err_msg =  "The requested element selection method \"" + hyper_reduction_element_selector + "\" is not in the rom application\n"
                err_msg += "Available options are: \"EmpiricalCubature\""
                raise Exception(err_msg)
        else:
            self.hyper_reduction_element_selector = None

    #### Internal functions ####
    def _CreateSolver(self):
        """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """
        ## Solver construction
        with open('RomParameters.json') as rom_parameters:
            rom_settings = KratosMultiphysics.Parameters(rom_parameters.read())
            self.project_parameters["solver_settings"].AddValue("rom_settings", rom_settings["rom_settings"])
        return solver_wrapper.CreateSolverByParameters(self.model, self.project_parameters["solver_settings"],self.project_parameters["problem_data"]["parallel_type"].GetString())

    def _GetSimulationName(self):
        return "::[ROM Simulation]:: "

    def ModifyAfterSolverInitialize(self):
        """Here is where the ROM_BASIS is imposed to each node"""
        super().ModifyAfterSolverInitialize()
        computing_model_part = self._solver.GetComputingModelPart()
        with open('RomParameters.json') as f:
            data = json.load(f)
            nodal_dofs = len(data["rom_settings"]["nodal_unknowns"])
            nodal_modes = data["nodal_modes"]
            counter = 0
            rom_dofs= self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()
            for node in computing_model_part.Nodes:
                aux = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
                for j in range(nodal_dofs):
                    Counter=str(node.Id)
                    for i in range(rom_dofs):
                        aux[j,i] = nodal_modes[Counter][j][i]
                node.SetValue(romapp.ROM_BASIS, aux ) # ROM basis
                counter+=1
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                self.ResidualUtilityObject = romapp.RomResidualsUtility(self._GetSolver().GetComputingModelPart(), self.project_parameters["solver_settings"]["rom_settings"], self._GetSolver().get_solution_scheme())

    def FinalizeSolutionStep(self):
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                print('\n\n\n\nGenerating matrix of residuals')
                ResMat = self.ResidualUtilityObject.GetResiduals()
                NP_ResMat = np.array(ResMat, copy=False)
                self.time_step_residual_matrix_container.append(NP_ResMat)
        super().FinalizeSolutionStep()

    def Finalize(self):
        super().Finalize()
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                OriginalNumberOfElements = self._GetSolver().GetComputingModelPart().NumberOfElements()
                self. hyper_reduction_element_selector.SetUp(self._ObtainBasis(), OriginalNumberOfElements)
                self.hyper_reduction_element_selector.Run()
                self._CreateHyperReducedModelPart()

    def _ObtainBasis(self):
        ### Building the Snapshot matrix ####
        for i in range (len(self.time_step_residual_matrix_container)):
            if i == 0:
                SnapshotMatrix = self.time_step_residual_matrix_container[i]
            else:
                SnapshotMatrix = np.c_[SnapshotMatrix,self.time_step_residual_matrix_container[i]]
        ### Taking the SVD ###  (randomized and truncated)
        u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(SnapshotMatrix, self.ECM_SVD_tolerance)
        return u


    def _CreateHyperReducedModelPart(self):
        ModelPartName = self._GetSolver().settings["model_import_settings"]["input_filename"].GetString()
        current_model = KratosMultiphysics.Model()
        computing_model_part = current_model.CreateModelPart("main")
        model_part_io = KratosMultiphysics.ModelPartIO(ModelPartName)
        model_part_io.ReadModelPart(computing_model_part)
        hyper_reduced_model_part_help =   current_model.CreateModelPart("Helping")

        #TODO optimize implementation
        with open('ElementsAndWeights.json') as f:
            HR_data = json.load(f)
            for key in HR_data["Elements"].keys():
                for node in computing_model_part.GetElement(int(key)+1).GetNodes():
                    hyper_reduced_model_part_help.AddNode(node,0)
            for key in HR_data["Conditions"].keys():
                for node in computing_model_part.GetCondition(int(key)+1).GetNodes():
                    hyper_reduced_model_part_help.AddNode(node,0)

        # The HROM model part. It will include two sub-model parts. One for caculation, another one for visualization
        HROM_Model_Part =  current_model.CreateModelPart("HROM_Model_Part")

        # Building the COMPUTE_HROM submodel part
        hyper_reduced_model_part = HROM_Model_Part.CreateSubModelPart("COMPUTE_HROM")

        # TODO implement the hyper-reduced model part creation in C++
        with open('ElementsAndWeights.json') as f:
            HR_data = json.load(f)
            for originalSubmodelpart in computing_model_part.SubModelParts:
                hyperReducedSubmodelpart = hyper_reduced_model_part.CreateSubModelPart(originalSubmodelpart.Name)
                print(f'originalSubmodelpart.Name {originalSubmodelpart.Name}')
                print(f'originalSubmodelpart.Elements {len(originalSubmodelpart.Elements)}')
                print(f'originalSubmodelpart.Conditions {len(originalSubmodelpart.Conditions)}')
                for originalNode in originalSubmodelpart.Nodes:
                    if originalNode in hyper_reduced_model_part_help.Nodes:
                        hyperReducedSubmodelpart.AddNode(originalNode,0)
                ## More eficient way to implement this is possible
                for originalElement in originalSubmodelpart.Elements:
                    for key in HR_data["Elements"].keys():
                        if originalElement.Id == int(key)+1:
                            hyperReducedSubmodelpart.AddElement(originalElement,0)
                            print(f'For the submodelpart {hyperReducedSubmodelpart.Name}, the element with the Id {originalElement.Id} is assigned the key {key}')
                for originalCondition in originalSubmodelpart.Conditions:
                    for key in HR_data["Conditions"].keys():
                        if originalCondition.Id == int(key)+1:
                            hyperReducedSubmodelpart.AddCondition(originalCondition,0)
                            print(f'For the submodelpart {hyperReducedSubmodelpart.Name}, the condition with the Id {originalCondition.Id} is assigned the key {key}')


        # Building the VISUALIZE_HROM submodel part
        print('Adding skin for visualization...')
        hyper_reduced_model_part2 = HROM_Model_Part.CreateSubModelPart("VISUALIZE_HROM")
        for condition in computing_model_part.Conditions:
            for node in condition.GetNodes():
                hyper_reduced_model_part2.AddNode(node, 0)
            hyper_reduced_model_part2.AddCondition(condition, 0)
        for node in computing_model_part.Nodes:
            hyper_reduced_model_part2.AddNode(node, 0)

        ## Creating the mdpa file using ModelPartIO object
        print('About to print ...')
        KratosMultiphysics.ModelPartIO("Hyper_Reduced_Model_Part", KratosMultiphysics.IO.WRITE| KratosMultiphysics.IO.MESH_ONLY ).WriteModelPart(HROM_Model_Part)
        print('\nHyper_Reduced_Model_Part.mdpa created!\n')
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting("Hyper_Reduced_Model_Part.time")


