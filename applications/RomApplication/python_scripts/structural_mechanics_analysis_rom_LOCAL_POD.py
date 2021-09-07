import KratosMultiphysics
import KratosMultiphysics.RomApplication as romapp
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.RomApplication.empirical_cubature_method_manySets import EmpiricalCubatureMethod
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom as solver_wrapper
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

import json
import numpy as np


class StructuralMechanicsAnalysisROM(StructuralMechanicsAnalysis):

    def __init__(self,model,project_parameters, hyper_reduction_element_selector = None):
        super().__init__(model,project_parameters)
        if hyper_reduction_element_selector != None :
            if hyper_reduction_element_selector == "EmpiricalCubature":
                self.hyper_reduction_element_selector = EmpiricalCubatureMethod()
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
        super().ModifyAfterSolverInitialize()
        computing_model_part = self._solver.GetComputingModelPart()
        with open('RomParameters.json') as f:
            data = json.load(f)
            nodal_dofs = len(data["rom_settings"]["nodal_unknowns"])
            nodal_modes = data["nodal_modes"]
            number_of_bases = len(nodal_modes['1'])
            POD_BASES = []
            rom_dofs=[]
            for i in range(number_of_bases):
                POD_BASES.append(romapp.RomBasis())
            for i in range(number_of_bases):
                rom_dofs.append(self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"][i].GetInt())
            for node in computing_model_part.Nodes:
                for k in range(number_of_bases):
                    aux = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs[k])
                    for j in range(nodal_dofs):
                        Counter=str(node.Id)
                        for i in range(rom_dofs[k]):
                            aux[j,i] = nodal_modes[Counter][str(k)][j][i]
                    POD_BASES[k].SetNodalBasis(node.Id, aux)
            distance_to_clusters = romapp.DistanceToClusters(number_of_bases)
            w = data["w"]
            z0 = data["z0"]
            for i in range(number_of_bases):
                for j in range(number_of_bases):
                    distance_to_clusters.SetZEntry(z0[f'{i}'][f'{j}'],i,j)
                    for k in range(number_of_bases):
                        entries = len(w[f'{i}'][f'{j}'][f'{k}'])
                        temp_vector = KratosMultiphysics.Vector(entries)
                        for entry in range(entries):
                            temp_vector[entry] = w[f'{i}'][f'{j}'][f'{k}'][entry]
                        distance_to_clusters.SetWEntry(temp_vector,i,j,k)
        Bases = romapp.RomBases()
        for i in range(number_of_bases):
            Bases.AddBasis(i, POD_BASES[i])
        self.Number_Of_Clusters = number_of_bases
        self._GetSolver().get_builder_and_solver().SetUpBases(Bases)
        self._GetSolver().get_builder_and_solver().SetUpDistances(distance_to_clusters)

        self.step_index = 0

        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                # for i in range(self.Number_Of_Clusters):
                #     self.time_step_residual_matrix_container.append([])
                self.ResidualUtilityObject = romapp.RomResidualsUtility(self._GetSolver().GetComputingModelPart(), self.project_parameters["solver_settings"]["rom_settings"], self._GetSolver().get_solution_scheme(), Bases)


    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self._GetSolver().get_builder_and_solver().UpdateCurrentCluster()


    def FinalizeSolutionStep(self):
        self._GetSolver().get_builder_and_solver().UpdateZMatrix()
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                print('\n\n\n\nGenerating matrix of residuals')
                ResMat = self.ResidualUtilityObject.GetResiduals(self._GetSolver().get_builder_and_solver().GetCurrentCluster())
                #print(np.shape(ResMat))
                NP_ResMat = np.array(ResMat, copy=False)
                #self.time_step_residual_matrix_container[self._GetSolver().get_builder_and_solver().GetCurrentCluster()].append(NP_ResMat)
                self.time_step_residual_matrix_container.append(NP_ResMat)
        super().FinalizeSolutionStep()


    def Finalize(self):
        super().Finalize()
        if self.hyper_reduction_element_selector != None:
            if self.hyper_reduction_element_selector.Name == "EmpiricalCubature":
                OriginalNumberOfElements = self._GetSolver().GetComputingModelPart().NumberOfElements()
                ModelPartName = self._GetSolver().settings["model_import_settings"]["input_filename"].GetString()
                LeftSingularVectorsOfResidualProjected = self._ObtainBasis()
                self. hyper_reduction_element_selector.SetUp(LeftSingularVectorsOfResidualProjected, OriginalNumberOfElements, ModelPartName)
                self.hyper_reduction_element_selector.Run()
                # for cluster_number in range(self.Number_Of_Clusters):
                #     self. hyper_reduction_element_selector.SetUp(self.time_step_residual_matrix_container[cluster_number], OriginalNumberOfElements, ModelPartName, cluster_number)
                #     self.hyper_reduction_element_selector.Run()


    def _ObtainBasis(self):
        ### Building the Snapshot matrix ####
        for i in range (len(self.time_step_residual_matrix_container)):
            if i == 0:
                SnapshotMatrix = self.time_step_residual_matrix_container[i]
            else:
                SnapshotMatrix = np.c_[SnapshotMatrix,self.time_step_residual_matrix_container[i]]
            self.time_step_residual_matrix_container[i] = 0
        U,_,_,_ = RandomizedSingularValueDecomposition().Calculate(SnapshotMatrix, 1e-6)
        return U.T




