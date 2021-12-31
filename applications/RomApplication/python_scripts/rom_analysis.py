import json
import importlib
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics.RomApplication import python_solvers_wrapper_rom
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

def CreateRomAnalysisInstance(cls, global_model, parameters):
    class RomAnalysis(cls):

        def __init__(self,global_model, parameters):
            super().__init__(global_model, parameters)

        def _CreateSolver(self):
            """ Create the Solver (and create and import the ModelPart if it is not alread in the model) """

            # Get the ROM settings from the RomParameters.json input file and set them in the "solver_settings" of the solver introducing the physics
            with open('RomParameters.json') as rom_parameters:
                rom_settings = KratosMultiphysics.Parameters(rom_parameters.read())
                self.project_parameters["solver_settings"].AddValue("rom_settings", rom_settings["rom_settings"])

            # Create the ROM solver
            return python_solvers_wrapper_rom.CreateSolverByParameters(
                self.model,
                self.project_parameters["solver_settings"],
                self.project_parameters["problem_data"]["parallel_type"].GetString())

        def _GetSimulationName(self):
            return "::[ROM Simulation]:: "

        def ModifyAfterSolverInitialize(self):
            """Here is where the ROM_BASIS is imposed to each node"""
            super().ModifyAfterSolverInitialize()

            # Get the model part where the ROM is to be applied
            computing_model_part = self._GetSolver().GetComputingModelPart()

            # Set ROM basis
            with open('RomParameters.json') as f:
                # Get the ROM data from RomParameters.json
                data = json.load(f)
                nodal_modes = data["nodal_modes"]
                nodal_dofs = len(data["rom_settings"]["nodal_unknowns"])
                rom_dofs = self.project_parameters["solver_settings"]["rom_settings"]["number_of_rom_dofs"].GetInt()

                # Set the nodal ROM basis
                aux = KratosMultiphysics.Matrix(nodal_dofs, rom_dofs)
                for node in computing_model_part.Nodes:
                    node_id = str(node.Id)
                    for j in range(nodal_dofs):
                        for i in range(rom_dofs):
                            aux[j,i] = nodal_modes[node_id][j][i]
                    node.SetValue(KratosROM.ROM_BASIS, aux)

                self.train_hrom = data["train_hrom"]
                if self.train_hrom:
                    self.hyper_reduction_element_selector = EmpiricalCubatureMethod()
                    self.time_step_residual_matrix_container = []
                    self.ECM_SVD_tolerance = data["ECM_SVD_tolerance"]
                if data["run_hrom"]:
                    for key in data["elements_and_weights"]["Elements"].keys():
                        computing_model_part.GetElement(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, data["elements_and_weights"]["Elements"][key])
                    for key in data["elements_and_weights"]["Conditions"].keys():
                        computing_model_part.GetCondition(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, data["elements_and_weights"]["Conditions"][key])

            # Hyper-reduction
            if self.train_hrom:
                self.ResidualUtilityObject = KratosROM.RomResidualsUtility(
                    computing_model_part,
                    self.project_parameters["solver_settings"]["rom_settings"],
                    self._GetSolver()._GetScheme())

        def FinalizeSolutionStep(self):
            if self.train_hrom:
                KratosMultiphysics.Logger.PrintInfo("RomAnalysis","Generating matrix of residuals.")
                res_mat = self.ResidualUtilityObject.GetResiduals()
                np_res_mat = np.array(res_mat, copy=False)
                self.time_step_residual_matrix_container.append(np_res_mat)

            super().FinalizeSolutionStep()

        def Finalize(self):
            super().Finalize()

            if self.train_hrom:
                self. hyper_reduction_element_selector.SetUp(self._ObtainBasis())
                self.hyper_reduction_element_selector.Run()
                self._AppendElementsAndWeightsToRomParameters()
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

        def _AppendElementsAndWeightsToRomParameters(self):
            original_number_of_elements = self._GetSolver().GetComputingModelPart().NumberOfElements()
            w = np.squeeze(self.hyper_reduction_element_selector.w)
            z = self.hyper_reduction_element_selector.z
            ### Saving Elements and conditions
            elements_and_weights = {}
            elements_and_weights["Elements"] = {}
            elements_and_weights["Conditions"] = {}
            #Only one element found !
            if type(z)==np.int64 or type(z)==np.int32:
                if z <=  original_number_of_elements-1:
                    elements_and_weights["Elements"][int(z)] = (float(w))
                else:
                    elements_and_weights["Conditions"][int(z)- original_number_of_elements] = (float(w))
            #Many elements found
            else:
                for j in range (0,len(z)):
                    if z[j] <=  original_number_of_elements -1:
                        elements_and_weights["Elements"][int(z[j])] = (float(w[j]))
                    else:
                        elements_and_weights["Conditions"][int(z[j])- original_number_of_elements] = (float(w[j]))

            with open('RomParameters.json','r') as f:
                updated_rom_parameters = json.load(f)
            with open('RomParameters.json','w') as f:
                updated_rom_parameters["train_hrom"] = False
                updated_rom_parameters["run_hrom"] = True
                updated_rom_parameters["elements_and_weights"] = elements_and_weights
                json.dump(updated_rom_parameters,f, indent = 4)
            print('\n\n RomParameters.json file was updated to include HROM info\n\n')


        def _CreateHyperReducedModelPart(self):
            ModelPartName = self._GetSolver().settings["model_import_settings"]["input_filename"].GetString()
            current_model = KratosMultiphysics.Model()
            computing_model_part = current_model.CreateModelPart("main")
            model_part_io = KratosMultiphysics.ModelPartIO(ModelPartName)
            model_part_io.ReadModelPart(computing_model_part)
            hyper_reduced_model_part_help =   current_model.CreateModelPart("Helping")

            #TODO optimize implementation
            with open('RomParameters.json') as f:
                HR_data = json.load(f)
                for key in HR_data["elements_and_weights"]["Elements"].keys():
                    for node in computing_model_part.GetElement(int(key)+1).GetNodes():
                        hyper_reduced_model_part_help.AddNode(node,0)
                for key in HR_data["elements_and_weights"]["Conditions"].keys():
                    for node in computing_model_part.GetCondition(int(key)+1).GetNodes():
                        hyper_reduced_model_part_help.AddNode(node,0)

            # The HROM model part. It will include two sub-model parts. One for caculation, another one for visualization
            HROM_Model_Part =  current_model.CreateModelPart("HROM_Model_Part")

            # Building the COMPUTE_HROM submodel part
            hyper_reduced_model_part = HROM_Model_Part.CreateSubModelPart("COMPUTE_HROM")

            # TODO implement the hyper-reduced model part creation in C++
            with open('RomParameters.json') as f:
                data = json.load(f)
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
                        for key in data["elements_and_weights"]["Elements"].keys():
                            if originalElement.Id == int(key)+1:
                                hyperReducedSubmodelpart.AddElement(originalElement,0)
                                print(f'For the submodelpart {hyperReducedSubmodelpart.Name}, the element with the Id {originalElement.Id} is assigned the key {key}')
                    for originalCondition in originalSubmodelpart.Conditions:
                        for key in data["elements_and_weights"]["Conditions"].keys():
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