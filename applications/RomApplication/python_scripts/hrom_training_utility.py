# Import Python modules
import json
import numpy as np

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

class HRomTrainingUtility(object):
    """Auxiliary utility for the HROM training.




    This class serves as a generic class to make a standard Kratos solver ROM-compatible.
    It extends the default parameters to include the \'rom_settings\' and overrides the
    creation of the builder and solver to use the ROM one.
    """

    def __init__(self, solver, custom_settings):
        # Validate and assign the HROM training settings
        settings = custom_settings["hrom_settings"].ValidateAndAssignDefaults(self.__GetHRomTrainingDefaultSettings())

        # Create the HROM element selector
        element_selection_type = settings["element_selection_type"].GetString()
        if element_selection_type == "empirical_cubature":
            self.hyper_reduction_element_selector = EmpiricalCubatureMethod()
            self.element_selection_svd_truncation_tolerance = settings["element_selection_svd_truncation_tolerance"].GetBool()
        else:
            err_msg = "\'{}\' HROM \'element_selection_type\' is not supported.".format(element_selection_type)
            raise Exception(err_msg)

        # Auxiliary member variables
        self.solver = solver
        self.rom_settings = settings["rom_settings"]
        self.echo_level = settings["echo_level"].GetInt()
        self.time_step_residual_matrix_container = []

    def AppendCurrentStepResiduals(self):
        # Get the computing model part from the solver implementing the problem physics
        computing_model_part = self.solver.GetComputingModelPart()

        # If not created yet, create the ROM residuals utility
        # Note that this ensures that the residuals utility is created in the first residuals append call
        # If not, it might happen that the solver scheme is created by the ROM residuals call rather than by the solver one
        if not hasattr(self, '__rom_residuals_utility'):
            self.__rom_residuals_utility = KratosROM.RomResidualsUtility(
                computing_model_part,
                self.rom_settings,
                self.solver.GetScheme())

            if self.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("HRomTrainingUtility","RomResidualsUtility created.")

        # Generate the matrix of residuals
        if self.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("HRomTrainingUtility","Generating matrix of residuals.")
        res_mat = self.__rom_residuals_utility.GetResiduals()
        np_res_mat = np.array(res_mat, copy=False)
        self.time_step_residual_matrix_container.append(np_res_mat)

    def CalculateAndSaveHRomWeights(self):
        # Calculate the residuals basis and compute the HROM weights from it
        residual_basis = self.__CalculateResidualBasis()
        self.hyper_reduction_element_selector.SetUp(residual_basis)
        self.hyper_reduction_element_selector.Run()

        # Save the HROM weights in the RomParameters.json
        # Note that in here we are assuming this naming convention for the ROM json file
        self.__AppendHRomWeightsToRomParameters()

    def CreateHRomModelParts(self):
        # Get model from solver
        current_model = self.solver.model
        model_part_name = self.solver.settings["model_import_settings"]["input_filename"].GetString()
        computing_model_part = self.solver.GetComputingModelPart()
        hyper_reduced_model_part_help = current_model.CreateModelPart("Helping")

        #FIXME: THIS IS A MESS. IMPLEMENT IT IN C++ USING A PointerVectorSet
        #FIXME: THIS IS A MESS. IMPLEMENT IT IN C++ USING A PointerVectorSet
        #FIXME: THIS IS A MESS. IMPLEMENT IT IN C++ USING A PointerVectorSet
        #FIXME: THIS IS A MESS. IMPLEMENT IT IN C++ USING A PointerVectorSet
        #FIXME: THIS IS A MESS. IMPLEMENT IT IN C++ USING A PointerVectorSet
        #FIXME: THIS IS A MESS. IMPLEMENT IT IN C++ USING A PointerVectorSet

        #TODO: optimize implementation
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

        # TODO: implement the hyper-reduced model part creation in C++
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

        #FIXME: THIS IS A MESS. IMPLEMENT IT IN C++ USING A PointerVectorSet
        #FIXME: THIS IS A MESS. IMPLEMENT IT IN C++ USING A PointerVectorSet
        #FIXME: THIS IS A MESS. IMPLEMENT IT IN C++ USING A PointerVectorSet
        #FIXME: THIS IS A MESS. IMPLEMENT IT IN C++ USING A PointerVectorSet
        #FIXME: THIS IS A MESS. IMPLEMENT IT IN C++ USING A PointerVectorSet
        #FIXME: THIS IS A MESS. IMPLEMENT IT IN C++ USING A PointerVectorSet

    @classmethod
    def __GetHRomTrainingDefaultSettings(cls):
        default_settings = KratosMultiphysics.Parameters("""{
            "element_selection_type": "empirical_cubature",
            "element_selection_svd_truncation_tolerance": 1.0e-6,
            "echo_level" : 0
        }""")
        return default_settings

    def __CalculateResidualBasis(self):
        # Set up the residual snapshots matrix
        n_steps = len(self.time_step_residual_matrix_container)
        residuals_snapshot_matrix = self.time_step_residual_matrix_container[0]
        for i in range(1,n_steps):
            residuals_snapshot_matrix = np.c_[residuals_snapshot_matrix,self.time_step_residual_matrix_container[i]]

        # Calculate the randomized and truncated SVD of the residual snapshots
        u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(
            residuals_snapshot_matrix,
            self.element_selection_svd_truncation_tolerance)

        return u

    def __AppendHRomWeightsToRomParameters(self):
        n_elements = self.solver.GetComputingModelPart().NumberOfElements()
        w = np.squeeze(self.hyper_reduction_element_selector.w)
        z = self.hyper_reduction_element_selector.z

        # Create dictionary with HROM weights
        hrom_weights = {}
        hrom_weights["Elements"] = {}
        hrom_weights["Conditions"] = {}

        if type(z)==np.int64 or type(z)==np.int32:
            # Only one element found !
            if z <= n_elements-1:
                hrom_weights["Elements"][int(z)] = float(w)
            else:
                hrom_weights["Conditions"][int(z)-n_elements] = float(w)
        else:
            # Many elements found
            for j in range (0,len(z)):
                if z[j] <=  n_elements -1:
                    hrom_weights["Elements"][int(z[j])] = float(w[j])
                else:
                    hrom_weights["Conditions"][int(z[j])-n_elements] = float(w[j])

        # Append weights to RomParameters.json
        # We first parse the current RomParameters.json to then append and edit the data
        with open('RomParameters.json','r') as f:
            updated_rom_parameters = json.load(f)
            updated_rom_parameters["train_hrom"] = False
            updated_rom_parameters["run_hrom"] = True
            updated_rom_parameters["hrom_weights"] = hrom_weights

        with open('RomParameters.json','w') as f:
            json.dump(updated_rom_parameters,f, indent = 4)

        if self.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("HRomTrainingUtility","\'RomParameters.json\' file updated with HROM weights.")
