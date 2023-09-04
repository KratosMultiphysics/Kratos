# Import Python modules
import json
import numpy as np
from pathlib import Path

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics.RomApplication.empirical_cubature_method import EmpiricalCubatureMethod
from KratosMultiphysics.RomApplication.randomized_singular_value_decomposition import RandomizedSingularValueDecomposition

class HRomTrainingUtility(object):
    """Auxiliary utility for the HROM training.
    This class encapsulates all the functions required for the HROM training.
    These are the ROM residuals collection, the calculation and storage of
    the HROM weights and the creation and export of the HROM model parts.
    """

    def __init__(self, solver, custom_settings):
        # Validate and assign the HROM training settings
        settings = custom_settings["hrom_settings"]
        settings.ValidateAndAssignDefaults(self.__GetHRomTrainingDefaultSettings())

        # Projection strategy (residual projection)
        settings["projection_strategy"] = custom_settings["projection_strategy"] # To be consistent with the solving strategy

        # Create the HROM element selector
        element_selection_type = settings["element_selection_type"].GetString()
        if element_selection_type == "empirical_cubature":
            self.hyper_reduction_element_selector = EmpiricalCubatureMethod()
            self.element_selection_svd_truncation_tolerance = settings["element_selection_svd_truncation_tolerance"].GetDouble()
        else:
            err_msg = "\'{}\' HROM \'element_selection_type\' is not supported.".format(element_selection_type)
            raise Exception(err_msg)

        # Auxiliary member variables
        self.solver = solver
        self.time_step_residual_matrix_container = []
        self.echo_level = settings["echo_level"].GetInt()
        self.rom_settings = custom_settings["rom_settings"]
        self.hrom_visualization_model_part = settings["create_hrom_visualization_model_part"].GetBool()
        self.projection_strategy = settings["projection_strategy"].GetString()
        self.hrom_output_format = settings["hrom_format"].GetString()
        self.include_minimum_condition = settings["include_minimum_condition"].GetBool()
        self.include_condition_parents = settings["include_condition_parents"].GetBool()

        # Retrieve list of model parts from settings
        self.include_conditions_model_parts_list = settings["include_conditions_model_parts_list"].GetStringArray()
        self.include_elements_model_parts_list = settings["include_elements_model_parts_list"].GetStringArray()
        self.include_nodal_neighbouring_elements_model_parts_list = settings["include_nodal_neighbouring_elements_model_parts_list"].GetStringArray()

        # Check if the model parts exist
        for model_part_name in self.include_conditions_model_parts_list:
            if not self.solver.model.HasModelPart(model_part_name):
                raise Exception('The model part named "' + model_part_name + '" does not exist in the model')

        for model_part_name in self.include_elements_model_parts_list:
            if not self.solver.model.HasModelPart(model_part_name):
                raise Exception('The model part named "' + model_part_name + '" does not exist in the model')

        # Rom settings files
        self.rom_basis_output_name = Path(custom_settings["rom_basis_output_name"].GetString())
        self.rom_basis_output_folder = Path(custom_settings["rom_basis_output_folder"].GetString())

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
                self.solver._GetScheme())

            if self.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("HRomTrainingUtility","RomResidualsUtility created.")

        # Generate the matrix of projected residuals
        if self.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("HRomTrainingUtility","Generating matrix of projected residuals.")
        if (self.projection_strategy=="galerkin"):
                res_mat = self.__rom_residuals_utility.GetProjectedResidualsOntoPhi()
        elif (self.projection_strategy=="petrov_galerkin"):
                res_mat = self.__rom_residuals_utility.GetProjectedResidualsOntoPsi()
        else:
            err_msg = f"Projection strategy \'{self.projection_strategy}\' for HROM is not supported."
            raise Exception(err_msg)

        np_res_mat = np.array(res_mat, copy=False)
        self.time_step_residual_matrix_container.append(np_res_mat)

    def CalculateAndSaveHRomWeights(self):
        # Calculate the residuals basis and compute the HROM weights from it
        residual_basis = self.__CalculateResidualBasis()
        n_conditions = self.solver.GetComputingModelPart().NumberOfConditions() # Conditions must be included as an extra restriction to enforce ECM to capture all BC's regions.
        self.hyper_reduction_element_selector.SetUp(residual_basis, constrain_sum_of_weights=True, constrain_conditions = False, number_of_conditions = n_conditions)
        self.hyper_reduction_element_selector.Run()

        # Save the HROM weights in the RomParameters.json
        # Note that in here we are assuming this naming convention for the ROM json file
        self.AppendHRomWeightsToRomParameters()

    def CreateHRomModelParts(self):
        # Get solver data
        model_part_name = self.solver.settings["model_part_name"].GetString()
        model_part_output_name = self.solver.settings["model_import_settings"]["input_filename"].GetString()
        # computing_model_part = self.solver.GetComputingModelPart()
        #computing_model_part = self.solver.GetComputingModelPart().GetRootModelPart() #TODO: DECIDE WHICH ONE WE SHOULD USE?¿?¿ MOST PROBABLY THE ROOT FOR THOSE CASES IN WHICH THE COMPUTING IS CUSTOM (e.g. CFD)
        aux_model = KratosMultiphysics.Model()
        computing_model_part = aux_model.CreateModelPart("main")
        model_part_io = KratosMultiphysics.ModelPartIO(model_part_output_name)
        model_part_io.ReadModelPart(computing_model_part)

        # Create a new model with the HROM main model part
        # This is intentionally done in order to completely emulate the origin model part
        aux_model = KratosMultiphysics.Model()
        hrom_main_model_part = aux_model.CreateModelPart(model_part_name)

        if self.hrom_output_format == "numpy":
            hrom_info = KratosMultiphysics.Parameters(json.JSONEncoder().encode(self.__CreateDictionaryWithRomElementsAndWeights()))
        elif self.hrom_output_format == "json":
            with (self.rom_basis_output_folder / self.rom_basis_output_name).with_suffix('.json').open('r') as f:
                rom_parameters = KratosMultiphysics.Parameters(f.read())
                hrom_info = rom_parameters["elements_and_weights"]

        # Get the weights and fill the HROM computing model part
        KratosROM.RomAuxiliaryUtilities.SetHRomComputingModelPart(hrom_info,computing_model_part,hrom_main_model_part)
        if self.echo_level > 0:
            KratosMultiphysics.Logger.PrintInfo("HRomTrainingUtility","HROM computing model part \'{}\' created.".format(hrom_main_model_part.FullName()))

        # Output the HROM model part in mdpa format
        hrom_output_name = "{}HROM".format(model_part_output_name)
        model_part_io = KratosMultiphysics.ModelPartIO(hrom_output_name, KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY | KratosMultiphysics.IO.SCIENTIFIC_PRECISION)
        model_part_io.WriteModelPart(hrom_main_model_part)
        KratosMultiphysics.kratos_utilities.DeleteFileIfExisting("{}.time".format(hrom_output_name))
        if self.echo_level > 0:
            KratosMultiphysics.Logger.PrintInfo("HRomTrainingUtility","HROM mesh written in \'{}.mdpa\'".format(hrom_output_name))

        #TODO: Make this optional
        #TODO: Move this out of here
        # Create the HROM visualization model parts
        if self.hrom_visualization_model_part:
            # Create the HROM visualization mesh from the origin model part
            hrom_visualization_model_part_name = "{}Visualization".format(hrom_main_model_part.Name)
            hrom_visualization_model_part = aux_model.CreateModelPart(hrom_visualization_model_part_name)
            KratosROM.RomAuxiliaryUtilities.SetHRomVolumetricVisualizationModelPart(computing_model_part, hrom_visualization_model_part)
            if self.echo_level > 0:
                KratosMultiphysics.Logger.PrintInfo("HRomTrainingUtility","HROM visualization model part \'{}\' created.".format(hrom_visualization_model_part.FullName()))

            print(hrom_visualization_model_part)

            # Write the HROM visualization mesh
            hrom_vis_output_name = "{}HROMVisualization".format(model_part_output_name)
            model_part_io = KratosMultiphysics.ModelPartIO(hrom_vis_output_name, KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY | KratosMultiphysics.IO.SCIENTIFIC_PRECISION)
            model_part_io.WriteModelPart(hrom_visualization_model_part)
            KratosMultiphysics.kratos_utilities.DeleteFileIfExisting("{}.time".format(hrom_vis_output_name))
            if self.echo_level > 0:
                KratosMultiphysics.Logger.PrintInfo("HRomTrainingUtility","HROM visualization mesh written in \'{}.mdpa\'".format(hrom_vis_output_name))

    @classmethod
    def __GetHRomTrainingDefaultSettings(cls):
        default_settings = KratosMultiphysics.Parameters("""{
            "hrom_format": "numpy",
            "element_selection_type": "empirical_cubature",
            "element_selection_svd_truncation_tolerance": 1.0e-6,
            "echo_level" : 0,
            "create_hrom_visualization_model_part" : true,
            "projection_strategy": "galerkin",
            "include_conditions_model_parts_list": [],
            "include_elements_model_parts_list": [],
            "include_nodal_neighbouring_elements_model_parts_list":[],
            "include_minimum_condition": false,
            "include_condition_parents": true
        }""")
        return default_settings

    def __CalculateResidualBasis(self):
        # Set up the residual snapshots matrix
        residuals_snapshot_matrix = self._GetResidualsProjectedMatrix()

        # Calculate the randomized and truncated SVD of the residual snapshots
        u,_,_,_ = RandomizedSingularValueDecomposition(COMPUTE_V=False).Calculate(
            residuals_snapshot_matrix,
            self.element_selection_svd_truncation_tolerance)

        return u



    def _GetResidualsProjectedMatrix(self):
        # Set up the residual snapshots matrix
        n_steps = len(self.time_step_residual_matrix_container)
        residuals_snapshot_matrix = self.time_step_residual_matrix_container[0]
        for i in range(1,n_steps):
            del self.time_step_residual_matrix_container[0] # Avoid having two matrices, numpy does not concatenate references.
            residuals_snapshot_matrix = np.c_[residuals_snapshot_matrix,self.time_step_residual_matrix_container[0]]
        return residuals_snapshot_matrix



    def AppendHRomWeightsToRomParameters(self):
        number_of_elements = self.solver.GetComputingModelPart().GetRootModelPart().NumberOfElements()
        weights = np.squeeze(self.hyper_reduction_element_selector.w)
        indexes = self.hyper_reduction_element_selector.z

        # Create dictionary with HROM weights (Only used for the expansion of the selected Conditions to include their parent Elements)
        hrom_weights = self.__CreateDictionaryWithRomElementsAndWeights(weights,indexes,number_of_elements)

        # Get the root model part
        root_model_part = self.solver.GetComputingModelPart().GetRootModelPart()

        # Add conditions
        for model_part_name in self.include_conditions_model_parts_list:
            # Check if the sub model part exists
            if self.solver.model.HasModelPart(model_part_name):
                conditions_to_include_model_part = self.solver.model.GetModelPart(model_part_name)

                # Call the GetConditionIdsNotInHRomModelPart function
                new_conditions = KratosROM.RomAuxiliaryUtilities.GetConditionIdsNotInHRomModelPart(
                    root_model_part, # The complete model part
                    conditions_to_include_model_part, # The model part containing the conditions to be included
                    hrom_weights)

                # Add the new conditions to the conditions dict with a null weight
                for condition_id in new_conditions:
                    hrom_weights["Conditions"][condition_id] = 0.0

                # If needed, update your weights and indexes using __AddSelectedConditionsWithZeroWeights function with the new_conditions
                weights, indexes = self.__AddSelectedConditionsWithZeroWeights(weights, indexes, new_conditions, number_of_elements)

        # Add elements
        for model_part_name in self.include_elements_model_parts_list:
            # Check if the sub model part exists
            if self.solver.model.HasModelPart(model_part_name):
                elements_to_include_model_part = self.solver.model.GetModelPart(model_part_name)

                # Call the GetElementIdsNotInHRomModelPart function
                new_elements = KratosROM.RomAuxiliaryUtilities.GetElementIdsNotInHRomModelPart(
                    root_model_part, # The complete model part
                    elements_to_include_model_part, # The model part containing the elements to be included
                    hrom_weights)

                # Add the new elements to the elements dict with a null weight
                for element_id in new_elements:
                    hrom_weights["Elements"][element_id] = 0.0

                # If needed, update your weights and indexes using __AddSelectedElementsWithZeroWeights function with the new_elements
                weights, indexes = self.__AddSelectedElementsWithZeroWeights(weights, indexes, new_elements)
        
        # Add nodal neighbouring elements
        for model_part_name in self.include_nodal_neighbouring_elements_model_parts_list:
            # Check if the sub model part exists
            if self.solver.model.HasModelPart(model_part_name):
                nodal_neighbours_model_part = self.solver.model.GetModelPart(model_part_name)

                # Call the GetNodalNeighbouringElementIdsNotInHRom function
                new_nodal_neighbours = KratosROM.RomAuxiliaryUtilities.GetNodalNeighbouringElementIdsNotInHRom(
                    root_model_part, # The complete model part
                    nodal_neighbours_model_part, # The model part containing the nodal neighbouring elements to be included
                    hrom_weights)

                # Add the new nodal neighbouring elements to the elements dict with a null weight
                for element_id in new_nodal_neighbours:
                    hrom_weights["Elements"][element_id] = 0.0

                # If needed, update your weights and indexes using __AddSelectedElementsWithZeroWeights function with the new_nodal_neighbours
                weights, indexes = self.__AddSelectedElementsWithZeroWeights(weights, indexes, new_nodal_neighbours)

        # If required, keep at least one condition per submodelpart
        # This might be required by those BCs involving the faces (e.g. slip BCs)
        if self.include_minimum_condition:
            # Get the HROM conditions to be added
            minimum_conditions = KratosROM.RomAuxiliaryUtilities.GetHRomMinimumConditionsIds(
                self.solver.GetComputingModelPart().GetRootModelPart(), #TODO: I think this one should be the root
                hrom_weights["Conditions"])

            # Add the selected conditions to the conditions dict with a null weight
            for cond_id in minimum_conditions:
                hrom_weights["Conditions"][cond_id] = 0.0
            weights, indexes = self.__AddSelectedConditionsWithZeroWeights(weights, indexes, minimum_conditions, number_of_elements)

        # If required, add the HROM conditions parent elements
        # Note that we add these with zero weight so their future assembly will have no effect
        if self.include_condition_parents:
            # Get the HROM condition parents from the current HROM weights
            missing_condition_parents = KratosROM.RomAuxiliaryUtilities.GetHRomConditionParentsIds(
                self.solver.GetComputingModelPart().GetRootModelPart(), #TODO: I think this one should be the root
                hrom_weights)

            # Add the missing parents to the elements dict with a null weight
            for parent_id in missing_condition_parents:
                hrom_weights["Elements"][parent_id] = 0.0
            weights, indexes = self.__AddSelectedElementsWithZeroWeights(weights,indexes, missing_condition_parents)

        if self.hrom_output_format=="numpy":
            element_indexes = np.where(indexes < number_of_elements)[0]
            condition_indexes = np.where(indexes >= number_of_elements)[0]
            np.save(self.rom_basis_output_folder / "HROM_ElementWeights.npy", weights[element_indexes])
            np.save(self.rom_basis_output_folder / "HROM_ConditionWeights.npy", weights[condition_indexes])
            np.save(self.rom_basis_output_folder / "HROM_ElementIds.npy", indexes[element_indexes])  # FIXME fix the -1 in the indexes of numpy and ids of Kratos
            np.save(self.rom_basis_output_folder / "HROM_ConditionIds.npy", indexes[condition_indexes] - number_of_elements)  # FIXME fix the -1 in the indexes of numpy and ids of Kratos

        elif self.hrom_output_format=="json":
            with (self.rom_basis_output_folder / self.rom_basis_output_name).with_suffix('.json').open('r') as f:
                updated_rom_parameters = json.load(f)
                updated_rom_parameters["elements_and_weights"] = hrom_weights  # TODO: Rename elements_and_weights to hrom_weights

            with (self.rom_basis_output_folder / self.rom_basis_output_name).with_suffix('.json').open('w') as f:
                json.dump(updated_rom_parameters, f, indent=4)

        if self.echo_level > 0 : KratosMultiphysics.Logger.PrintInfo("HRomTrainingUtility","\'RomParameters.json\' file updated with HROM weights.")

    def __AddSelectedElementsWithZeroWeights(self, original_weights,original_elements, elements_to_add):

        added_elements_ids = np.array(elements_to_add)
        updated_elements = np.r_[original_elements, added_elements_ids]
        updated_weights = np.r_[original_weights, np.zeros(added_elements_ids.shape)]

        return updated_weights, updated_elements

    def __AddSelectedConditionsWithZeroWeights(self, original_weights, original_conditions, conditions_to_add, number_of_elements):

        added_conditions_ids = np.array(conditions_to_add)+number_of_elements
        updated_conditions = np.r_[original_conditions, added_conditions_ids]
        updated_weights = np.r_[original_weights, np.zeros(added_conditions_ids.shape)]

        return updated_weights, updated_conditions

    def __CreateDictionaryWithRomElementsAndWeights(self, weights = None, indexes=None, number_of_elements = None):

        if number_of_elements is None:
            number_of_elements = self.solver.GetComputingModelPart().NumberOfElements()
        if weights is None:
            weights = np.r_[np.load(self.rom_basis_output_folder / "HROM_ElementWeights.npy"), np.load(self.rom_basis_output_folder / "HROM_ConditionWeights.npy")]
        if indexes is None:
            indexes = np.r_[np.load(self.rom_basis_output_folder / "HROM_ElementIds.npy"), np.load(self.rom_basis_output_folder / "HROM_ConditionIds.npy") + number_of_elements]
        hrom_weights = {}
        hrom_weights["Elements"] = {}
        hrom_weights["Conditions"] = {}

        if type(indexes)==np.int64 or type(indexes)==np.int32:
            # Only one element found !
            if indexes <= number_of_elements-1:
                hrom_weights["Elements"][int(indexes)] = float(weights)
            else:
                hrom_weights["Conditions"][int(indexes)-number_of_elements] = float(weights)
        else:
            # Many elements found
            for j in range (len(indexes)):
                if indexes[j] <=  number_of_elements -1:
                    hrom_weights["Elements"][int(indexes[j])] = float(weights[j])
                else:
                    hrom_weights["Conditions"][int(indexes[j])-number_of_elements] = float(weights[j])

        return hrom_weights

