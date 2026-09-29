import os
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
from KratosMultiphysics.RomApplication.hrom_training_utility import HRomTrainingUtility

if kratos_utilities.CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication

@KratosUnittest.skipIfApplicationsNotAvailable("FluidDynamicsApplication")
class TestHromTrainingUtilityRom(KratosUnittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.work_folder = "hrom_training_utility_test_files/"
        parameters_filename = "ProjectParameters.json"

        with KratosUnittest.WorkFolderScope(cls.work_folder, __file__):
            with open(parameters_filename, 'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            cls.model = KratosMultiphysics.Model()
            cls.simulation = rom_testing_utilities.SetUpSimulationInstance(cls.model, parameters)
            cls.simulation.Initialize()
            cls.simulation.time = cls.simulation._GetSolver().AdvanceInTime(cls.simulation.time)
            cls.simulation.InitializeSolutionStep()

            # Assuming test_empirical_cubature_method.py passes, use pre-existing weights (w) and elements (z)
            # to ensure consistency and avoid duplication, as np.linalg.svd() can vary with thread settings (number of threads).
            cls.simulation._RomAnalysis__hrom_training_utility.hyper_reduction_element_selector.w = np.load("input_w.npy")
            cls.simulation._RomAnalysis__hrom_training_utility.hyper_reduction_element_selector.z = np.load("input_z.npy")
            cls.simulation._RomAnalysis__hrom_training_utility.AppendHRomWeightsToRomParameters()
            cls.simulation._RomAnalysis__hrom_training_utility.CreateHRomModelParts()

            # Load the expected and obtained model parts
            cls.expected_model_part = cls.model.CreateModelPart("expected")
            KratosMultiphysics.ModelPartIO("Expected_couette_flow_testHROM").ReadModelPart(cls.expected_model_part)

            cls.obtained_model_part = cls.model.CreateModelPart("obtained")
            KratosMultiphysics.ModelPartIO("../../../FluidDynamicsApplication/tests/CouetteFlowTest/couette_flow_testHROM").ReadModelPart(cls.obtained_model_part)

    def test_model_parts_and_numpy_arrays(self):
        # Testing Nodes in Main Model Part
        expected_nodes = sorted([node.Id for node in self.expected_model_part.Nodes])
        obtained_nodes = sorted([node.Id for node in self.obtained_model_part.Nodes])
        self.assertListEqual(expected_nodes, obtained_nodes, "Node IDs in the main model part do not match")

        for node_id in expected_nodes:
            expected_node = self.expected_model_part.Nodes[node_id]
            obtained_node = self.obtained_model_part.Nodes[node_id]
            self.assertAlmostEqual(expected_node.X, obtained_node.X, msg=f"X coordinate mismatch for node {node_id}")
            self.assertAlmostEqual(expected_node.Y, obtained_node.Y, msg=f"Y coordinate mismatch for node {node_id}")
            self.assertAlmostEqual(expected_node.Z, obtained_node.Z, msg=f"Z coordinate mismatch for node {node_id}")

        # Testing Elements in Main Model Part
        expected_elements = sorted([element.Id for element in self.expected_model_part.Elements])
        obtained_elements = sorted([element.Id for element in self.obtained_model_part.Elements])
        self.assertListEqual(expected_elements, obtained_elements, "Element IDs in the main model part do not match")

        # Testing Nodes in Submodel Part
        expected_sub_mp = self.expected_model_part.GetSubModelPart("FluidParts_Fluid")
        obtained_sub_mp = self.obtained_model_part.GetSubModelPart("FluidParts_Fluid")

        expected_nodes = sorted([node.Id for node in expected_sub_mp.Nodes])
        obtained_nodes = sorted([node.Id for node in obtained_sub_mp.Nodes])
        self.assertListEqual(expected_nodes, obtained_nodes, "Node IDs in 'FluidParts_Fluid' submodel part do not match")

        # Testing Elements in Submodel Part
        expected_elements = sorted([element.Id for element in expected_sub_mp.Elements])
        obtained_elements = sorted([element.Id for element in obtained_sub_mp.Elements])
        self.assertListEqual(expected_elements, obtained_elements, "Element IDs in 'FluidParts_Fluid' submodel part do not match")

        # Define file names
        file_names = [
            "ConditionIds", "ConditionWeights", "ElementIds", "ElementWeights"
        ]

        # Loop through file names and assert expected vs obtained
        for name in file_names:
            expected = np.load(f"{self.work_folder}rom_data/Expected_HROM_{name}.npy")
            obtained = np.load(f"{self.work_folder}rom_data/HROM_{name}.npy")

            # Check if both arrays are empty
            if expected.size == 0 and obtained.size == 0:
                self.assertTrue(True, f"Both expected and obtained {name} are correctly empty.")
            else:
                # Use Kratos' custom assert to compare arrays as matrices
                self.assertVectorAlmostEqual(expected, obtained, msg=f"Mismatch in {name}")


    def test_hrom_zero_weight_append_edge_cases(self):
        # Grab the utility from the simulation
        hrom_util = self.simulation._RomAnalysis__hrom_training_utility

        # Temporarily change the output folder to avoid polluting the main test's files
        from pathlib import Path
        edge_case_dir = Path(self.work_folder) / "rom_edge_case_data"
        edge_case_dir.mkdir(exist_ok=True)
        hrom_util.rom_basis_output_folder = edge_case_dir

        # Mock the ECM output: Select exactly two elements (indices 0 and 1) and NO conditions
        hrom_util.hyper_reduction_element_selector.z = np.array([0, 1])
        hrom_util.hyper_reduction_element_selector.w = np.array([1.5, 2.5])

        # Force the utility to append extra zero-weight entities
        hrom_util.include_minimum_condition = True
        hrom_util.include_condition_parents = True

        # Execute the method (this verifies the shape mismatch bug is fixed)
        hrom_util.AppendHRomWeightsToRomParameters()

        # Load the newly generated numpy files to assert their correctness
        cond_ids = np.load(edge_case_dir / "HROM_ConditionIds.npy")
        cond_weights = np.load(edge_case_dir / "HROM_ConditionWeights.npy")
        elem_ids = np.load(edge_case_dir / "HROM_ElementIds.npy")
        elem_weights = np.load(edge_case_dir / "HROM_ElementWeights.npy")

        # --- ASSERTIONS ---

        # A. Conditions edge case: Ensure conditions were safely added to the empty list
        self.assertGreater(len(cond_ids), 0, "Conditions should have been added by include_minimum_condition.")
        self.assertEqual(len(cond_ids), len(cond_weights), "Condition IDs and Weights shapes do not match.")
        for w in cond_weights:
            self.assertEqual(w, 0.0, "Appended conditions must have exactly 0.0 weight.")

        # B. Elements edge case: Ensure elements were appended correctly
        self.assertGreater(len(elem_ids), 2, "Parent elements should have been appended to the initial 2 elements.")
        self.assertEqual(len(elem_ids), len(elem_weights), "Element IDs and Weights shapes do not match.")

        # C. Verify the original non-zero elements retained their mapped IDs and weights
        elem_id_0 = hrom_util.numpy_index_to_element_id_mapping[0]
        elem_id_1 = hrom_util.numpy_index_to_element_id_mapping[1]

        idx_0 = np.where(elem_ids == elem_id_0)[0][0]
        idx_1 = np.where(elem_ids == elem_id_1)[0][0]

        self.assertEqual(elem_weights[idx_0], 1.5, "Original element weight was corrupted.")
        self.assertEqual(elem_weights[idx_1], 2.5, "Original element weight was corrupted.")

        # D. Ensure all other appended elements have exactly 0.0 weight
        for i, w in enumerate(elem_weights):
            if i not in [idx_0, idx_1]:
                self.assertEqual(w, 0.0, "Appended elements must have exactly 0.0 weight.")


    @classmethod
    def tearDownClass(cls):
        # Specific files and paths to delete
        specific_files = [
            "../../../FluidDynamicsApplication/tests/CouetteFlowTest/couette_flow_testHROM.mdpa",
            "rom_data/HROM_ConditionIds.npy",
            "rom_data/HROM_ConditionWeights.npy",
            "rom_data/HROM_ElementIds.npy",
            "rom_data/HROM_ElementWeights.npy",
            "rom_edge_case_data/HROM_ConditionIds.npy",
            "rom_edge_case_data/HROM_ConditionWeights.npy",
            "rom_edge_case_data/HROM_ElementIds.npy",
            "rom_edge_case_data/HROM_ElementWeights.npy"
        ]

        # Iterating over the specific file list and deleting each
        for file_path in specific_files:
            full_path = os.path.join(cls.work_folder, file_path)
            kratos_utilities.DeleteFileIfExisting(full_path)

        # Clean up the edge case directory itself
        edge_case_dir = os.path.join(cls.work_folder, "rom_edge_case_data")
        kratos_utilities.DeleteDirectoryIfExisting(edge_case_dir)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()