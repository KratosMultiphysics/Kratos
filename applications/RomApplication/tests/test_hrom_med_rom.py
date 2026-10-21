import os
import shutil
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication as KratosROM

if kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication", "MedApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication
    import KratosMultiphysics.MedApplication

# Get the absolute path to the MED test file
KRATOS_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
MED_TEST_FILE = os.path.join(KRATOS_ROOT, "applications", "MedApplication", "tests", "med_files", "tetrahedral_10N.med")

@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "MedApplication")
class TestHromMedRom(KratosUnittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.work_folder = "hrom_med_test_files/"

    def test_hrom_create_entities_from_geometries_modeler(self):
        """Test the HRomCreateEntitiesFromGeometriesModeler directly."""
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory

            # Create a model with geometries (simulating MED import)
            model = KratosMultiphysics.Model()
            model_part = model.CreateModelPart("TestPart")
            io = KratosMultiphysics.MedApplication.MedModelPartIO(MED_TEST_FILE, KratosMultiphysics.IO.READ)
            io.ReadModelPart(model_part)

            # Create entities from geometries using the standard modeler
            modelers = KratosModelParametersFactory(model).ConstructListOfItems(
                KratosMultiphysics.Parameters("""
                [{
                    "name" : "Modelers.KratosMultiphysics.CreateEntitiesFromGeometriesModeler",
                    "parameters" : {
                        "elements_list" : [{
                            "model_part_name" : "TestPart",
                            "element_name" : "SmallDisplacementElement3D10N"
                        }],
                        "conditions_list" : [{
                            "model_part_name" : "TestPart.load",
                            "condition_name" : "SmallDisplacementSurfaceLoadCondition3D6N"
                        }]
                    }
                }]
                """)
            )

            for modeler in modelers:
                modeler.SetupGeometryModel()
            for modeler in modelers:
                modeler.PrepareGeometryModel()
            for modeler in modelers:
                modeler.SetupModelPart()

            # Check entities were created
            self.assertGreater(model_part.NumberOfElements(), 0)
            self.assertGreater(model_part.NumberOfConditions(), 0)

            # Check AssignPropertiesToEntitiesIfMissing
            KratosROM.RomAuxiliaryUtilities.AssignPropertiesToEntitiesIfMissing(model_part)

            for element in model_part.Elements:
                self.assertIsNotNone(element.Properties)

            for condition in model_part.Conditions:
                self.assertIsNotNone(condition.Properties)

    def test_hrom_create_entities_from_geometries_modeler_with_id_assignment(self):
        """Test the HRomCreateEntitiesFromGeometriesModeler with ID assignment from file."""
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory

            # Create a model with geometries (simulating MED import)
            model = KratosMultiphysics.Model()
            model_part = model.CreateModelPart("TestPart")
            io = KratosMultiphysics.MedApplication.MedModelPartIO(MED_TEST_FILE, KratosMultiphysics.IO.READ)
            io.ReadModelPart(model_part)

            # Create a geometry IDs file (mapping geometry IDs to entity IDs) in the work folder
            geometry_ids_file = "geometry_ids_test.csv"
            with open(geometry_ids_file, "w") as f:
                for element in model_part.Elements:
                    f.write(f"{element.GetGeometry().Id},{element.Id}\n")
                for condition in model_part.Conditions:
                    f.write(f"{condition.GetGeometry().Id},{condition.Id}\n")

# Create entities from geometries using the HRom modeler with ID assignment
            modelers = KratosModelParametersFactory(model).ConstructListOfItems(
                KratosMultiphysics.Parameters("""
                [{
                    "name" : "Modelers.KratosMultiphysics.RomApplication.HRomCreateEntitiesFromGeometriesModeler",
                    "parameters" : {
                        "elements_list" : [{
                            "model_part_name" : "TestPart",
                            "element_name" : "SmallDisplacementElement3D10N"
                        }],
                        "conditions_list" : [{
                            "model_part_name" : "TestPart.load",
                            "condition_name" : "SmallDisplacementSurfaceLoadCondition3D6N"
                        }],
                        "assign_ids_from_file" : true,
                        "geometry_ids_file" : "%s"
                    }
                }]
                """ % geometry_ids_file)
            )

            for modeler in modelers:
                modeler.SetupGeometryModel()
            for modeler in modelers:
                modeler.PrepareGeometryModel()
            for modeler in modelers:
                modeler.SetupModelPart()

            # Check entities were created
            self.assertGreater(model_part.NumberOfElements(), 0)
            self.assertGreater(model_part.NumberOfConditions(), 0)

            # Check AssignPropertiesToEntitiesIfMissing
            KratosROM.RomAuxiliaryUtilities.AssignPropertiesToEntitiesIfMissing(model_part)

            for element in model_part.Elements:
                self.assertIsNotNone(element.Properties)

            for condition in model_part.Conditions:
                self.assertIsNotNone(condition.Properties)

            # Clean up
            kratos_utilities.DeleteFileIfExisting(geometry_ids_file)

    def test_add_geometries_utility(self):
        """Test the AddGeometries utility function."""
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory

            # Create origin model part with geometries and entities
            model = KratosMultiphysics.Model()
            origin_mp = model.CreateModelPart("OriginPart")
            io = KratosMultiphysics.MedApplication.MedModelPartIO(MED_TEST_FILE, KratosMultiphysics.IO.READ)
            io.ReadModelPart(origin_mp)

            modelers = KratosModelParametersFactory(model).ConstructListOfItems(
                KratosMultiphysics.Parameters("""
                [{
                    "name" : "Modelers.KratosMultiphysics.CreateEntitiesFromGeometriesModeler",
                    "parameters" : {
                        "elements_list" : [{
                            "model_part_name" : "OriginPart",
                            "element_name" : "SmallDisplacementElement3D10N"
                        }],
                        "conditions_list" : [{
                            "model_part_name" : "OriginPart.load",
                            "condition_name" : "SmallDisplacementSurfaceLoadCondition3D6N"
                        }]
                    }
                }]
                """)
            )

            for modeler in modelers:
                modeler.SetupGeometryModel()
            for modeler in modelers:
                modeler.PrepareGeometryModel()
            for modeler in modelers:
                modeler.SetupModelPart()

            # Create destination model part (simulating HROM model part) with only some entities
            dest_model = KratosMultiphysics.Model()
            dest_mp = dest_model.CreateModelPart("DestPart")

            # Copy some nodes from origin
            for node in origin_mp.Nodes:
                dest_mp.CreateNewNode(node.Id, node.X, node.Y, node.Z)

            # Copy a subset of elements and conditions (simulating HROM selection)
            element_ids = [e.Id for e in origin_mp.Elements][:10]  # First 10 elements
            condition_ids = [c.Id for c in origin_mp.Conditions][:5]  # First 5 conditions

            for elem_id in element_ids:
                orig_elem = origin_mp.GetElement(elem_id)
                # Get the node IDs of the geometry
                node_ids = [n.Id for n in orig_elem.GetGeometry()]
                dest_mp.CreateNewElement("SmallDisplacementElement3D10N", elem_id, node_ids, orig_elem.Properties)

            for cond_id in condition_ids:
                orig_cond = origin_mp.GetCondition(cond_id)
                node_ids = [n.Id for n in orig_cond.GetGeometry()]
                dest_mp.CreateNewCondition("SmallDisplacementSurfaceLoadCondition3D6N", cond_id, node_ids, orig_cond.Properties)

            # Add geometries from origin to destination
            KratosROM.RomAuxiliaryUtilities.AddGeometries(origin_mp, dest_mp)

            # Check that geometries were added
            self.assertGreater(dest_mp.NumberOfGeometries(), 0)

            # Check that the geometries match the entities
            for element in dest_mp.Elements:
                self.assertIsNotNone(element.GetGeometry())

            for condition in dest_mp.Conditions:
                self.assertIsNotNone(condition.GetGeometry())

    def test_assign_properties_to_entities_if_missing(self):
        """Test the AssignPropertiesToEntitiesIfMissing utility function."""
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory

            # Create a model with geometries (simulating MED import)
            model = KratosMultiphysics.Model()
            model_part = model.CreateModelPart("TestPart")
            io = KratosMultiphysics.MedApplication.MedModelPartIO(MED_TEST_FILE, KratosMultiphysics.IO.READ)
            io.ReadModelPart(model_part)

            # Create entities from geometries using the standard modeler
            modelers = KratosModelParametersFactory(model).ConstructListOfItems(
                KratosMultiphysics.Parameters("""
                [{
                    "name" : "Modelers.KratosMultiphysics.CreateEntitiesFromGeometriesModeler",
                    "parameters" : {
                        "elements_list" : [{
                            "model_part_name" : "TestPart",
                            "element_name" : "SmallDisplacementElement3D10N"
                        }],
                        "conditions_list" : [{
                            "model_part_name" : "TestPart.load",
                            "condition_name" : "SmallDisplacementSurfaceLoadCondition3D6N"
                        }]
                    }
                }]
                """)
            )

            for modeler in modelers:
                modeler.SetupGeometryModel()
            for modeler in modelers:
                modeler.PrepareGeometryModel()
            for modeler in modelers:
                modeler.SetupModelPart()

            # Entities created from geometries have no properties
            for element in model_part.Elements:
                self.assertIsNone(element.Properties)

            for condition in model_part.Conditions:
                self.assertIsNone(condition.Properties)

            # Assign properties
            KratosROM.RomAuxiliaryUtilities.AssignPropertiesToEntitiesIfMissing(model_part)

            # Check that all entities now have properties
            for element in model_part.Elements:
                self.assertIsNotNone(element.Properties)

            for condition in model_part.Conditions:
                self.assertIsNotNone(condition.Properties)

            # Check that submodelparts also get properties
            for smp in model_part.SubModelParts:
                for element in smp.Elements:
                    self.assertIsNotNone(element.Properties)
                for condition in smp.Conditions:
                    self.assertIsNotNone(condition.Properties)

    def test_run_hrom_on_hrom_model_part(self):
        """Trains a ROM and an HROM with the RomManager on the MED case and checks that the HROM run on the HROM model part
        (written during the HROM training) gives the solution of the HROM run on the full model part."""
        self._RunHromOnHromModelPart("RomManager")

    def test_run_hrom_on_hrom_model_part_without_nodes_in_a_sub_model_part(self):
        """As test_run_hrom_on_hrom_model_part, with the HROM elements selected among the neighbours of the load, so that
        the HROM model part has no node of the fixed sub model part ("fix"), which is used by a boundary condition process."""
        self._RunHromOnHromModelPart("RomManagerWithoutFixNodes", initial_candidate_elements_model_part_list=["Structure.load"])

    def _RunHromOnHromModelPart(self, folder_name, initial_candidate_elements_model_part_list=None):
        work_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), self.work_folder, folder_name)
        self.addCleanup(kratos_utilities.DeleteDirectoryIfExisting, work_folder)
        os.makedirs(work_folder, exist_ok=True)
        # The HROM model part is written next to the input MED file, so the case is run on a copy of it
        for file_name in ["ProjectParameters.json", "Materials.json", "rom_manager_parameters.json"]:
            shutil.copy(os.path.join(os.path.dirname(work_folder), file_name), work_folder)
        shutil.copy(MED_TEST_FILE, work_folder)

        with KratosUnittest.WorkFolderScope(work_folder, __file__):
            from KratosMultiphysics.RomApplication.rom_manager import RomManager

            displacements = {}
            def CustomizeSimulation(cls, global_model, parameters, mu=None):
                class CustomSimulation(cls):
                    def Finalize(self):
                        model_part = self._GetSolver().GetComputingModelPart().GetRootModelPart()
                        displacements[CustomizeSimulation.label] = {node.Id: list(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)) for node in model_part.Nodes}
                        super().Finalize()
                return CustomSimulation(global_model, parameters)

            with open("rom_manager_parameters.json", 'r') as parameter_file:
                general_rom_manager_parameters = KratosMultiphysics.Parameters(parameter_file.read())
            if initial_candidate_elements_model_part_list is not None:
                general_rom_manager_parameters["HROM"]["initial_candidate_elements_model_part_list"].SetStringArray(initial_candidate_elements_model_part_list)
            rom_manager = RomManager(project_parameters_name="ProjectParameters.json", general_rom_manager_parameters=general_rom_manager_parameters, CustomizeSimulation=CustomizeSimulation)
            mu = [[1.0]]
            CustomizeSimulation.label = "Fit"
            rom_manager.Fit(mu)
            self.assertTrue(os.path.isfile("tetrahedral_10NHROM.med"))
            self.assertTrue(os.path.isfile("rom_data/HROM_GeometryIds.csv"))

            for label, use_full_model_part in [("full model part", True), ("HROM model part", False)]:
                CustomizeSimulation.label = label
                rom_manager.RunHROM(mu_run=mu, mu_train=mu, use_full_model_part=use_full_model_part)

            # The HROM model part has only some of the nodes, on which both solutions must coincide
            full, hrom = displacements["full model part"], displacements["HROM model part"]
            self.assertLess(len(hrom), len(full))
            expected = np.array([full[node_id] for node_id in hrom])
            obtained = np.array([hrom[node_id] for node_id in hrom])
            self.assertLess(np.linalg.norm(obtained - expected) / np.linalg.norm(expected), 1.0e-10)

    @classmethod
    def tearDownClass(cls):
        pass

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()