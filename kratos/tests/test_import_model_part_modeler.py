import os

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.modelers.import_mdpa_modeler import ImportMDPAModeler

def GetFilePath(filename):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), filename)

class TestImportMDPAModeler(KratosUnittest.TestCase):
    @classmethod
    def tearDown(cls):
        # Clean up temporary files
        KratosUtilities.DeleteFileIfExisting(GetFilePath("test_model_part_io_write.out.mdpa"))
        KratosUtilities.DeleteFileIfExisting(GetFilePath("test_model_part_io_write.out.time"))
        KratosUtilities.DeleteFileIfExisting(GetFilePath("test_model_part_io_write.time"))
        KratosUtilities.DeleteFileIfExisting(GetFilePath("test_model_part_io_write_mesh_only.out.mdpa"))
        KratosUtilities.DeleteFileIfExisting(GetFilePath("test_model_part_io_write_mesh_only.out.time"))
        KratosUtilities.DeleteFileIfExisting(GetFilePath("test_model_part_io_write_mesh_only.time"))

    def testImportMDPAModeler(self):
        # Set up the import model part modeler
        model = KratosMultiphysics.Model()
        settings = KratosMultiphysics.Parameters('''{
            "input_filename" : "",
            "model_part_name" : "Main"
        }''')
        settings["input_filename"].SetString(GetFilePath("test_files/mdpa_files/test_model_part_io_read"))
        import_mdpa_modeler = ImportMDPAModeler(model, settings)

        # Get the model part created by the modeler
        model_part = model.GetModelPart(settings["model_part_name"].GetString())

        # Add nodal solution step variable data
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # Call the modeler methods
        import_mdpa_modeler.SetupGeometryModel()
        import_mdpa_modeler.PrepareGeometryModel()
        import_mdpa_modeler.SetupModelPart()

        # Check results
        self.assertEqual(model_part.NumberOfSubModelParts(), 2)

        self.assertEqual(model_part.NumberOfTables(), 1)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfNodes(), 6)
        self.assertEqual(model_part.NumberOfProperties(), 1)
        self.assertEqual(model_part.NumberOfGeometries(), 9)
        self.assertEqual(model_part.NumberOfElements(), 4)
        self.assertEqual(model_part.NumberOfConditions(), 5)

        self.assertEqual(model_part[KratosMultiphysics.AMBIENT_TEMPERATURE], 250.0)
        self.assertEqual(model_part[KratosMultiphysics.DISPLACEMENT_X], 2.1)
        self.assertEqual(model_part[KratosMultiphysics.DISPLACEMENT_Y], 3.2)
        self.assertEqual(model_part[KratosMultiphysics.DISPLACEMENT_Z], 4.3)
        self.assertEqual(model_part[KratosMultiphysics.VELOCITY_X], 3.8)
        self.assertEqual(model_part[KratosMultiphysics.VELOCITY_Y], 4.9)
        self.assertEqual(model_part[KratosMultiphysics.VELOCITY_Z], 0.0)

        self.assertTrue(model_part.GetNode(1).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertTrue(model_part.GetNode(2).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertFalse(model_part.GetNode(3).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertFalse(model_part.GetNode(972).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertTrue(model_part.GetNode(973).IsFixed(KratosMultiphysics.DISPLACEMENT_X))
        self.assertTrue(model_part.GetNode(974).IsFixed(KratosMultiphysics.DISPLACEMENT_X))

        self.assertTrue(model_part.GetNode(1).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertTrue(model_part.GetNode(2).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertFalse(model_part.GetNode(3).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertFalse(model_part.GetNode(972).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertTrue(model_part.GetNode(973).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))
        self.assertTrue(model_part.GetNode(974).IsFixed(KratosMultiphysics.DISPLACEMENT_Y))

        self.assertTrue(model_part.GetNode(1).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertTrue(model_part.GetNode(2).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertFalse(model_part.GetNode(3).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertFalse(model_part.GetNode(972).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertTrue(model_part.GetNode(973).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))
        self.assertTrue(model_part.GetNode(974).IsFixed(KratosMultiphysics.DISPLACEMENT_Z))

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.1)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.2)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X),0.0)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.0)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.0)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.000973)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y), 0.000974)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z), 0.0)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.01)
        self.assertEqual(model_part.GetNode(2).GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.01)
        self.assertEqual(model_part.GetNode(3).GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.0)
        self.assertEqual(model_part.GetNode(972).GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.0)
        self.assertEqual(model_part.GetNode(973).GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.01)
        self.assertEqual(model_part.GetNode(974).GetSolutionStepValue(KratosMultiphysics.VISCOSITY), 0.01)

        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_X), 1.1)
        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y), 2.2)
        self.assertEqual(model_part.GetNode(1).GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z), 3.3)

        self.assertTrue(model_part.HasSubModelPart("Inlets"))

        inlets_model_part = model_part.GetSubModelPart("Inlets")

        self.assertEqual(inlets_model_part.NumberOfTables(), 1)
        self.assertEqual(inlets_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlets_model_part.NumberOfNodes(), 3)
        self.assertEqual(inlets_model_part.NumberOfGeometries(), 0)
        self.assertEqual(inlets_model_part.NumberOfElements(), 1)
        self.assertEqual(inlets_model_part.NumberOfConditions(), 3)
        self.assertEqual(inlets_model_part.NumberOfSubModelParts(), 2)
        self.assertTrue(inlets_model_part.HasSubModelPart("Inlet1"))
        self.assertTrue(inlets_model_part.HasSubModelPart("Inlet2"))

        inlet1_model_part = inlets_model_part.GetSubModelPart("Inlet1")

        self.assertEqual(inlet1_model_part.NumberOfTables(), 0)
        self.assertEqual(inlet1_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlet1_model_part.NumberOfNodes(), 2)
        self.assertEqual(inlet1_model_part.NumberOfGeometries(), 0)
        self.assertEqual(inlet1_model_part.NumberOfElements(), 0)
        self.assertEqual(inlet1_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlet1_model_part.NumberOfSubModelParts(), 0)

        inlet2_model_part = inlets_model_part.GetSubModelPart("Inlet2")

        self.assertEqual(inlet2_model_part.NumberOfTables(), 0)
        self.assertEqual(inlet2_model_part.NumberOfProperties(), 0)
        self.assertEqual(inlet2_model_part.NumberOfNodes(), 0)
        self.assertEqual(inlet2_model_part.NumberOfGeometries(), 0)
        self.assertEqual(inlet2_model_part.NumberOfElements(), 0)
        self.assertEqual(inlet2_model_part.NumberOfConditions(), 2)
        self.assertEqual(inlet2_model_part.NumberOfSubModelParts(), 0)

        self.assertTrue(model_part.HasSubModelPart("Outlet"))

        outlet_model_part = model_part.GetSubModelPart("Outlet")

        self.assertEqual(outlet_model_part.NumberOfTables(), 0)
        self.assertEqual(outlet_model_part.NumberOfProperties(), 1)
        self.assertEqual(outlet_model_part.NumberOfNodes(), 0)
        self.assertEqual(outlet_model_part.NumberOfGeometries(), 0)
        self.assertEqual(outlet_model_part.NumberOfElements(), 0)
        self.assertEqual(outlet_model_part.NumberOfConditions(), 1)
        self.assertEqual(outlet_model_part.NumberOfSubModelParts(), 0)

        properties_1 = model_part.GetProperties()[1]
        #Bools
        self.assertTrue(properties_1[KratosMultiphysics.IS_RESTARTED])
        self.assertFalse(properties_1[KratosMultiphysics.COMPUTE_DYNAMIC_TANGENT])
        self.assertFalse(properties_1[KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX])
        #Double
        self.assertEqual(properties_1[KratosMultiphysics.DENSITY], 3.4E-5)
        #Array3
        self.assertEqual(properties_1[KratosMultiphysics.VOLUME_ACCELERATION][0], 0.00)
        self.assertEqual(properties_1[KratosMultiphysics.VOLUME_ACCELERATION][1], 0.00)
        self.assertEqual(properties_1[KratosMultiphysics.VOLUME_ACCELERATION][2], 9.8)
        #Matrix3x3
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][0,0], 0)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][0,1], 0.27)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][0,2], 0.27)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][1,0], 0.087)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][1,1], 0)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][1,2], 0.27)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][2,0], 0.075)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][2,1], 0.23)
        self.assertEqual(properties_1[KratosMultiphysics.LOCAL_INERTIA_TENSOR][2,2], 0)

        #SubModelPartData
        self.assertTrue(inlets_model_part[KratosMultiphysics.IS_RESTARTED])
        self.assertTrue(inlets_model_part[KratosMultiphysics.COMPUTE_LUMPED_MASS_MATRIX])
        self.assertFalse(inlets_model_part[KratosMultiphysics.COMPUTE_DYNAMIC_TANGENT])

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
