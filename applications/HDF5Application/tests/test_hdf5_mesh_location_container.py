from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

class TestMeshLocationContainer(KratosUnittest.TestCase):
    def test_SetHasGet(self):
        mesh_location = KratosHDF5.MeshLocationContainer()
        self.assertFalse(mesh_location.Has(-1, 1))
        mesh_location.Set(-1, 1, "test")
        self.assertTrue(mesh_location.Has(-1, 1))
        self.assertEqual(mesh_location.Get(-1, 1), "test")

        mesh_location.Set(-1, 1, "test1")
        self.assertTrue(mesh_location.Has(-1, 1))
        self.assertEqual(mesh_location.Get(-1, 1), "test1")

        mesh_location.Set(-1, 1, "test2")
        self.assertTrue(mesh_location.Has(-1, 1))
        self.assertEqual(mesh_location.Get(-1, 1), "test2")

        self.assertFalse(mesh_location.Has(-1, 2))
        mesh_location.Set(-1, 2, "test3")
        self.assertFalse(mesh_location.Has(1, 2))
        mesh_location.Set(1, 2, "test4")
        self.assertTrue(mesh_location.Has(-1, 1))
        self.assertEqual(mesh_location.Get(-1, 1), "test2")
        self.assertTrue(mesh_location.Has(1, 2))
        self.assertEqual(mesh_location.Get(1, 2), "test4")

    def test_AddHasGetParameters(self):
        parameters = Kratos.Parameters("""{}""")
        self.assertFalse(KratosHDF5.HasProcessId(parameters))
        KratosHDF5.AddProcessId(parameters, -1, 2)
        self.assertTrue(KratosHDF5.HasProcessId(parameters))
        self.assertTrue(parameters.IsEquivalentTo(Kratos.Parameters("""{"__hdf5_process_id": [-1, 2]}""")))
        self.assertEqual(KratosHDF5.GetProcessId(parameters), (-1, 2))

        with self.assertRaises(RuntimeError):
            KratosHDF5.AddProcessId(parameters, 1, 2)

    def test_ModelPart(self):
        mesh_location = KratosHDF5.MeshLocationContainer()
        mesh_location.Set(-1, 1, "test1")
        mesh_location.Set(-1, 2, "test2")
        mesh_location.Set(-1, 3, "test3")
        mesh_location.Set(-1, 4, "test4")
        mesh_location.Set(1, 1, "test11")
        mesh_location.Set(1, 2, "test12")
        mesh_location.Set(1, 3, "test13")
        mesh_location.Set(1, 4, "test14")

        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        self.assertFalse(KratosHDF5.HasMeshLocationContainer(model_part))
        KratosHDF5.SetMeshLocationContainer(model_part, mesh_location)
        self.assertTrue(KratosHDF5.HasMeshLocationContainer(model_part))
        self.assertEqual(KratosHDF5.GetMeshLocationContainer(model_part), mesh_location)

    def test_Serialization(self):
        self.addCleanup(kratos_utils.DeleteFileIfExisting, str((Path(__file__).parent / "test_serialization.rest").absolute()))

        mesh_location = KratosHDF5.MeshLocationContainer()
        mesh_location.Set(-1, 1, "test1")
        mesh_location.Set(-1, 2, "test2")
        mesh_location.Set(-1, 3, "test3")
        mesh_location.Set(-1, 4, "test4")
        mesh_location.Set(1, 1, "test11")
        mesh_location.Set(1, 2, "test12")
        mesh_location.Set(1, 3, "test13")
        mesh_location.Set(1, 4, "test14")

        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        KratosHDF5.SetMeshLocationContainer(model_part, mesh_location)

        file_serializer = Kratos.FileSerializer("test_serialization")
        file_serializer.Save("test_process_data", model)
        del file_serializer

        read_model = Kratos.Model()
        file_serializer = Kratos.FileSerializer("test_serialization")
        file_serializer.Load("test_process_data", read_model)

        self.assertTrue(read_model.HasModelPart("test"))
        read_model_part = read_model["test"]
        self.assertTrue(KratosHDF5.HasMeshLocationContainer(read_model_part))

        read_mesh_locations = KratosHDF5.GetMeshLocationContainer(read_model_part)
        self.assertNotEqual(mesh_location, read_mesh_locations)
        self.assertEqual(str(mesh_location), str(read_mesh_locations))

if __name__ == "__main__":
    KratosUnittest.main()