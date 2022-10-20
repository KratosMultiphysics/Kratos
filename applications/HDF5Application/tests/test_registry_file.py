# --- Kratos Imports ---
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.HDF5Application as HDF5

# --- STD Imports ---
import pathlib


class TestRegistryFile(KratosUnittest.TestCase):

    @property
    def test_file_path(self) -> pathlib.Path:
        return pathlib.Path("test_registry_file.log")

    def test_RegistryFile(self) -> None:
        registry = HDF5.RegistryFile(self.test_file_path)
        model = KratosMultiphysics.Model()

        extractor = lambda model: KratosMultiphysics.Parameters('["1st","2nd"]')
        registry.SetExtractor(extractor)

        registry.Clear()
        self.assertEqual(len(registry), 0)

        for i in range(1, 3):
            registry.Push(model)
            self.assertEqual(len(registry), i)
            for item in registry:
                self.assertTrue(isinstance(item, KratosMultiphysics.Parameters))
                self.assertEqual(item.WriteJsonString(), '["1st","2nd"]')

        registry.Clear()
        self.assertEqual(len(registry), 0)

    def test_Extractor(self) -> None:
        def extractor(model: KratosMultiphysics.Model) -> KratosMultiphysics.Parameters:
            model_part = model.GetModelPart("test")
            process_info = model_part.ProcessInfo
            return KratosMultiphysics.Parameters(f"""[
                "example/relative/directory/example_file_name_{process_info[KratosMultiphysics.STEP]}.h5",
                {process_info[KratosMultiphysics.STEP]},
                {process_info[KratosMultiphysics.TIME]},
                true
            ]""")

        registry = HDF5.RegistryFile(self.test_file_path, extractor)
        registry.Clear()

        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("test")

        for i_step in range(10):
            model_part.CloneTimeStep()
            model_part.ProcessInfo[KratosMultiphysics.STEP] = i_step
            model_part.ProcessInfo[KratosMultiphysics.TIME] = i_step / 10.0
            registry.Push(model)

        self.assertEqual(len(registry), 10)
        for index, item in enumerate(registry):
            self.assertTrue(isinstance(item, KratosMultiphysics.Parameters))
            self.assertTrue(item.IsArray())
            self.assertEqual(item.size(), 4)

            self.assertTrue(item[0].IsString())
            self.assertEqual(item[0].GetString(), f"example/relative/directory/example_file_name_{index}.h5")

            self.assertTrue(item[1].IsInt())
            self.assertEqual(item[1].GetInt(), index)

            self.assertTrue(item[2].IsDouble())
            self.assertEqual(item[2].GetDouble(), index / 10.0)

            self.assertTrue(item[3].IsBool())
            self.assertEqual(item[3].GetBool(), True)


if __name__ == "__main__":
    KratosUnittest.main()
