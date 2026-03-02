import os

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils

dependencies_are_available = KratosUtils.CheckIfApplicationsAvailable("FluidDynamicsApplication")
if dependencies_are_available:
    from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis


def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestSerializer(KratosUnittest.TestCase):
    def _prepare_fluid_test(self):
        # Define a model and load the parameters
        self.pre_serialized_model = KratosMultiphysics.Model()
        with open(GetFilePath("auxiliar_files_for_python_unittest/parameters_files/test_serializer.json"),'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        file_name = parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(GetFilePath(file_name))
        # First the model is initialized
        self.pre_serialized_simulation = FluidDynamicsAnalysis(self.pre_serialized_model, parameters)
        self.pre_serialized_simulation.Initialize()

        # Before serializing the model, main model part is set to RESTARTED
        self.main_model_part_name = parameters["solver_settings"]["model_part_name"].GetString()
        self.pre_serialized_model.GetModelPart(self.main_model_part_name).ProcessInfo.SetValue(KratosMultiphysics.IS_RESTARTED,True)
        self.serialized_model = KratosMultiphysics.StreamSerializer()
        self.serialized_model.Save("ModelSerialization",self.pre_serialized_model)

        with open(GetFilePath("auxiliar_files_for_python_unittest/parameters_files/test_serializer.json"),'r') as parameter_file:
            self.project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
        # Parameters are read again and input type set to use_input_model_part since the serialized model already has the mdpa loaded
        self.project_parameters["solver_settings"]["model_import_settings"]["input_type"].SetString("use_input_model_part")

        # Deserialize and store the new model
        self.current_model = KratosMultiphysics.Model()
        self.serialized_model.Load("ModelSerialization",self.current_model)

    def _check_results(self):
        pre_serialized_model_part = self.pre_serialized_model.GetModelPart(self.main_model_part_name)
        pre_serialized_pressure_results = [node.GetSolutionStepValue(KratosMultiphysics.PRESSURE) for node in pre_serialized_model_part.Nodes]
        pre_serialized_velocity_results = [node.GetSolutionStepValue(KratosMultiphysics.VELOCITY) for node in pre_serialized_model_part.Nodes]

        serialized_model_part = self.current_model.GetModelPart(self.main_model_part_name)
        serialized_pressure_results = [node.GetSolutionStepValue(KratosMultiphysics.PRESSURE) for node in serialized_model_part.Nodes]
        serialized_velocity_results = [node.GetSolutionStepValue(KratosMultiphysics.VELOCITY) for node in serialized_model_part.Nodes]

        # Comparing results before and after serializing
        for pre_serialized_result, serialized_result in zip(pre_serialized_pressure_results,serialized_pressure_results):
            self.assertAlmostEqual(pre_serialized_result, serialized_result)
        for pre_serialized_result, serialized_result in zip(pre_serialized_velocity_results,serialized_velocity_results):
            for value_pre_seralized, value_serialized in zip(pre_serialized_result, serialized_result):
                self.assertAlmostEqual(value_pre_seralized, value_serialized)

    @KratosUnittest.skipUnless(dependencies_are_available,"FluidDynamicsApplication is not available")
    def test_serializer_fluid_analysis(self):
        self._prepare_fluid_test()
        # Solving simulation before serializing to later check the results
        self.pre_serialized_simulation.RunSolutionLoop()
        self.pre_serialized_simulation.Finalize()
        # Solving simulation after serializing
        self.serialized_simulation = FluidDynamicsAnalysis(self.current_model, self.project_parameters)
        self.serialized_simulation.Run()
        self._check_results()

    def test_serializer_loading(self):
        # Creating a model with nodes and variables
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        model_part.CreateSubModelPart("Inlets")
        model_part.CreateSubModelPart("Temp")
        model_part.CreateNewNode(1,0.0,0.0,0.0)
        other = current_model.CreateModelPart("Other")
        other.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        other.CreateNewNode(1,0.0,0.0,0.0)

        # Serializing model
        serialized_model = KratosMultiphysics.StreamSerializer()
        serialized_model.Save("ModelSerialization", current_model)

        # Loading model several times
        first_model = KratosMultiphysics.Model()
        serialized_model.Load("ModelSerialization", first_model)
        second_model = KratosMultiphysics.Model()
        serialized_model.LoadFromBeginning("ModelSerialization", second_model)
        third_model = KratosMultiphysics.Model()
        serialized_model.LoadFromBeginning("ModelSerialization", third_model)

        self.assertTrue(third_model["Main"].HasNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE))
        self.assertTrue(third_model["Other"].HasNodalSolutionStepVariable(KratosMultiphysics.PRESSURE))
        self.assertTrue(third_model.HasModelPart("Main.Inlets"))
        self.assertTrue(third_model.HasModelPart("Main.Temp"))
        self.assertTrue(1 in third_model["Main"].Nodes)
        self.assertTrue(1 in third_model["Other"].Nodes)

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def test_DataOnly(self):
        self.addCleanup(KratosUtils.DeleteFileIfExisting, "serializer_data_only.rest")
        import KratosMultiphysics.StructuralMechanicsApplication

        KratosUtils.DeleteFileIfExisting("serializer_data_only.rest")

        def create_model_part(factor = 1) -> KratosMultiphysics.Model:
            model = KratosMultiphysics.Model()
            model_part = model.CreateModelPart("test")
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.RECOVERED_STRESS)
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE)

            model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
            model_part.ProcessInfo[KratosMultiphysics.PRESSURE] = 4.1 * factor
            model_part.ProcessInfo[KratosMultiphysics.VELOCITY] = KratosMultiphysics.Array3([1.6 * factor, 2.4 * factor, 6.1 * factor])
            model_part.ProcessInfo[KratosMultiphysics.RECOVERED_STRESS] = KratosMultiphysics.Vector([7.1 * factor, 1.3 * factor, 2.0 * factor, 2.9 * factor, 2.7 * factor, 6.1 * factor])
            model_part.ProcessInfo[KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE] = KratosMultiphysics.Matrix([[4.2 * factor, 6.1 * factor, 8.3 * factor], [7.4 * factor, 8.0 * factor, 9.9 * factor]])

            # create nodes
            for i in range(5):
                node: KratosMultiphysics.Node = model_part.CreateNewNode(i + 1, i + 1, 0, 0)
                node.SetValue(KratosMultiphysics.PRESSURE, (10.0 + i) * factor)
                node.SetValue(KratosMultiphysics.VELOCITY, KratosMultiphysics.Array3([(1.1 + i) * factor, (2.1 + i) * factor, (3.1 + i) * factor]))
                node.SetValue(KratosMultiphysics.RECOVERED_STRESS, KratosMultiphysics.Vector([(1.1 + i) * factor, (1.2 + i) * factor, (2.3 + i) * factor, (2.1 + i) * factor, (2.4 + i) * factor, (6.4 + i) * factor]))
                node.SetValue(KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE, KratosMultiphysics.Matrix([[(1.1 + i) * factor, (1.2 + i) * factor, (2.3 + i) * factor], [(2.1 + i) * factor, (2.4 + i) * factor, (6.4 + i) * factor]]))

                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, (10.0 + (i + 1)) * 3.1 * factor)
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, KratosMultiphysics.Array3([(1.1 + (i + 1)) * 3.1 * factor, (2.1 + (i + 1)) * 3.1 * factor, (3.1 + (i + 1)) * 3.1 * factor]))
                node.SetSolutionStepValue(KratosMultiphysics.RECOVERED_STRESS, KratosMultiphysics.Vector([(1.1 + (i + 1)) * 3.1 * factor, (1.2 + (i + 1)) * 3.1 * factor, (2.3 + (i + 1)) * 3.1 * factor, (2.1 + (i + 1)) * 3.1 * factor, (2.4 + (i + 1)) * 3.1 * factor, (6.4 + (i + 1)) * 3.1 * factor]))
                node.SetSolutionStepValue(KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE, KratosMultiphysics.Matrix([[(1.1 + (i + 1)) * 3.1 * factor, (1.2 + (i + 1)) * 3.1 * factor, (2.3 + (i + 1)) * 3.1 * factor], [(2.1 + (i + 1)) * 3.1 * factor, (2.4 + (i + 1)) * 3.1 * factor, (6.4 + (i + 1)) * 3.1 * factor]]))

            sub_mp_1 = model_part.CreateSubModelPart("sub_mp_1.sub_mp2")
            sub_mp_1.AddNodes([1, 2])

            material_parameters = KratosMultiphysics.Parameters("""
            {
                "properties" : [{
                    "model_part_name" : "test",
                    "properties_id"   : 1,
                    "Material"        : {
                        "constitutive_law" : {
                            "name" : "LinearElasticPlaneStress2DLaw"
                        },
                        "Variables"        : {
                            "DENSITY"                : 1.4,
                            "YOUNG_MODULUS"          : 2.4,
                            "VELOCITY"               : [1.2, 3.4, 56.9],
                            "RECOVERED_STRESS"       : [5.6, 10.3, 11.4, 67.8]
                        },
                        "Tables"           : {}
                    }
                }]
            }
            """)
            KratosMultiphysics.ReadMaterialsUtility(model).ReadMaterials(material_parameters)

            properties = model_part.GetProperties(1)

            # create beam element
            for i in range(model_part.NumberOfNodes() - 1):
                element: KratosMultiphysics.Element = model_part.CreateNewElement("CrBeamElement3D2N", i + 1, [i + 1, i + 2], properties)
                element.SetValue(KratosMultiphysics.PRESSURE, (10.0 + i) * factor * 2)
                element.SetValue(KratosMultiphysics.VELOCITY, KratosMultiphysics.Array3([(1.1 + i) * factor * 2, (2.1 + i) * factor * 2, (3.1 + i) * factor * 2]))
                element.SetValue(KratosMultiphysics.RECOVERED_STRESS, KratosMultiphysics.Vector([(1.1 + i) * factor * 2, (1.2 + i) * factor * 2, (2.3 + i) * factor * 2, (2.1 + i) * factor * 2, (2.4 + i) * factor * 2, (6.4 + i) * factor * 2]))
                element.SetValue(KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE, KratosMultiphysics.Matrix([[(1.1 + i) * factor * 2, (1.2 + i) * factor * 2, (2.3 + i) * factor * 2], [(2.1 + i) * factor * 2, (2.4 + i) * factor * 2, (6.4 + i) * factor * 2]]))

            return model

        def comparator(model_1: KratosMultiphysics.Model, model_2: KratosMultiphysics.Model, factor: int):
            mp1 = model_1["test"]
            mp2 = model_2["test"]

            # process info check
            self.assertEqual(mp1.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE], mp2.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

            self.assertAlmostEqual(mp1.ProcessInfo[KratosMultiphysics.PRESSURE], mp2.ProcessInfo[KratosMultiphysics.PRESSURE] * factor)
            self.assertVectorAlmostEqual(mp1.ProcessInfo[KratosMultiphysics.VELOCITY], mp2.ProcessInfo[KratosMultiphysics.VELOCITY] * factor)
            self.assertVectorAlmostEqual(mp1.ProcessInfo[KratosMultiphysics.RECOVERED_STRESS], mp2.ProcessInfo[KratosMultiphysics.RECOVERED_STRESS] * factor)
            self.assertMatrixAlmostEqual(mp1.ProcessInfo[KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE], mp2.ProcessInfo[KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE] * factor)

            for node_1, node_2 in zip(mp1.Nodes, mp2.Nodes):
                self.assertAlmostEqual(node_1.GetValue(KratosMultiphysics.PRESSURE), node_2.GetValue(KratosMultiphysics.PRESSURE) * factor)
                self.assertVectorAlmostEqual(node_1.GetValue(KratosMultiphysics.VELOCITY), node_2.GetValue(KratosMultiphysics.VELOCITY) * factor)
                self.assertVectorAlmostEqual(node_1.GetValue(KratosMultiphysics.RECOVERED_STRESS), node_2.GetValue(KratosMultiphysics.RECOVERED_STRESS) * factor)
                self.assertMatrixAlmostEqual(node_1.GetValue(KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE), node_2.GetValue(KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE) * factor)

                self.assertAlmostEqual(node_1.GetSolutionStepValue(KratosMultiphysics.PRESSURE), node_2.GetSolutionStepValue(KratosMultiphysics.PRESSURE) * factor)
                self.assertVectorAlmostEqual(node_1.GetSolutionStepValue(KratosMultiphysics.VELOCITY), node_2.GetSolutionStepValue(KratosMultiphysics.VELOCITY) * factor)
                self.assertVectorAlmostEqual(node_1.GetSolutionStepValue(KratosMultiphysics.RECOVERED_STRESS), node_2.GetSolutionStepValue(KratosMultiphysics.RECOVERED_STRESS) * factor)
                self.assertMatrixAlmostEqual(node_1.GetSolutionStepValue(KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE), node_2.GetSolutionStepValue(KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE) * factor)

            for element_1, element_2 in zip(mp1.Elements, mp2.Elements):
                self.assertAlmostEqual(element_1.GetValue(KratosMultiphysics.PRESSURE), element_2.GetValue(KratosMultiphysics.PRESSURE) * factor)
                self.assertVectorAlmostEqual(element_1.GetValue(KratosMultiphysics.VELOCITY), element_2.GetValue(KratosMultiphysics.VELOCITY) * factor)
                self.assertVectorAlmostEqual(element_1.GetValue(KratosMultiphysics.RECOVERED_STRESS), element_2.GetValue(KratosMultiphysics.RECOVERED_STRESS) * factor)
                self.assertMatrixAlmostEqual(element_1.GetValue(KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE), element_2.GetValue(KratosMultiphysics.NORMAL_SHAPE_DERIVATIVE) * factor)

        model_1 = create_model_part(factor = 3)
        save_data_serializer = KratosMultiphysics.FileSerializer("serializer_data_only", KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE, True)
        save_data_serializer.Save("test", model_1["test"])

        # forcing the buffer to be written to the file.
        del save_data_serializer

        model_2 = create_model_part(factor = 6)
        comparator(model_1, model_2, 0.5)

        load_data_serializer = KratosMultiphysics.FileSerializer("serializer_data_only", KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE, True)
        load_data_serializer.Load("test", model_2["test"])
        comparator(model_1, model_2, 1.0)


if __name__ == '__main__':
    KratosUnittest.main()
