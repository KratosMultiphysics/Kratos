import os
import pickle
import numpy

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


    def test_TensorAdaptorSerialization(self):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("Test")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        for i in range(3):
            node = model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetValue(KratosMultiphysics.PRESSURE, i + 1.0)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 10.0 * (i + 1))

        ta = KratosMultiphysics.TensorAdaptors
        variable_ta = ta.VariableTensorAdaptor(model_part.Nodes, KratosMultiphysics.PRESSURE)
        variable_ta.CollectData()
        historical_ta = ta.HistoricalVariableTensorAdaptor(model_part.Nodes, KratosMultiphysics.TEMPERATURE)
        historical_ta.CollectData()
        combined = ta.DoubleCombinedTensorAdaptor([variable_ta, historical_ta], False, False, False)
        combined.CollectData()

        serializer = KratosMultiphysics.StreamSerializer()
        serializer.Set(KratosMultiphysics.Serializer.SHALLOW_GLOBAL_POINTERS_SERIALIZATION)
        serializer.Save("Model", model)
        serializer.Save("Variable", variable_ta)
        serializer.Save("Combined", combined)

        loaded_model = KratosMultiphysics.Model()
        loaded_variable = ta.VariableTensorAdaptor()
        loaded_combined = ta.DoubleCombinedTensorAdaptor()
        serializer.Load("Model", loaded_model)
        serializer.Load("Variable", loaded_variable)
        serializer.Load("Combined", loaded_combined)

        self.assertVectorAlmostEqual(loaded_variable.data, variable_ta.data)
        self.assertVectorAlmostEqual(loaded_combined.data, combined.data)
        children = loaded_combined.GetTensorAdaptors()
        self.assertIsInstance(children[0], ta.VariableTensorAdaptor)
        self.assertIsInstance(children[1], ta.HistoricalVariableTensorAdaptor)

        # the container must point into the loaded model, not a copy: StoreData() changes the
        # loaded nodes and leaves the original model untouched.
        loaded_variable.data[:] = [7.0, 8.0, 9.0]
        loaded_variable.StoreData()
        loaded_model_part = loaded_model.GetModelPart("Test")
        self.assertEqual([node.GetValue(KratosMultiphysics.PRESSURE) for node in loaded_model_part.Nodes], [7.0, 8.0, 9.0])
        self.assertEqual([node.GetValue(KratosMultiphysics.PRESSURE) for node in model_part.Nodes], [1.0, 2.0, 3.0])

    def test_NestedCombinedTensorAdaptorSerialization(self):
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("Test")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        for i in range(3):
            node = model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetValue(KratosMultiphysics.PRESSURE, i + 1.0)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 10.0 * (i + 1))

        # 3 levels: [[[nodal PRESSURE, historical TEMPERATURE], node position], historical TEMPERATURE],
        # every level collects and stores recursively.
        ta = KratosMultiphysics.TensorAdaptors
        def CreateNested(model_part):
            level_3 = ta.DoubleCombinedTensorAdaptor([ta.VariableTensorAdaptor(model_part.Nodes, KratosMultiphysics.PRESSURE), ta.HistoricalVariableTensorAdaptor(model_part.Nodes, KratosMultiphysics.TEMPERATURE)], True, True, False)
            level_2 = ta.DoubleCombinedTensorAdaptor([level_3, ta.NodePositionTensorAdaptor(model_part.Nodes, KratosMultiphysics.Configuration.Current)], True, True, False)
            return ta.DoubleCombinedTensorAdaptor([level_2, ta.HistoricalVariableTensorAdaptor(model_part.Nodes, KratosMultiphysics.TEMPERATURE)], True, True, False)

        combined = CreateNested(model_part)
        combined.CollectData()
        original_data = numpy.array(combined.data)

        serializer = KratosMultiphysics.StreamSerializer()
        serializer.Set(KratosMultiphysics.Serializer.SHALLOW_GLOBAL_POINTERS_SERIALIZATION)
        serializer.Save("Model", model)
        serializer.Save("Combined", combined)

        loaded_model = KratosMultiphysics.Model()
        loaded = ta.DoubleCombinedTensorAdaptor()
        serializer.Load("Model", loaded_model)
        serializer.Load("Combined", loaded)

        self.assertVectorAlmostEqual(loaded.data, original_data)
        level_1 = loaded.GetTensorAdaptors()
        self.assertEqual([type(child) for child in level_1], [ta.DoubleCombinedTensorAdaptor, ta.HistoricalVariableTensorAdaptor])
        level_2 = level_1[0].GetTensorAdaptors()
        self.assertEqual([type(child) for child in level_2], [ta.DoubleCombinedTensorAdaptor, ta.NodePositionTensorAdaptor])
        level_3 = level_2[0].GetTensorAdaptors()
        self.assertEqual([type(child) for child in level_3], [ta.VariableTensorAdaptor, ta.HistoricalVariableTensorAdaptor])

        # the recursive StoreData() reaches the loaded model through all levels, not the original model
        loaded.data[:] = -original_data
        loaded.StoreData()
        check = CreateNested(loaded_model.GetModelPart("Test"))
        check.CollectData()
        self.assertVectorAlmostEqual(check.data, -original_data)
        check = CreateNested(model_part)
        check.CollectData()
        self.assertVectorAlmostEqual(check.data, original_data)

    def test_TensorAdaptorLoadKeepsNumpyViews(self):
        # A same-size load overwrites the existing buffer in place, so numpy arrays sharing it stay
        # valid and see the loaded values, and a buffer owned by numpy is written back into.
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("Test")
        for i in range(3):
            model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0).SetValue(KratosMultiphysics.PRESSURE, i + 1.0)

        ta = KratosMultiphysics.TensorAdaptors
        variable_ta = ta.VariableTensorAdaptor(model_part.Nodes, KratosMultiphysics.PRESSURE)
        variable_ta.CollectData()
        view = variable_ta.data

        np_array = numpy.array([4.0, 5.0, 6.0])
        nd_data_ta = ta.DoubleTensorAdaptor(model_part.Nodes, KratosMultiphysics.DoubleNDData(np_array, copy=False), copy=False)

        serializer = KratosMultiphysics.StreamSerializer()
        serializer.Set(KratosMultiphysics.Serializer.SHALLOW_GLOBAL_POINTERS_SERIALIZATION)
        serializer.Save("Model", model)
        serializer.Save("Variable", variable_ta)
        serializer.Save("NDData", nd_data_ta)

        variable_ta.data[:] = 0.0
        np_array[:] = 0.0

        loaded_model = KratosMultiphysics.Model()
        serializer.Load("Model", loaded_model)
        serializer.Load("Variable", variable_ta)
        serializer.Load("NDData", nd_data_ta)

        self.assertVectorAlmostEqual(view, [1.0, 2.0, 3.0])
        self.assertVectorAlmostEqual(variable_ta.data, [1.0, 2.0, 3.0])
        self.assertVectorAlmostEqual(np_array, [4.0, 5.0, 6.0])
        self.assertVectorAlmostEqual(nd_data_ta.data, [4.0, 5.0, 6.0])

    def test_TensorAdaptorRepeatedLoadIntoSameAdaptor(self):
        # A non-DataOnly load into a container-bound adaptor reloads its nodes' historical and
        # non-historical data in place, so every repeated load must restore the saved state.
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("Test")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.SetBufferSize(2)
        for i in range(3):
            node = model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetValue(KratosMultiphysics.PRESSURE, i + 1.0)
            node.SetValue(KratosMultiphysics.VELOCITY, KratosMultiphysics.Array3([i, 2.0 * i, 3.0 * i]))
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0, 10.0 * (i + 1))
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, 20.0 * (i + 1))
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, KratosMultiphysics.Array3([i, -i, 0.5 * i]))

        def GetNodalState():
            return [(node.GetValue(KratosMultiphysics.PRESSURE), list(node.GetValue(KratosMultiphysics.VELOCITY)),
                     node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0), node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1),
                     list(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0))) for node in model_part.Nodes]

        ta = KratosMultiphysics.TensorAdaptors
        variable_ta = ta.VariableTensorAdaptor(model_part.Nodes, KratosMultiphysics.VELOCITY)
        variable_ta.CollectData()
        original_data = numpy.array(variable_ta.data)
        original_state = GetNodalState()

        saved = KratosMultiphysics.StreamSerializer()
        saved.Save("Variable", variable_ta)

        for _ in range(3):
            variable_ta.data[:] = 0.0
            for node in model_part.Nodes:
                node.SetValue(KratosMultiphysics.PRESSURE, -1.0)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0, -1.0)
            pickle.loads(pickle.dumps(saved)).Load("Variable", variable_ta) # a fresh copy of the saved stream
            self.assertVectorAlmostEqual(variable_ta.data.ravel(), original_data.ravel())
            self.assertEqual(GetNodalState(), original_state)

        # the adaptor still points to the model part's nodes
        self.assertEqual(model_part.NumberOfNodes(), 3)
        variable_ta.data[:] = 7.0
        variable_ta.StoreData()
        self.assertEqual([list(node.GetValue(KratosMultiphysics.VELOCITY)) for node in model_part.Nodes], [[7.0] * 3] * 3)

    def test_TensorAdaptorLoadIntoNodesWithoutHistoricalVariables(self):
        # The nodes of a model part without historical variables share an empty variables list. An in-place load of
        # nodes with historical variables must raise without adding them to it, as the other nodes can't hold them.
        model = KratosMultiphysics.Model()
        saved_model_part = model.CreateModelPart("Saved")
        saved_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        saved_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LOCAL_AXES_MATRIX)
        model_part = model.CreateModelPart("Test")
        sub_model_part = model_part.CreateSubModelPart("Sub")
        for i in range(3):
            node = saved_model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(KratosMultiphysics.LOCAL_AXES_MATRIX, KratosMultiphysics.Matrix(3, 3, 1.0))
            sub_model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            model_part.CreateNewNode(i + 4, float(i), 1.0, 0.0) # not in the archive

        ta = KratosMultiphysics.TensorAdaptors
        serializer = KratosMultiphysics.StreamSerializer()
        serializer.Save("Variable", ta.VariableTensorAdaptor(saved_model_part.Nodes, KratosMultiphysics.VELOCITY))
        with self.assertRaisesRegex(RuntimeError, "The archive does not match this VariablesList"):
            serializer.Load("Variable", ta.VariableTensorAdaptor(sub_model_part.Nodes, KratosMultiphysics.VELOCITY))
        self.assertFalse(model_part.HasNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE))
        self.assertFalse(model_part.HasNodalSolutionStepVariable(KratosMultiphysics.LOCAL_AXES_MATRIX))
        model.DeleteModelPart("Test") # destructs every node's data with the (unchanged) empty list

    def test_UninitializedTensorAdaptor(self):
        # A default-constructed TA is only a load target: its default storage's PointerWrapper
        # raises on any data access, including the ones made through str(), Check() and CollectData().
        ta = KratosMultiphysics.TensorAdaptors
        for tensor_adaptor_type in [ta.DoubleTensorAdaptor, ta.VariableTensorAdaptor, ta.ConnectivityIdsTensorAdaptor]:
            tensor_adaptor = tensor_adaptor_type()
            operations = {
                "data": lambda: tensor_adaptor.data,
                "str": lambda: str(tensor_adaptor),
                "Shape": lambda: tensor_adaptor.Shape(),
                "Check": lambda: tensor_adaptor.Check(),
                "CollectData": lambda: tensor_adaptor.CollectData(),
                "Save": lambda: KratosMultiphysics.StreamSerializer().Save("TA", tensor_adaptor)}
            if tensor_adaptor_type is ta.DoubleTensorAdaptor:
                del operations["CollectData"] # the base class raises its own "can only be used for data storage" error
            for name, operation in operations.items():
                with self.subTest(type=tensor_adaptor_type.__name__, operation=name):
                    with self.assertRaisesRegex(RuntimeError, "Uninitialized PointerWrapper"):
                        operation()

if __name__ == '__main__':
    KratosUnittest.main()
