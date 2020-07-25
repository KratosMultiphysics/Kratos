import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication import CoSimIO
from KratosMultiphysics import kratos_utilities as kratos_utils
from KratosMultiphysics.testing.utilities import GetPython3Command

import os
import subprocess


class TestCoSimIOPyExposure_aux_tests(KratosUnittest.TestCase):

    def test_InfoFromParameters(self):
        params = KM.Parameters("""{
            "some_bool"      : true,
            "the_string"     : "file",
            "another_string" : "i25mmk",
            "int_for_it"     : 12,
            "double_val"     : 0.223,
            "array_val_ign"  : [0.923, 1.8]
        }""")

        info = CoSimIO.InfoFromParameters(params)
        self.assertEqual(info.Size(), 5)
        self.assertTrue(info.Has("some_bool"))
        self.assertTrue(info.Has("the_string"))
        self.assertTrue(info.Has("another_string"))
        self.assertTrue(info.Has("int_for_it"))
        self.assertTrue(info.Has("double_val"))
        # only int, bool, double & string can be converted, others are ignored
        self.assertFalse(info.Has("array_val_ign"))

        self.assertTrue(info.GetBool("some_bool"))
        self.assertEqual(info.GetString("the_string"), "file")
        self.assertEqual(info.GetString("another_string"), "i25mmk")
        self.assertEqual(info.GetInt("int_for_it"), 12)
        self.assertAlmostEqual(info.GetDouble("double_val"), 0.223)


class TestCoSimIOPyExposure(KratosUnittest.TestCase):

    def test_Connect_Disconnect(self):
        connection_settings = CoSimIO.Info()
        connection_settings.SetString("connection_name", "c_d_test")
        connection_settings.SetInt("echo_level", 0)
        info = CoSimIO.Connect(connection_settings)
        self.assertEqual(info.GetInt("connection_status"), CoSimIO.ConnectionStatus.Connected)

        RunPythonInSubProcess("connect_disconnect")

        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", "c_d_test")
        info = CoSimIO.Disconnect(disconnect_settings)
        self.assertEqual(info.GetInt("connection_status"), CoSimIO.ConnectionStatus.Disconnected)

    def test_Export_Import_Data_raw_values(self):
        connection_settings = CoSimIO.Info()
        connection_settings.SetString("connection_name", "im_exp_data")
        connection_settings.SetInt("echo_level", 0)
        info = CoSimIO.Connect(connection_settings)
        self.assertEqual(info.GetInt("connection_status"), CoSimIO.ConnectionStatus.Connected)

        values = [1.0, 2.5, 3.3, -9.4]

        export_info = CoSimIO.Info()
        export_info.SetString("connection_name", "im_exp_data")
        export_info.SetString("identifier", "data_exchange_1")
        CoSimIO.ExportData(export_info, values)

        RunPythonInSubProcess("import_export_data")

        import_info = CoSimIO.Info()
        import_info.SetString("connection_name", "im_exp_data")
        import_info.SetString("identifier", "data_exchange_2")
        imported_values = CoSimIO.ImportData(import_info)

        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", "im_exp_data")

        info = CoSimIO.Disconnect(disconnect_settings)
        self.assertEqual(info.GetInt("connection_status"), CoSimIO.ConnectionStatus.Disconnected)

        # checking the values after disconnecting to avoid deadlock
        self.assertVectorAlmostEqual(KM.Vector(values), KM.Vector(imported_values))

    def test_Export_Import_Data_ModelPart_scalar_node_historical(self):
        model = KM.Model()

        model_part = model.CreateModelPart("for_test")
        model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        model_part.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        for i in range(5):
            node = model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0) # using same coord, doesn't matter for this test
            node.SetSolutionStepValue(KM.PRESSURE, i*1.7)

        self.__ExportImportDataOnModelPart(model_part, KM.PRESSURE, KM.TEMPERATURE, CoSimIO.DataLocation.NodeHistorical)

        # checking the values after disconnecting to avoid deadlock
        for i, node in enumerate(model_part.Nodes):
            self.assertAlmostEqual(node.GetSolutionStepValue(KM.TEMPERATURE), i*1.7)

    def test_Export_Import_Data_ModelPart_vector_node_historical(self):
        model = KM.Model()

        model_part = model.CreateModelPart("for_test")
        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        for i in range(5):
            node = model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0) # using same coord, doesn't matter for this test
            node.SetSolutionStepValue(KM.DISPLACEMENT, [i*1.7, i+1.1, i**1.2])

        self.__ExportImportDataOnModelPart(model_part, KM.DISPLACEMENT, KM.VELOCITY, CoSimIO.DataLocation.NodeHistorical)

        # checking the values after disconnecting to avoid deadlock
        for i, node in enumerate(model_part.Nodes):
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(KM.VELOCITY), KM.Vector([i*1.7, i+1.1, i**1.2]))

    def test_Export_Import_Data_ModelPart_scalar_node_non_historical(self):
        model = KM.Model()

        model_part = model.CreateModelPart("for_test")
        model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        model_part.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        for i in range(5):
            node = model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0) # using same coord, doesn't matter for this test
            node.SetValue(KM.PRESSURE, i*1.6)

        self.__ExportImportDataOnModelPart(model_part, KM.PRESSURE, KM.TEMPERATURE, CoSimIO.DataLocation.NodeNonHistorical)

        # checking the values after disconnecting to avoid deadlock
        for i, node in enumerate(model_part.Nodes):
            self.assertAlmostEqual(node.GetValue(KM.TEMPERATURE), i*1.6)

    def test_Export_Import_Data_ModelPart_vector_node_non_historical(self):
        model = KM.Model()

        model_part = model.CreateModelPart("for_test")
        for i in range(5):
            node = model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0) # using same coord, doesn't matter for this test
            node.SetValue(KM.DISPLACEMENT, [i*1.7, i+1.99, i**1.2])

        self.__ExportImportDataOnModelPart(model_part, KM.DISPLACEMENT, KM.VELOCITY, CoSimIO.DataLocation.NodeNonHistorical)

        # checking the values after disconnecting to avoid deadlock
        for i, node in enumerate(model_part.Nodes):
            self.assertVectorAlmostEqual(node.GetValue(KM.VELOCITY), KM.Vector([i*1.7, i+1.99, i**1.2]))

    def test_Export_Import_Data_ModelPart_scalar_element(self):
        model = KM.Model()

        model_part = model.CreateModelPart("for_test")
        for i in range(5):
            model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0) # using same coord, doesn't matter for this test
        props = model_part.CreateNewProperties(0)
        for i in range(4):
            element = model_part.CreateNewElement("Element2D2N", i+1, [i+1, i+2], props)
            element.SetValue(KM.PRESSURE, i*1.6)

        self.__ExportImportDataOnModelPart(model_part, KM.PRESSURE, KM.TEMPERATURE, CoSimIO.DataLocation.Element)

        # checking the values after disconnecting to avoid deadlock
        for i, elem in enumerate(model_part.Elements):
            self.assertAlmostEqual(elem.GetValue(KM.TEMPERATURE), i*1.6)

    def test_Export_Import_Data_ModelPart_vector_element(self):
        model = KM.Model()

        model_part = model.CreateModelPart("for_test")
        for i in range(5):
            model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0) # using same coord, doesn't matter for this test
        props = model_part.CreateNewProperties(0)
        for i in range(4):
            element = model_part.CreateNewElement("Element2D2N", i+1, [i+1, i+2], props)
            element.SetValue(KM.DISPLACEMENT, [i*1.7, i+1.99, i**1.2])

        self.__ExportImportDataOnModelPart(model_part, KM.DISPLACEMENT, KM.VELOCITY, CoSimIO.DataLocation.Element)

        # checking the values after disconnecting to avoid deadlock
        for i, elem in enumerate(model_part.Elements):
            self.assertVectorAlmostEqual(elem.GetValue(KM.VELOCITY), KM.Vector([i*1.7, i+1.99, i**1.2]))

    def test_Export_Import_Data_ModelPart_scalar_condition(self):
        model = KM.Model()

        model_part = model.CreateModelPart("for_test")
        for i in range(5):
            model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0) # using same coord, doesn't matter for this test
        props = model_part.CreateNewProperties(0)
        for i in range(5):
            condition = model_part.CreateNewCondition("PointCondition2D1N", i+1, [i+1], props)
            condition.SetValue(KM.PRESSURE, i*(-11.6))

        self.__ExportImportDataOnModelPart(model_part, KM.PRESSURE, KM.TEMPERATURE, CoSimIO.DataLocation.Condition)

        # checking the values after disconnecting to avoid deadlock
        for i, cond in enumerate(model_part.Conditions):
            self.assertAlmostEqual(cond.GetValue(KM.TEMPERATURE), i*(-11.6))

    def test_Export_Import_Data_ModelPart_vector_condition(self):
        model = KM.Model()

        model_part = model.CreateModelPart("for_test")
        for i in range(5):
            model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0) # using same coord, doesn't matter for this test
        props = model_part.CreateNewProperties(0)
        for i in range(5):
            condition = model_part.CreateNewCondition("PointCondition2D1N", i+1, [i+1], props)
            condition.SetValue(KM.DISPLACEMENT, [i*1.7, i+1.25, i**1.2])

        self.__ExportImportDataOnModelPart(model_part, KM.DISPLACEMENT, KM.VELOCITY, CoSimIO.DataLocation.Condition)

        # checking the values after disconnecting to avoid deadlock
        for i, cond in enumerate(model_part.Conditions):
            self.assertVectorAlmostEqual(cond.GetValue(KM.VELOCITY), KM.Vector([i*1.7, i+1.25, i**1.2]))

    def test_Export_Import_Data_ModelPart_scalar_model_part(self):
        model = KM.Model()

        model_part = model.CreateModelPart("for_test")
        model_part.SetValue(KM.PRESSURE, 333.896)

        self.__ExportImportDataOnModelPart(model_part, KM.PRESSURE, KM.TEMPERATURE, CoSimIO.DataLocation.ModelPart)

        # checking the values after disconnecting to avoid deadlock
        self.assertAlmostEqual(model_part.GetValue(KM.TEMPERATURE), 333.896)

    def test_Export_Import_Data_ModelPart_vector_model_part(self):
        model = KM.Model()

        model_part = model.CreateModelPart("for_test")
        model_part.SetValue(KM.DISPLACEMENT, [14.3, -333.896, 987.2])

        self.__ExportImportDataOnModelPart(model_part, KM.DISPLACEMENT, KM.VELOCITY, CoSimIO.DataLocation.ModelPart)

        # checking the values after disconnecting to avoid deadlock
        self.assertVectorAlmostEqual(model_part.GetValue(KM.VELOCITY), KM.Vector([14.3, -333.896, 987.2]))

    def test_Export_Import_Mesh(self):
        model = KM.Model()

        model_part = model.CreateModelPart("for_test")
        model_part_returned = model.CreateModelPart("for_test_2")

        for i in range(15):
            model_part.CreateNewNode(i+1, i*1.1, 1.1**i, 0.0)
        props = model_part.CreateNewProperties(0)
        for i in range(9): # this leaves some hanging nodes
            model_part.CreateNewElement("Element2D2N", i+1, [i+1, i+2], props)

        connection_settings = CoSimIO.Info()
        connection_settings.SetString("connection_name", "im_exp_mesh")
        connection_settings.SetInt("echo_level", 0)
        info = CoSimIO.Connect(connection_settings)
        self.assertEqual(info.GetInt("connection_status"), CoSimIO.ConnectionStatus.Connected)

        export_info = CoSimIO.Info()
        export_info.SetString("connection_name", "im_exp_mesh")
        export_info.SetString("identifier", "mesh_exchange_1")
        CoSimIO.ExportMesh(export_info, model_part)

        RunPythonInSubProcess("import_export_mesh")

        import_info = CoSimIO.Info()
        import_info.SetString("connection_name", "im_exp_mesh")
        import_info.SetString("identifier", "mesh_exchange_2")
        CoSimIO.ImportMesh(import_info, model_part_returned)

        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", "im_exp_mesh")

        info = CoSimIO.Disconnect(disconnect_settings)
        self.assertEqual(info.GetInt("connection_status"), CoSimIO.ConnectionStatus.Disconnected)

        # checking the values after disconnecting to avoid deadlock
        self.assertEqual(model_part.NumberOfNodes(), model_part_returned.NumberOfNodes())
        self.assertEqual(model_part.NumberOfElements(), model_part_returned.NumberOfElements())

        for node_orig, node_exchanged in zip(model_part.Nodes, model_part_returned.Nodes):
            self.assertAlmostEqual(node_orig.X0, node_exchanged.X0)
            self.assertAlmostEqual(node_orig.Y0, node_exchanged.Y0)
            self.assertAlmostEqual(node_orig.Z0, node_exchanged.Z0)

        for elem_orig, elem_exchanged in zip(model_part.Elements, model_part_returned.Elements):
            self.assertEqual(len(elem_orig.GetNodes()), len(elem_exchanged.GetNodes()))
            for node_orig, node_exchanged in zip(elem_orig.GetNodes(), elem_exchanged.GetNodes()):
                self.assertAlmostEqual(node_orig.X0, node_exchanged.X0)
                self.assertAlmostEqual(node_orig.Y0, node_exchanged.Y0)
                self.assertAlmostEqual(node_orig.Z0, node_exchanged.Z0)


    def __ExportImportDataOnModelPart(self, model_part, var_export, var_import, data_location):
        connection_settings = CoSimIO.Info()
        connection_settings.SetString("connection_name", "im_exp_data")
        connection_settings.SetInt("echo_level", 0)
        info = CoSimIO.Connect(connection_settings)
        self.assertEqual(info.GetInt("connection_status"), CoSimIO.ConnectionStatus.Connected)

        export_info = CoSimIO.Info()
        export_info.SetString("connection_name", "im_exp_data")
        export_info.SetString("identifier", "data_exchange_1")
        CoSimIO.ExportData(export_info, model_part, var_export, data_location)

        RunPythonInSubProcess("import_export_data")

        import_info = CoSimIO.Info()
        import_info.SetString("connection_name", "im_exp_data")
        import_info.SetString("identifier", "data_exchange_2")
        CoSimIO.ImportData(import_info, model_part, var_import, data_location)

        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", "im_exp_data")

        info = CoSimIO.Disconnect(disconnect_settings)
        self.assertEqual(info.GetInt("connection_status"), CoSimIO.ConnectionStatus.Disconnected)

def RunPythonInSubProcess(python_script_name):
    if not python_script_name.endswith(".py"):
        python_script_name += ".py"

    py_cmd = GetPython3Command()

    cmd_list = [py_cmd, os.path.join("co_sim_io_py_exposure_aux_files", python_script_name)]
    subprocess.run(cmd_list, check=True, shell=os.name=="nt") # crashes the calling script too, otherwise the error is silent (using shell in Win)


if __name__ == '__main__':
    KratosUnittest.main()
