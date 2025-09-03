import typing
import math
import numpy as np
import xml.etree.ElementTree as ET
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess
from test_vtk_output_process import SetupModelPart2D, SetupModelPart3D

class TestVtuOutputBase:
    @classmethod
    def SetSolution(cls):
        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, 0,[node.X * 2, node.Y * 3, node.Z * 4])
            node.SetSolutionStepValue(Kratos.VELOCITY, 0,[node.X * 5, node.Y * 6, node.Z * 7])
            node.SetSolutionStepValue(Kratos.PRESSURE, 0, node.Id * 8.0)

        for i_elem, elem in enumerate(cls.model_part.Elements):
            elem.SetValue(Kratos.DETERMINANT, [i_elem*0.189, i_elem * 1.236, i_elem * 2.365])

        for i_cond, cond in enumerate(cls.model_part.Conditions):
            cond.SetValue(Kratos.DENSITY, i_cond * 4.362)
            cond.SetValue(Kratos.YOUNG_MODULUS, i_cond * 5.326)

    @classmethod
    def setUpClass(cls, output_prefix: str, setup_method, output_sub_model_parts: bool) -> None:
        cls.model =  Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.STEP] = 1
        cls.model_part.ProcessInfo[Kratos.TIME] = 1.0
        cls.output_prefix = output_prefix
        cls.output_sub_model_parts = output_sub_model_parts
        setup_method(cls.model_part)
        cls.SetSolution()

    def WriteVtu(self, output_format: Kratos.VtuOutput.WriterFormat):
        vtu_output = Kratos.VtuOutput(self.model_part, True, output_format, 9, echo_level=0, output_sub_model_parts=self.output_sub_model_parts)
        vtu_output.AddVariable(Kratos.PRESSURE, Kratos.Globals.DataLocation.NodeHistorical)
        vtu_output.AddVariable(Kratos.VELOCITY, Kratos.Globals.DataLocation.NodeHistorical)
        vtu_output.AddVariable(Kratos.DISPLACEMENT, Kratos.Globals.DataLocation.NodeHistorical)
        vtu_output.AddVariable(Kratos.DETERMINANT, Kratos.Globals.DataLocation.Element)
        vtu_output.AddVariable(Kratos.DENSITY, Kratos.Globals.DataLocation.Condition)
        vtu_output.AddVariable(Kratos.YOUNG_MODULUS, Kratos.Globals.DataLocation.Condition)

        a = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.PRESSURE, True)
        a *= 3
        vtu_output.AddContainerExpression("hist_exp", a)

        a = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.DETERMINANT)
        a *= 3
        vtu_output.AddContainerExpression("elem_exp", a)

        ta_1 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        ta_1.CollectData()
        ta_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.DETERMINANT)
        ta_2.CollectData()

        ta_1.data *= 3
        ta_2.data *= 3

        vtu_output.AddTensorAdaptor("hist_ta", ta_1)
        vtu_output.AddTensorAdaptor("elem_ta", ta_2)

        with kratos_unittest.WorkFolderScope("./auxiliar_files_for_python_unittest/vtk_output_process_ref_files", __file__, True):
            if output_format == Kratos.VtuOutput.ASCII:
                output_file_prefix = "ascii" + self.output_prefix + "/Main"
            else:
                output_file_prefix = "binary" + self.output_prefix + "/Main"
            vtu_output.PrintOutput("temp/" + output_file_prefix)
            self.Check("temp/" + output_file_prefix,  output_file_prefix)

    def test_WriteMeshAscii(self):
        self.WriteVtu(Kratos.VtuOutput.ASCII)

    def test_WriteMeshBinary(self):
        self.WriteVtu(Kratos.VtuOutput.BINARY)

    def Check(self, output_prefix, reference_prefix):
        def check_file(output_file_name: str, reference_file_name: str):
            ## Settings string in json format
            params = Kratos.Parameters("""{
                "reference_file_name" : "",
                "output_file_name"    : "",
                "comparison_type"     : "deterministic"
            }""")
            params["reference_file_name"].SetString(reference_file_name)
            params["output_file_name"].SetString(output_file_name)
            CompareTwoFilesCheckProcess(params).Execute()

        for file_path in Path(reference_prefix).iterdir():
            self.assertTrue((Path(output_prefix) / file_path.name).is_file())
            check_file(f"{output_prefix}/{file_path.name}", str(file_path))
        check_file(f"{output_prefix}.pvd", f"{reference_prefix}.pvd")


        kratos_utils.DeleteDirectoryIfExisting("temp")

class TestVtuOutput(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_location = Kratos.Globals.DataLocation

        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.STEP)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.CONSTITUTIVE_MATRIX)

        cls.model_part.ProcessInfo[Kratos.STEP] = 1
        cls.model_part.ProcessInfo[Kratos.TIME] = 1

        # creating a circular pie with origin (0, 0)
        radius = 1.0
        number_of_elements_or_conditions = 50
        angle = 2 * math.pi / number_of_elements_or_conditions

        cls.model_part.CreateNewNode(1, 0, 0, 0) # origin node
        for i in range(number_of_elements_or_conditions):
            theta = i * angle
            cls.model_part.CreateNewNode(i + 2, radius * math.cos(theta), radius * math.sin(theta), 0.0)

        prop = cls.model_part.CreateNewProperties(1)
        for i in range(number_of_elements_or_conditions):
            cls.model_part.CreateNewElement("Element2D3N", i + 1, [1, i + 2, (i + 1) % (number_of_elements_or_conditions) + 2], prop)
            cls.model_part.CreateNewCondition("LineCondition2D2N", i + 1, [i + 2, (i + 1) % (number_of_elements_or_conditions) + 2], prop)

        # create the sub model part structure
        for j in range(4):
            l_1 = j // 9
            l_2 = (j // 3) % 3
            l_3 = j % 3
            model_part = cls.model_part.CreateSubModelPart(f"mp_{l_1}.mp_{l_2}.mp_{l_3}")
            # fill the sub range
            for i in range(number_of_elements_or_conditions):
                if (j % 2 == 0):
                    if (i % ((l_1 + 1) * (l_2 + 1) * (l_3 + 1)) == 0):
                        model_part.AddCondition(cls.model_part.GetCondition(i + 1))
                else:
                    if (i % ((l_1 + 1) * (l_2 + 1) * (l_3 + 1)) == 1):
                        model_part.AddElement(cls.model_part.GetElement(i + 1))

        # now recursively fill the nodes
        def fill_nodes(model_part: Kratos.ModelPart):
            list_of_node_ids: 'list[int]' = []
            for condition in model_part.Conditions:
                for node in condition.GetGeometry():
                    list_of_node_ids.append(node.Id)
            for element in model_part.Elements:
                for node in element.GetGeometry():
                    list_of_node_ids.append(node.Id)

            model_part.AddNodes(list_of_node_ids)

            for sub_model_part in model_part.SubModelParts:
                fill_nodes(sub_model_part)

        fill_nodes(cls.model_part)

        def add_variables(container, setter):
            for entity in container:
                setter(entity, Kratos.STEP, entity.Id * 2)
                setter(entity, Kratos.PRESSURE, entity.Id)
                setter(entity, Kratos.DISPLACEMENT, Kratos.Array3([entity.Id, entity.Id + 1, entity.Id + 2]))
                setter(entity, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, Kratos.Vector([entity.Id, entity.Id + 1, entity.Id + 2, entity.Id + 3, entity.Id + 4]))
                setter(entity, Kratos.CONSTITUTIVE_MATRIX, Kratos.Matrix([[entity.Id, entity.Id + 1], [entity.Id + 2, entity.Id + 3], [entity.Id + 4, entity.Id + 5]]))

        add_variables(cls.model_part.Nodes, lambda x, y, z: x.SetSolutionStepValue(y, z))
        add_variables(cls.model_part.Nodes, lambda x, y, z: x.SetValue(y, z))
        add_variables(cls.model_part.Conditions, lambda x, y, z: x.SetValue(y, z))
        add_variables(cls.model_part.Elements, lambda x, y, z: x.SetValue(y, z))

    def test_PointVariableAddition(self):
        vtu_output = Kratos.VtuOutput(self.model_part, output_sub_model_parts=True)
        vtu_output.AddVariable(Kratos.PRESSURE, self.data_location.NodeHistorical)
        vtu_output.AddVariable(Kratos.PRESSURE, self.data_location.Condition)
        vtu_output.AddVariable(Kratos.PRESSURE, self.data_location.Element)
        vtu_output.AddIntegrationPointVariable(Kratos.PRESSURE, self.data_location.Condition)
        vtu_output.AddIntegrationPointVariable(Kratos.PRESSURE, self.data_location.Element)

        with self.assertRaises(RuntimeError):
            vtu_output.AddVariable(Kratos.PRESSURE, self.data_location.NodeHistorical)

        with self.assertRaises(RuntimeError):
            vtu_output.AddVariable(Kratos.PRESSURE, self.data_location.NodeNonHistorical)

        ta = Kratos.TensorAdaptors.DoubleTensorAdaptor(self.model_part.Nodes, Kratos.DoubleNDData([51, 5]))
        ta.Check()
        with self.assertRaises(RuntimeError):
            vtu_output.AddTensorAdaptor("PRESSURE", ta)
        vtu_output.AddTensorAdaptor("PRESSURE_1", ta)

        exp = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(exp, Kratos.PRESSURE, is_historical=False)
        with self.assertRaises(RuntimeError):
            vtu_output.AddContainerExpression("PRESSURE", exp)
        with self.assertRaises(RuntimeError):
            vtu_output.AddContainerExpression("PRESSURE_1", exp)
        vtu_output.AddContainerExpression("PRESSURE_2", exp)

        with self.assertRaises(RuntimeError):
            vtu_output.AddTensorAdaptor("PRESSURE_1", ta)
        with self.assertRaises(RuntimeError):
            vtu_output.AddTensorAdaptor("PRESSURE_2", ta)
        vtu_output.AddTensorAdaptor("PRESSURE_3", ta)

        ta = Kratos.TensorAdaptors.DoubleTensorAdaptor(self.model_part.Conditions, Kratos.DoubleNDData([50,2,4]))
        ta.Check()
        with self.assertRaises(RuntimeError):
            vtu_output.AddTensorAdaptor("PRESSURE", ta)
        vtu_output.AddTensorAdaptor("PRESSURE_1", ta)
        vtu_output.AddTensorAdaptor("PRESSURE_2", ta)
        vtu_output.AddTensorAdaptor("PRESSURE_3", ta)

        ta = Kratos.TensorAdaptors.DoubleTensorAdaptor(self.model_part.Elements, Kratos.DoubleNDData([50,2,5]))
        ta.Check()
        with self.assertRaises(RuntimeError):
            vtu_output.AddTensorAdaptor("PRESSURE", ta)
        vtu_output.AddTensorAdaptor("PRESSURE_1", ta)
        vtu_output.AddTensorAdaptor("PRESSURE_2", ta)
        vtu_output.AddTensorAdaptor("PRESSURE_3", ta)

        exp = Kratos.Expression.ConditionExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(exp, Kratos.PRESSURE)
        with self.assertRaises(RuntimeError):
            vtu_output.AddContainerExpression("PRESSURE_1", exp)
        vtu_output.AddContainerExpression("PRESSURE_4", exp)

        exp = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(exp, Kratos.PRESSURE)
        with self.assertRaises(RuntimeError):
            vtu_output.AddContainerExpression("PRESSURE_2", exp)
        vtu_output.AddContainerExpression("PRESSURE_4", exp)

        vtu_output.PrintOutput("temp/vtu_output/variable_test")

        model_part = vtu_output.GetModelPart()
        unstructured_grid_list = self.GetUnstructuredGridList(model_part, model_part.GetCommunicator().GetDataCommunicator(), True)
        for model_part, use_nodes, container in unstructured_grid_list:
            vtu_file_name = TestVtuOutput.GetUnstructuredGridName((model_part, use_nodes, container), "temp/vtu_output/variable_test", 0, model_part.GetCommunicator().GetDataCommunicator())
            if use_nodes:
                TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNodes(model_part, use_nodes, container), container, "binary",
                                        {
                                            "PRESSURE": (1, "Float64", None)
                                        },
                                        {
                                            "PRESSURE": (1, "Float64", None)
                                        })
            else:
                TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNodes(model_part, use_nodes, container), container, "binary",
                                        {},
                                        {
                                            "PRESSURE": (1, "Float64", None)
                                        })

            if model_part.FullName() == "test":
                if use_nodes:
                    TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNodes(model_part, use_nodes, container), container, "binary",
                                            {
                                                "PRESSURE": (1, "Float64", None),
                                                "PRESSURE_1": (5, "Float64", None),
                                                "PRESSURE_2": (1, "Float64", None),
                                                "PRESSURE_3": (5, "Float64", None)
                                            },
                                            {
                                                "PRESSURE": (1, "Float64", None)
                                            })
                elif isinstance(container, Kratos.ConditionsArray):
                    TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNodes(model_part, use_nodes, container), container, "binary",
                                            {},
                                            {
                                                "PRESSURE": (1, "Float64", None),
                                                "PRESSURE_1": (8, "Float64", None),
                                                "PRESSURE_2": (8, "Float64", None),
                                                "PRESSURE_3": (8, "Float64", None)
                                            })
                elif isinstance(container, Kratos.ElementsArray):
                    TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNodes(model_part, use_nodes, container), container, "binary",
                                            {
                                                "PRESSURE": (1, "Float64", None),
                                                "PRESSURE_1": (5, "Float64", None),
                                                "PRESSURE_2": (1, "Float64", None),
                                                "PRESSURE_3": (5, "Float64", None)
                                            },
                                            {
                                                "PRESSURE": (1, "Float64", None),
                                                "PRESSURE_1": (10, "Float64", None),
                                                "PRESSURE_2": (10, "Float64", None),
                                                "PRESSURE_3": (10, "Float64", None)
                                            })
            else:
                if use_nodes:
                    TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNodes(model_part, use_nodes, container), container, "binary",
                                            {
                                                "PRESSURE": (1, "Float64", None)
                                            },
                                            {
                                                "PRESSURE": (1, "Float64", None)
                                            })
                else:
                   TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNodes(model_part, use_nodes, container), container, "binary",
                                            {
                                            },
                                            {
                                                "PRESSURE": (1, "Float64", None)
                                            })

        kratos_utils.DeleteDirectoryIfExisting("temp/vtu_output/variable_test")

    def test_CellVariableAddition(self):
        vtu_output = Kratos.VtuOutput(self.model_part, binary_output=Kratos.VtuOutput.ASCII, output_sub_model_parts=True)
        vtu_output.AddVariable(Kratos.PRESSURE, self.data_location.Condition)
        vtu_output.AddVariable(Kratos.PRESSURE, self.data_location.NodeHistorical)
        vtu_output.AddVariable(Kratos.PRESSURE, self.data_location.Element)
        vtu_output.AddIntegrationPointVariable(Kratos.PRESSURE, self.data_location.Condition)
        vtu_output.AddIntegrationPointVariable(Kratos.PRESSURE, self.data_location.Element)

        vtu_output.AddVariable(Kratos.DISPLACEMENT, self.data_location.Element)
        vtu_output.AddVariable(Kratos.DISPLACEMENT, self.data_location.NodeHistorical)
        vtu_output.AddVariable(Kratos.DISPLACEMENT, self.data_location.Condition)
        vtu_output.AddIntegrationPointVariable(Kratos.DISPLACEMENT, self.data_location.Condition)
        vtu_output.AddIntegrationPointVariable(Kratos.DISPLACEMENT, self.data_location.Element)

        vtu_output.PrintOutput("temp/vtu_output/time_step_test")
        with self.assertRaises(RuntimeError):
            vtu_output.PrintOutput("temp/vtu_output/time_step_test")
        vtu_output.GetModelPart().ProcessInfo[Kratos.TIME] += 1e-9
        vtu_output.PrintOutput("temp/vtu_output/time_step_test")

        model_part = vtu_output.GetModelPart()
        unstructured_grid_list = self.GetUnstructuredGridList(model_part, model_part.GetCommunicator().GetDataCommunicator(), True)
        for step in range(2):
            for model_part, use_nodes, container in unstructured_grid_list:
                vtu_file_name = TestVtuOutput.GetUnstructuredGridName((model_part, use_nodes, container), "temp/vtu_output/time_step_test", step, model_part.GetCommunicator().GetDataCommunicator())
                if use_nodes:
                    TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNodes(model_part, use_nodes, container), container, "ascii",
                                            {
                                                "PRESSURE": (1, "Float64", None),
                                                "DISPLACEMENT": (3, "Float64", None)
                                            },
                                            {
                                                "PRESSURE": (1, "Float64", None),
                                                "DISPLACEMENT": (3, "Float64", None)
                                            })
                else:
                    TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNodes(model_part, use_nodes, container), container, "ascii",
                                            {
                                            },
                                            {
                                                "PRESSURE": (1, "Float64", None),
                                                "DISPLACEMENT": (3, "Float64", None)
                                            })

        kratos_utils.DeleteDirectoryIfExisting("temp/vtu_output/time_step_test")

    @staticmethod
    def GetUnstructuredGridList(model_part: Kratos.ModelPart, data_communicator: Kratos.DataCommunicator, recursively: bool) -> 'list[tuple[Kratos.ModelPart, bool, typing.Optional[typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]]]]':
        unstructured_grid_list: 'list[tuple[Kratos.ModelPart, bool, typing.Optional[typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]]]]' = []
        TestVtuOutput.__GetUnstructuredGridList(unstructured_grid_list, model_part, data_communicator, recursively)
        return unstructured_grid_list

    @staticmethod
    def GetUnstructuredGridName(unstructured_grid: 'tuple[Kratos.ModelPart, bool, typing.Optional[typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]]]', prefix: str, step: int, data_communicator: Kratos.DataCommunicator) -> str:
        if data_communicator.IsDistributed():
            rank_suffix = f"_{data_communicator.Rank()}"
        else:
            rank_suffix = ""

        if unstructured_grid[2] is None:
            entity_type = "nodes"
        elif isinstance(unstructured_grid[2], Kratos.ConditionsArray):
            entity_type = "conditions"
        elif isinstance(unstructured_grid[2], Kratos.ElementsArray):
            entity_type = "elements"
        return f"{prefix}/{unstructured_grid[0].FullName()}_{entity_type}_{step}{rank_suffix}.vtu"

    @staticmethod
    def GetNodes(model_part: Kratos.ModelPart, use_model_part_nodes: bool, container: 'typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]') -> 'list[Kratos.Node]':
        if use_model_part_nodes:
            return model_part.Nodes
        else:
            temp = []
            for entity in container:
                for node in entity.GetGeometry():
                    temp.append(node.Id)
            return [model_part.GetNode(node_id) for node_id in list(sorted(set(temp)))]

    @staticmethod
    def CheckVtuFile(
        test_class: kratos_unittest.TestCase,
        vtu_file_name: str,
        nodes: 'list[Kratos.Node]',
        container: 'typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray, list[int]]',
        data_format : str,
        point_data_fields: 'dict[str, tuple[int, str, typing.Any]]',
        cell_data_fields: 'dict[str, tuple[int, str, typing.Any]]'):

        def check_data_array(xml_element: ET.Element, number_of_components: int, name: str, data_type: str, data_format: str, data: 'typing.Optional[np.ndarray]' = None, tolerance = 1e-9):
            test_class.assertEqual(xml_element.tag, "DataArray")
            test_class.assertEqual(xml_element.get("NumberOfComponents"), str(number_of_components))
            test_class.assertEqual(xml_element.get("format"), data_format)
            test_class.assertEqual(xml_element.get("type"), data_type)
            test_class.assertEqual(xml_element.get("Name"), name)
            if data_format == "ascii" and not data is None:
                np_array = np.fromstring(xml_element.text, sep=" ")
                test_class.assertTrue(np.linalg.norm(data.ravel() - np_array) < tolerance)

        test_class.assertTrue(Path(vtu_file_name).is_file(), f"The file {vtu_file_name} not found.")
        tree = ET.parse(vtu_file_name)
        root = tree.getroot()

        test_class.assertEqual(root.tag, "VTKFile")
        test_class.assertEqual(root.get("type"), "UnstructuredGrid")
        test_class.assertEqual(root.get("version"), "0.1")

        unstructured_grid = root.find("UnstructuredGrid")
        piece = unstructured_grid.find("Piece")

        test_class.assertEqual(piece.get("NumberOfCells"), f"{len(container)}", f"Vtu file name = {vtu_file_name}")
        np_positions = np.zeros((len(nodes), 3), dtype=np.float64)
        for i, node in enumerate(nodes):
            np_positions[i, 0] = node.X0
            np_positions[i, 1] = node.Y0
            np_positions[i, 2] = node.Z0

        test_class.assertEqual(piece.get("NumberOfPoints"), f"{np_positions.shape[0]}", f"Vtu file name = {vtu_file_name}")

        points = piece.find("Points")
        check_data_array(points.find("DataArray"), 3, "Position", "Float64", data_format, np_positions)

        cells = piece.find("Cells")
        data_arrays = cells.findall("DataArray")
        np_offsets = np.zeros(len(container))
        total_connectivities = 0
        for i, entity in enumerate(container):
            total_connectivities += len(entity.GetGeometry())
            np_offsets[i] = total_connectivities
        check_data_array(data_arrays[1], 1, "offsets", "Int32", data_format, np_offsets)

        # now checking for point data
        point_data = piece.find("PointData")
        for data_field_name, (data_field_number_of_components, data_field_data_type, data) in point_data_fields.items():
            found = False
            for point_data_array in point_data.findall("DataArray"):
                if point_data_array.get("Name") == data_field_name:
                    check_data_array(point_data_array, data_field_number_of_components, data_field_name, data_field_data_type, data_format, data)
                    found = True
                    break
            test_class.assertTrue(found, f"Point data field \"{data_field_name}\" not found in the \"{vtu_file_name}\"")

        # now checking for cell data
        cell_data = piece.find("CellData")
        if container is not None:
            for data_field_name, (data_field_number_of_components, data_field_data_type, data) in cell_data_fields.items():
                found = False
                for point_data_array in cell_data.findall("DataArray"):
                    if point_data_array.get("Name") == data_field_name:
                        check_data_array(point_data_array, data_field_number_of_components, data_field_name, data_field_data_type, data_format, data)
                        found = True
                        break
                test_class.assertTrue(found, f"Cell data field \"{data_field_name}\" not found in the \"{vtu_file_name}\"")

    @staticmethod
    def CheckPvdFile(test_class: kratos_unittest.TestCase, pvt_file_name: str):
        test_class.assertTrue(Path(pvt_file_name).is_file(), f"The file {pvt_file_name} not found.")
        tree = ET.parse(pvt_file_name)
        root = tree.getroot()

    @staticmethod
    def __GetUnstructuredGridList(unstructured_grid_list: 'list[tuple[Kratos.ModelPart, bool, typing.Optional[typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]]]]', model_part: Kratos.ModelPart, data_communicator: Kratos.DataCommunicator, recursively: bool) -> None:
        availability = data_communicator.MaxAll(Kratos.Array3([model_part.NumberOfNodes(), model_part.NumberOfConditions(), model_part.NumberOfElements()]))
        has_nodes = availability[0] > 0
        has_conditions = availability[1] > 0
        has_elements = availability[2] > 0

        if has_elements:
            unstructured_grid_list.append((model_part, has_nodes, model_part.Elements))

        if has_conditions:
            unstructured_grid_list.append((model_part, not has_elements and has_nodes, model_part.Conditions))

        if not has_elements and not has_conditions:
            unstructured_grid_list.append((model_part, has_nodes, None))

        if recursively:
            for sub_model_part in model_part.SubModelParts:
                TestVtuOutput.__GetUnstructuredGridList(unstructured_grid_list, sub_model_part, data_communicator, recursively)

class TestVtuOutput2D(TestVtuOutputBase, kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass("2D", SetupModelPart2D, output_sub_model_parts = True)

class TestVtuOutput3D(TestVtuOutputBase, kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        # the method SetupModelPart3D does not create sub model parts
        # with nodes which do not include nodes from its conditions or elements. It uses
        # some random nodes. Hence sub_model_part output is disabled.
        super().setUpClass("3D", SetupModelPart3D, output_sub_model_parts = False)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.INFO)
    kratos_unittest.main()