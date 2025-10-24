import math, typing
from pathlib import Path
import xml.etree.ElementTree as ET

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

# import from the tests
with kratos_unittest.WorkFolderScope("../../tests", __file__, True):
    from test_vtk_output_process import SetupModelPart2D, SetupModelPart3D

class TestVtuOutputBase:
    @classmethod
    def SetSolution(cls):
        node: Kratos.Node
        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, 0,[node.X * 2, node.Y * 3, node.Z * 4])
            node.SetSolutionStepValue(Kratos.VELOCITY, 0,[node.X * 5, node.Y * 6, node.Z * 7])
            node.SetSolutionStepValue(Kratos.PRESSURE, 0, node.Id * 8.0)

        elem: Kratos.Element
        for i_elem, elem in enumerate(cls.model_part.Elements):
            elem.SetValue(Kratos.DETERMINANT, [i_elem*0.189, i_elem * 1.236, i_elem * 2.365])

        cond: Kratos.Condition
        for i_cond, cond in enumerate(cls.model_part.Conditions):
            cond.SetValue(Kratos.DENSITY, i_cond * 4.362)
            cond.SetValue(Kratos.YOUNG_MODULUS, i_cond * 5.326)
            cond.SetValue(Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, Kratos.Vector())

    @classmethod
    def setUpClass(cls, output_prefix: str, setup_method, output_sub_model_parts: bool) -> None:
        cls.model =  Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.STEP] = 0
        cls.model_part.ProcessInfo[Kratos.TIME] = 1.0
        cls.output_prefix = output_prefix
        cls.output_sub_model_parts = output_sub_model_parts
        setup_method(cls.model_part)
        cls.SetSolution()

    def WriteVtu(self, output_format: Kratos.Future.VtuOutput.WriterFormat, error_check = False):
        vtu_output = Kratos.Future.VtuOutput(self.model_part, Kratos.Configuration.Initial, output_format, 9, echo_level=0, output_sub_model_parts=self.output_sub_model_parts, write_ids=True)
        vtu_output.AddVariable(Kratos.PRESSURE, Kratos.Globals.DataLocation.NodeHistorical)
        vtu_output.AddVariable(Kratos.VELOCITY, Kratos.Globals.DataLocation.NodeHistorical)
        vtu_output.AddVariable(Kratos.DISPLACEMENT, Kratos.Globals.DataLocation.NodeHistorical)
        vtu_output.AddVariable(Kratos.DETERMINANT, Kratos.Globals.DataLocation.Element)
        vtu_output.AddVariable(Kratos.DENSITY, Kratos.Globals.DataLocation.Condition)
        vtu_output.AddVariable(Kratos.YOUNG_MODULUS, Kratos.Globals.DataLocation.Condition)
        if error_check:
            # adds a vector with zero size to check for the error. Because, in vtu lib,
            # if the behavior is undefined if the "NumberOfComponents = 0".
            vtu_output.AddVariable(Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, Kratos.Globals.DataLocation.Condition)

        ta_1 = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        ta_1.CollectData()
        ta_2 = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.DETERMINANT)
        ta_2.CollectData()

        ta_1.data *= 3
        ta_2.data *= 3

        vtu_output.AddTensorAdaptor("hist_ta", ta_1)
        vtu_output.AddTensorAdaptor("elem_ta", ta_2)

        with kratos_unittest.WorkFolderScope("vtk_output_process_ref_files", __file__, True):
            output_file_prefix = output_format.name.lower() + self.output_prefix + "/Main"
            vtu_output.PrintOutput("temp/" + output_file_prefix)
            self.Check("temp/" + output_file_prefix,  output_file_prefix)

    def test_WriteMeshAscii(self):
        self.WriteVtu(Kratos.Future.VtuOutput.ASCII)

    def test_WriteMeshBinary(self):
        self.WriteVtu(Kratos.Future.VtuOutput.BINARY)

    def test_WriteMeshRaw(self):
        self.WriteVtu(Kratos.Future.VtuOutput.RAW)

    def test_WriteMeshCompressedRaw(self):
        self.WriteVtu(Kratos.Future.VtuOutput.COMPRESSED_RAW)

    def test_WriteMeshAsciiWithError(self):
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, "vtk_output_process_ref_files/temp")
        with self.assertRaises(RuntimeError):
            self.WriteVtu(Kratos.Future.VtuOutput.ASCII, True)

    def test_WriteMeshBinaryWithError(self):
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, "vtk_output_process_ref_files/temp")
        with self.assertRaises(RuntimeError):
            self.WriteVtu(Kratos.Future.VtuOutput.BINARY, True)

    def test_WriteMeshRawWithError(self):
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, "vtk_output_process_ref_files/temp")
        with self.assertRaises(RuntimeError):
            self.WriteVtu(Kratos.Future.VtuOutput.RAW, True)

    def test_WriteMeshCompressedRawWithError(self):
        self.addCleanup(kratos_utils.DeleteDirectoryIfExisting, "vtk_output_process_ref_files/temp")
        with self.assertRaises(RuntimeError):
            self.WriteVtu(Kratos.Future.VtuOutput.COMPRESSED_RAW, True)

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
            self.assertTrue((Path(output_prefix) / file_path.name).is_file(), msg=f"\"{(Path(output_prefix) / file_path.name)}\" is not a file.")
            check_file(f"{output_prefix}/{file_path.name}", str(file_path))
        check_file(f"{output_prefix}.pvd", f"{reference_prefix}.pvd")

        kratos_utils.DeleteDirectoryIfExistingAndEmpty("temp")

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

    def setUp(self):
        self.model_part.ProcessInfo[Kratos.STEP] = 0

    def test_GetOutputContainerList(self):
        unstructured_grid_list = self.GetUnstructuredGridList(self.model_part, self.model_part.GetRootModelPart().GetCommunicator().GetDataCommunicator(), True)

        vtu_output = Kratos.Future.VtuOutput(self.model_part, output_sub_model_parts=True, output_format=Kratos.Future.VtuOutput.BINARY)
        list_of_containers = vtu_output.GetOutputContainerList()

        list_of_ref_containers = []
        for mp, is_nodes_used, cells_container in unstructured_grid_list:
            if is_nodes_used:
                list_of_ref_containers.append(mp.Nodes)
                list_of_ref_containers.append(mp.GetCommunicator().LocalMesh().Nodes)

            if cells_container is not None:
                list_of_ref_containers.append(cells_container)
                if isinstance(cells_container, Kratos.ConditionsArray):
                    list_of_ref_containers.append(mp.GetCommunicator().LocalMesh().Conditions)
                elif isinstance(cells_container, Kratos.ElementsArray):
                    list_of_ref_containers.append(mp.GetCommunicator().LocalMesh().Elements)

        self.assertEqual(list_of_containers, list_of_ref_containers)

    def test_PointVariableAddition(self):
        vtu_output = Kratos.Future.VtuOutput(self.model_part, output_sub_model_parts=True, output_format=Kratos.Future.VtuOutput.BINARY)
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

        with self.assertRaises(RuntimeError):
            vtu_output.AddTensorAdaptor("PRESSURE_1", ta)
        vtu_output.AddTensorAdaptor("PRESSURE_2", ta)

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

        vtu_output.PrintOutput("temp/vtu_output/variable_test")

        model_part = vtu_output.GetModelPart()
        unstructured_grid_list = self.GetUnstructuredGridList(model_part, model_part.GetCommunicator().GetDataCommunicator(), True)
        list_of_vtu_file_names: 'list[str]' = []
        for model_part, use_nodes, container in unstructured_grid_list:
            vtu_file_name = TestVtuOutput.GetUnstructuredGridName((model_part, use_nodes, container), "temp/vtu_output/variable_test", 0, model_part.GetCommunicator().GetDataCommunicator())
            list_of_vtu_file_names.append(vtu_file_name)
            if model_part.FullName() == "test":
                if use_nodes:
                    TestVtuOutput.CheckVtuFile(self, vtu_file_name, model_part.NumberOfNodes(), len(container), "binary",
                                            {
                                                "PRESSURE": (1, "Float64"),
                                                "PRESSURE_1": (5, "Float64"),
                                                "PRESSURE_2": (5, "Float64")
                                            },
                                            {
                                                "PRESSURE": (1, "Float64")
                                            })
                elif isinstance(container, Kratos.ConditionsArray):
                    TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNumberOfNodes(container), len(container), "binary",
                                            {},
                                            {
                                                "PRESSURE": (1, "Float64"),
                                                "PRESSURE_1": (8, "Float64"),
                                                "PRESSURE_2": (8, "Float64"),
                                                "PRESSURE_3": (8, "Float64")
                                            })
                elif isinstance(container, Kratos.ElementsArray):
                    TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNumberOfNodes(container), len(container), "binary",
                                            {
                                                "PRESSURE": (1, "Float64"),
                                                "PRESSURE_1": (5, "Float64"),
                                                "PRESSURE_2": (1, "Float64"),
                                                "PRESSURE_3": (5, "Float64")
                                            },
                                            {
                                                "PRESSURE": (1, "Float64"),
                                                "PRESSURE_1": (10, "Float64"),
                                                "PRESSURE_2": (10, "Float64"),
                                                "PRESSURE_3": (10, "Float64")
                                            })
            else:
                if use_nodes:
                    TestVtuOutput.CheckVtuFile(self, vtu_file_name, model_part.NumberOfNodes(), len(container), "binary",
                                            {
                                                "PRESSURE": (1, "Float64")
                                            },
                                            {
                                                "PRESSURE": (1, "Float64")
                                            })
                else:
                   TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNumberOfNodes(container), len(container), "binary",
                                            {
                                            },
                                            {
                                                "PRESSURE": (1, "Float64")
                                            })
            # check gauss
            list_of_vtu_file_names.append(TestVtuOutput.CheckGaussVtuFile(self, model_part, container, "temp/vtu_output/variable_test", 0, "binary", model_part.GetCommunicator().GetDataCommunicator(),
                                        {
                                            "PRESSURE": (1, "Float64")
                                        }))
        TestVtuOutput.CheckPvdFile(self, "temp/vtu_output/variable_test.pvd", list_of_vtu_file_names, [1 + 1e-9])

    def test_CellVariableAddition(self):
        vtu_output = Kratos.Future.VtuOutput(self.model_part, output_format=Kratos.Future.VtuOutput.ASCII, output_sub_model_parts=True)
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
        vtu_output.GetModelPart().ProcessInfo[Kratos.TIME] += 1e-9
        vtu_output.GetModelPart().ProcessInfo[Kratos.STEP] += 1
        vtu_output.PrintOutput("temp/vtu_output/time_step_test")

        model_part = vtu_output.GetModelPart()
        unstructured_grid_list = self.GetUnstructuredGridList(model_part, model_part.GetCommunicator().GetDataCommunicator(), True)
        list_of_vtu_file_names: 'list[str]' = []
        for step in range(2):
            for model_part, use_nodes, container in unstructured_grid_list:
                vtu_file_name = TestVtuOutput.GetUnstructuredGridName((model_part, use_nodes, container), "temp/vtu_output/time_step_test", step, model_part.GetCommunicator().GetDataCommunicator())
                list_of_vtu_file_names.append(vtu_file_name)
                if use_nodes:
                    TestVtuOutput.CheckVtuFile(self, vtu_file_name, model_part.NumberOfNodes(), len(container), "ascii",
                                            {
                                                "PRESSURE": (1, "Float64"),
                                                "DISPLACEMENT": (3, "Float64")
                                            },
                                            {
                                                "PRESSURE": (1, "Float64"),
                                                "DISPLACEMENT": (3, "Float64")
                                            })
                else:
                    TestVtuOutput.CheckVtuFile(self, vtu_file_name, TestVtuOutput.GetNumberOfNodes(container), len(container), "ascii",
                                            {
                                            },
                                            {
                                                "PRESSURE": (1, "Float64"),
                                                "DISPLACEMENT": (3, "Float64")
                                            })

                # check gauss
                list_of_vtu_file_names.append(TestVtuOutput.CheckGaussVtuFile(self, model_part, container, "temp/vtu_output/time_step_test", step, "ascii", model_part.GetCommunicator().GetDataCommunicator(),
                                            {
                                                "DISPLACEMENT": (3, "Float64")
                                            }))

        TestVtuOutput.CheckPvdFile(self, "temp/vtu_output/time_step_test.pvd", list_of_vtu_file_names, [1, 1 + 1e-9])

    @staticmethod
    def GetUnstructuredGridList(model_part: Kratos.ModelPart, data_communicator: Kratos.DataCommunicator, recursively: bool) -> 'list[tuple[Kratos.ModelPart, bool, typing.Optional[typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]]]]':
        unstructured_grid_list: 'list[tuple[Kratos.ModelPart, bool, typing.Optional[typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]]]]' = []
        TestVtuOutput.__GetUnstructuredGridList(unstructured_grid_list, model_part, data_communicator, recursively)
        return sorted(unstructured_grid_list, key = lambda x: x[0].FullName())

    @staticmethod
    def GetUnstructuredGridName(unstructured_grid: 'tuple[Kratos.ModelPart, bool, typing.Optional[typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]]]', prefix: str, step: int, data_communicator: Kratos.DataCommunicator, entity_suffix = "s") -> str:
        if data_communicator.IsDistributed():
            rank_suffix = f"_{data_communicator.Rank()}"
        else:
            rank_suffix = ""

        if unstructured_grid[2] is None:
            entity_type = "node"
        elif isinstance(unstructured_grid[2], Kratos.ConditionsArray):
            entity_type = "condition"
        elif isinstance(unstructured_grid[2], Kratos.ElementsArray):
            entity_type = "element"
        return f"{prefix}/{unstructured_grid[0].FullName()}_{entity_type}{entity_suffix}_{step}{rank_suffix}.vtu"

    @staticmethod
    def GetNodes(model_part: Kratos.ModelPart, use_model_part_nodes: bool, container: 'typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]') -> 'list[Kratos.Node]':
        if use_model_part_nodes:
            return model_part.Nodes
        else:
            temp = []
            for entity in container:
                for node in entity.GetGeometry():
                    temp.append(node.Id)
            return [model_part.GetRootModelPart().GetNode(node_id) for node_id in list(sorted(set(temp)))]

    def GetNumberOfNodes(container: 'typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]') -> int:
        temp = []
        for entity in container:
            for node in entity.GetGeometry():
                temp.append(node.Id)
        return len(list(set(temp)))

    @staticmethod
    def CheckDataArray(test_class: kratos_unittest.TestCase, xml_element: ET.Element, number_of_components: int, name: str, data_type: str, data_format: 'typing.Optional[str]' = None):
        test_class.assertEqual(xml_element.get("NumberOfComponents"), str(number_of_components))
        test_class.assertEqual(xml_element.get("type"), data_type)
        test_class.assertEqual(xml_element.get("Name"), name)
        if data_format is not None:
            test_class.assertEqual(xml_element.get("format"), data_format)

    @staticmethod
    def CheckVtuFile(
        test_class: kratos_unittest.TestCase,
        vtu_file_name: str,
        number_of_points: int,
        number_of_cells: 'typing.Optional[int]',
        data_format : str,
        point_data_fields: 'dict[str, tuple[int, str]]',
        cell_data_fields: 'dict[str, tuple[int, str]]'):

        test_class.assertTrue(Path(vtu_file_name).is_file(), f"The file {vtu_file_name} not found.")
        tree = ET.parse(vtu_file_name)
        root = tree.getroot()

        test_class.assertEqual(root.tag, "VTKFile")
        test_class.assertEqual(root.get("type"), "UnstructuredGrid")
        test_class.assertEqual(root.get("version"), "0.1")

        unstructured_grid = root.find("UnstructuredGrid")
        piece = unstructured_grid.find("Piece")

        if number_of_cells is None:
            local_number_of_cells = 0
        else:
            local_number_of_cells = number_of_cells

        test_class.assertEqual(piece.get("NumberOfPoints"), f"{number_of_points}", f"Vtu file name = {vtu_file_name}")
        test_class.assertEqual(piece.get("NumberOfCells"), f"{local_number_of_cells}", f"Vtu file name = {vtu_file_name}")

        points = piece.find("Points")
        TestVtuOutput.CheckDataArray(test_class, points.find("DataArray"), 3, "Position", "Float64", data_format)

        cells = piece.find("Cells")
        if not number_of_cells is None:
            data_arrays = cells.findall("DataArray")
            TestVtuOutput.CheckDataArray(test_class, data_arrays[0], 1, "connectivity", "Int32", data_format)
            TestVtuOutput.CheckDataArray(test_class, data_arrays[1], 1, "offsets", "Int32", data_format)
            TestVtuOutput.CheckDataArray(test_class, data_arrays[2], 1, "types", "UInt8", data_format)

        # now checking for point data
        point_data = piece.find("PointData")
        for data_field_name, (data_field_number_of_components, data_field_data_type) in point_data_fields.items():
            found = False
            for point_data_array in point_data.findall("DataArray"):
                if point_data_array.get("Name") == data_field_name:
                    TestVtuOutput.CheckDataArray(test_class, point_data_array, data_field_number_of_components, data_field_name, data_field_data_type, data_format)
                    found = True
                    break
            test_class.assertTrue(found, f"Point data field \"{data_field_name}\" not found in the \"{vtu_file_name}\"")

        if not number_of_cells is None:
            # now checking for cell data
            cell_data = piece.find("CellData")
            for data_field_name, (data_field_number_of_components, data_field_data_type) in cell_data_fields.items():
                found = False
                for point_data_array in cell_data.findall("DataArray"):
                    if point_data_array.get("Name") == data_field_name:
                        TestVtuOutput.CheckDataArray(test_class, point_data_array, data_field_number_of_components, data_field_name, data_field_data_type, data_format)
                        found = True
                        break
                test_class.assertTrue(found, f"Cell data field \"{data_field_name}\" not found in the \"{vtu_file_name}\"")
        kratos_utils.DeleteFileIfExisting(vtu_file_name)

    @staticmethod
    def CheckPvdFile(test_class: kratos_unittest.TestCase, pvd_file_name: str, vtu_file_name_list: 'list[str]', time_step_list: 'list[float]'):
        test_class.assertTrue(Path(pvd_file_name).is_file(), f"The file {pvd_file_name} not found.")
        tree = ET.parse(pvd_file_name)
        root = tree.getroot()

        test_class.assertEqual(root.tag, "VTKFile")
        test_class.assertEqual(root.get("type"), "Collection")
        test_class.assertEqual(root.get("version"), "1.0")

        collection = root.find("Collection")
        datasets = collection.findall("DataSet")
        test_class.assertEqual(len(datasets), len(vtu_file_name_list), f"file name = {pvd_file_name}, list_of_time_steps = {time_step_list}, list_of_vtu_files = \n" + "\n\t".join(vtu_file_name_list))
        for i, dataset in enumerate(datasets):
            relative_path = Path(vtu_file_name_list[i]).absolute().relative_to(Path(pvd_file_name).absolute().parent)
            test_class.assertEqual(dataset.get("file"), str(relative_path))
            test_class.assertEqual(dataset.get("name"), relative_path.name[:relative_path.name.rfind("_")])
            test_class.assertEqual(dataset.get("part"), str(i % (len(vtu_file_name_list) // len(time_step_list))))
            test_class.assertEqual(dataset.get("timestep"), f"{time_step_list[i // (len(vtu_file_name_list) // len(time_step_list))]:0.9e}", f"file name = {relative_path}")
        kratos_utils.DeleteFileIfExisting(pvd_file_name)

    @staticmethod
    def CheckGaussVtuFile(test_class: kratos_unittest.TestCase, model_part: Kratos.ModelPart, container, prefix: str, step_id: int, output_type: str, data_communicator: Kratos.DataCommunicator, point_fields):
        number_of_points = 0
        for entity in container:
            if len(entity.GetGeometry()) == 4:
                number_of_points += 4
            elif len(entity.GetGeometry()) == 3:
                number_of_points += 1
            elif len(entity.GetGeometry()) == 2:
                number_of_points += 1
            else:
                raise RuntimeError("Unsupported geometry")

        if isinstance(container, Kratos.ConditionsArray):
            suffix = "condition"
        elif isinstance(container, Kratos.ElementsArray):
            suffix = "element"
        else:
            raise RuntimeError("Unsupported container type.")

        if data_communicator.IsDistributed():
            rank_suffix = f"_{data_communicator.Rank()}"
        else:
            rank_suffix = ""

        file_name = f"{prefix}/{model_part.FullName()}_{suffix}_gauss_{step_id}{rank_suffix}.vtu"
        TestVtuOutput.CheckVtuFile(
            test_class,
            file_name,
            number_of_points, None, output_type, point_fields, {})
        kratos_utils.DeleteFileIfExisting(file_name)
        return file_name

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

    @classmethod
    def tearDownClass(cls):
        kratos_utils.DeleteDirectoryIfExistingAndEmpty("temp")


if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.INFO)
    kratos_unittest.main()