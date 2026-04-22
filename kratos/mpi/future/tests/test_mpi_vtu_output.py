import math
import KratosMultiphysics as Kratos
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.kratos_utilities as kratos_utils
import xml.etree.ElementTree as ET

with kratos_unittest.WorkFolderScope("../../../future/tests", __file__, True):
    import test_vtu_output

from pathlib import Path
class TestMPIVtuOutput(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(cls.model_part, Kratos.Testing.GetDefaultDataCommunicator())
        cls.__InitializeModelPart(cls.model_part)

    @staticmethod
    def __CreateEntities(model_part: Kratos.ModelPart, my_pid: int, num_proc: int) -> None:
        my_num_quad = 2 # user-defined.
        num_local_nodes = 3 * my_num_quad
        num_ghost_nodes = 3
        local_start_index = num_local_nodes * my_pid + 1
        ghost_start_index = local_start_index + num_local_nodes
        local_node_ids = list(range(local_start_index, local_start_index + num_local_nodes))
        ghost_node_ids = list(range(ghost_start_index, ghost_start_index + num_ghost_nodes))

        partition_index = dict()
        for i in local_node_ids:
            partition_index[i] = my_pid
        if (my_pid == num_proc - 1): # Connect ring start and ring end.
            ghost_node_ids = [1, 2, 3]
        for i in ghost_node_ids:
            partition_index[i] = (my_pid + 1) % num_proc
        node_ids = local_node_ids + ghost_node_ids
        # Create nodes.
        for i in node_ids:
            radius = 0.5 + 0.5 * ((i - 1) % 3) / 2.0
            phase = 2.0 * math.pi * ((i - 1) // 3) / float(my_num_quad * num_proc)
            x = radius * math.cos(phase)
            y = radius * math.sin(phase)
            model_part.CreateNewNode(i, x, y, 0.0)
        # Create elements and conditions.
        for i in range(0, num_local_nodes, 3):
            prop_id = 1
            prop = model_part.GetProperties()[prop_id]
            # First triangle.
            eid = local_start_index + 3 * (i // 3)
            nids = [node_ids[i], node_ids[i + 1], node_ids[i + 4]]
            model_part.CreateNewElement("Element2D3N", eid, nids, prop)
            model_part.CreateNewCondition("SurfaceCondition3D3N", eid, nids, prop)
            # Second triangle.
            eid = eid + 1
            nids = [node_ids[i], node_ids[i + 4], node_ids[i + 3]]
            model_part.CreateNewElement("Element2D3N", eid, nids, prop)
            model_part.CreateNewCondition("SurfaceCondition3D3N", eid, nids, prop)
            # Quad.
            eid = eid + 1
            nids = [node_ids[i + 1], node_ids[i + 2], node_ids[i + 5], node_ids[i + 4]]
            model_part.CreateNewElement("Element2D4N", eid, nids, prop)
            model_part.CreateNewCondition("SurfaceCondition3D4N", eid, nids, prop)
        if my_pid == 0:
            # Here we create a special condition that only exists on the first
            # process. This is to test the collective write when at least one
            # process has an empty set.
            model_part.CreateNewCondition("LineCondition2D2N", eid + 1, [node_ids[i + 1], node_ids[i + 2]], prop)

        return partition_index

    @staticmethod
    def __InitializeModelPart(model_part: Kratos.ModelPart):
        # Add variables.
        model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT) # array_1d
        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE) # double
        model_part.AddNodalSolutionStepVariable(Kratos.ACTIVATION_LEVEL) # int
        model_part.AddNodalSolutionStepVariable(Kratos.DETERMINANTS_OF_JACOBIAN_PARENT) # vector
        model_part.AddNodalSolutionStepVariable(Kratos.CONSTITUTIVE_MATRIX) # matrix
        model_part.AddNodalSolutionStepVariable(Kratos.PARTITION_INDEX)

        # Make a mesh out of two structured rings (inner triangles, outer quads).
        data_communicator: Kratos.DataCommunicator = model_part.GetCommunicator().GetDataCommunicator()
        num_proc = data_communicator.Size()
        my_pid = data_communicator.Rank()

        # the last rank is kept empty for empty rank check
        if my_pid != num_proc - 1:
            partition_index = TestMPIVtuOutput.__CreateEntities(model_part, my_pid, num_proc - 1)

        model_part.SetBufferSize(1)
        TestMPIVtuOutput.__SetAllVariables(model_part)

        # Write PARTITION_INDEX
        for node in model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PARTITION_INDEX, partition_index[node.Id])

        KratosMPI.ParallelFillCommunicator(model_part.GetRootModelPart(), data_communicator).Execute()
        communicator: Kratos.Communicator = model_part.GetCommunicator()
        communicator.SynchronizeNodalSolutionStepsData()
        communicator.SynchronizeNonHistoricalVariable(Kratos.DISPLACEMENT)
        communicator.SynchronizeNonHistoricalVariable(Kratos.PRESSURE)
        communicator.SynchronizeNonHistoricalVariable(Kratos.ACTIVATION_LEVEL)
        communicator.SynchronizeNonHistoricalVariable(Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
        communicator.SynchronizeNonHistoricalVariable(Kratos.CONSTITUTIVE_MATRIX)
        communicator.SynchronizeNodalFlags()

        # create some sub_model_parts
        sub_model_part = model_part.CreateSubModelPart("sub_1")
        for condition in model_part.Conditions:
            if condition.Id % 2 == 0:
                sub_model_part.AddCondition(condition)
        for element in model_part.Elements:
            if element.Id % 2 == 1:
                sub_model_part.AddElement(element)

        # create some sub_model_parts
        sub_model_part = model_part.CreateSubModelPart("sub_2.sub_1")
        for condition in model_part.Conditions:
            if condition.Id % 3 == 2:
                sub_model_part.AddCondition(condition)

        node_ids: 'list[int]' = []
        for condition in sub_model_part.Conditions:
            for node in condition.GetGeometry():
                node_ids.append(node.Id)
        sub_model_part.AddNodes(node_ids)

        # Set some process info variables.
        model_part.ProcessInfo[Kratos.TIME] = 1.2345 # float

    @staticmethod
    def __SetVariable(container, variable, setter, offset = 0):
        for entity in container:
            entity_id = entity.Id + offset
            if isinstance(variable, Kratos.IntegerVariable):
                setter(entity, variable, entity_id + 1)
            elif isinstance(variable, Kratos.DoubleVariable):
                setter(entity, variable, entity_id + 1)
            elif isinstance(variable, Kratos.Array1DVariable3):
                setter(entity, variable, [entity_id + 1, entity_id + 2, entity_id + 3])
            elif isinstance(variable, Kratos.Array1DVariable4):
                setter(entity, variable, [entity_id + 1, entity_id + 2, entity_id + 3, entity_id + 4])
            elif isinstance(variable, Kratos.Array1DVariable6):
                setter(entity, variable, [entity_id + 1, entity_id + 2, entity_id + 3, entity_id + 4, entity_id + 5])
            elif isinstance(variable, Kratos.Array1DVariable9):
                setter(entity, variable, [entity_id + 1, entity_id + 2, entity_id + 3, entity_id + 4, entity_id + 5, entity_id + 6])
            elif isinstance(variable, Kratos.VectorVariable):
                setter(entity, variable, [entity_id + 1, entity_id + 2, entity_id + 3, entity_id + 4, entity_id + 5, entity_id + 6, entity_id + 7])
            elif isinstance(variable, Kratos.MatrixVariable):
                setter(entity, variable, Kratos.Matrix([[entity_id + 1, entity_id + 2], [entity_id + 3, entity_id + 4], [entity.Id + 5, entity_id + 6]]))

    @staticmethod
    def __SetAllVariables(model_part: Kratos.ModelPart, offset = 0):
        TestMPIVtuOutput.__SetVariable(model_part.Nodes, Kratos.DISPLACEMENT, lambda x, y, z : x.SetSolutionStepValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Nodes, Kratos.PRESSURE, lambda x, y, z : x.SetSolutionStepValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Nodes, Kratos.ACTIVATION_LEVEL, lambda x, y, z : x.SetSolutionStepValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Nodes, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, lambda x, y, z : x.SetSolutionStepValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Nodes, Kratos.CONSTITUTIVE_MATRIX, lambda x, y, z : x.SetSolutionStepValue(y, z), offset)

        TestMPIVtuOutput.__SetVariable(model_part.Nodes, Kratos.DISPLACEMENT, lambda x, y, z : x.SetValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Nodes, Kratos.PRESSURE, lambda x, y, z : x.SetValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Nodes, Kratos.ACTIVATION_LEVEL, lambda x, y, z : x.SetValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Nodes, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, lambda x, y, z : x.SetValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Nodes, Kratos.CONSTITUTIVE_MATRIX, lambda x, y, z : x.SetValue(y, z), offset)

        TestMPIVtuOutput.__SetVariable(model_part.Conditions, Kratos.DISPLACEMENT, lambda x, y, z : x.SetValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Conditions, Kratos.PRESSURE, lambda x, y, z : x.SetValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Conditions, Kratos.ACTIVATION_LEVEL, lambda x, y, z : x.SetValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Conditions, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, lambda x, y, z : x.SetValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Conditions, Kratos.CONSTITUTIVE_MATRIX, lambda x, y, z : x.SetValue(y, z), offset)

        TestMPIVtuOutput.__SetVariable(model_part.Elements, Kratos.DISPLACEMENT, lambda x, y, z : x.SetValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Elements, Kratos.PRESSURE, lambda x, y, z : x.SetValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Elements, Kratos.ACTIVATION_LEVEL, lambda x, y, z : x.SetValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Elements, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, lambda x, y, z : x.SetValue(y, z), offset)
        TestMPIVtuOutput.__SetVariable(model_part.Elements, Kratos.CONSTITUTIVE_MATRIX, lambda x, y, z : x.SetValue(y, z), offset)

    def setUp(self):
        # Set some process info variables.
        self.model_part.ProcessInfo[Kratos.STEP] = 0 # float
        self.model_part.ProcessInfo[Kratos.TIME] = 1.2345 # float

    def __OutputTest(self, output_type: str):
        if output_type == "ascii":
            vtu_output = Kratos.Future.VtuOutput(self.model_part, output_format=Kratos.Future.VtuOutput.ASCII, output_sub_model_parts=True, echo_level=0)
        elif output_type == "binary":
            vtu_output = Kratos.Future.VtuOutput(self.model_part, output_format=Kratos.Future.VtuOutput.BINARY, output_sub_model_parts=True, echo_level=0)

        for data_location in [Kratos.Globals.DataLocation.NodeNonHistorical, Kratos.Globals.DataLocation.Condition, Kratos.Globals.DataLocation.Element]:
            vtu_output.AddVariable(Kratos.DISPLACEMENT, data_location)
            vtu_output.AddVariable(Kratos.PRESSURE, data_location)
            vtu_output.AddVariable(Kratos.ACTIVATION_LEVEL, data_location)
            vtu_output.AddVariable(Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, data_location)
            vtu_output.AddVariable(Kratos.CONSTITUTIVE_MATRIX, data_location)

        for data_location in [Kratos.Globals.DataLocation.Condition, Kratos.Globals.DataLocation.Element]:
            vtu_output.AddIntegrationPointVariable(Kratos.DISPLACEMENT, data_location)
            vtu_output.AddIntegrationPointVariable(Kratos.PRESSURE, data_location)
            vtu_output.AddIntegrationPointVariable(Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, data_location)
            vtu_output.AddIntegrationPointVariable(Kratos.CONSTITUTIVE_MATRIX, data_location)

        for i in range(2):
            self.model_part.ProcessInfo[Kratos.STEP] = i
            self.model_part.ProcessInfo[Kratos.TIME] += i
            self.__SetAllVariables(self.model_part, i)

            def AddTensorAdaptor(ta_name_prefix, ta_type, container, variable):
                ta = ta_type(container, variable)
                ta.Check()
                ta.CollectData()
                if i == 0:
                    vtu_output.AddTensorAdaptor(f"ta_{ta_name_prefix}_{variable.Name()}", ta)
                else:
                    vtu_output.ReplaceTensorAdaptor(f"ta_{ta_name_prefix}_{variable.Name()}", ta)

            AddTensorAdaptor("hist", Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, self.model_part.Nodes, Kratos.PRESSURE)
            AddTensorAdaptor("hist", Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, self.model_part.Nodes, Kratos.DISPLACEMENT)
            AddTensorAdaptor("hist", Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, self.model_part.Nodes, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
            AddTensorAdaptor("hist", Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, self.model_part.Nodes, Kratos.CONSTITUTIVE_MATRIX)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Nodes, Kratos.PRESSURE)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Nodes, Kratos.DISPLACEMENT)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Nodes, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Nodes, Kratos.CONSTITUTIVE_MATRIX)

            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Conditions, Kratos.PRESSURE)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Conditions, Kratos.DISPLACEMENT)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Conditions, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Conditions, Kratos.CONSTITUTIVE_MATRIX)

            AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Conditions, Kratos.PRESSURE)
            AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Conditions, Kratos.DISPLACEMENT)
            AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Conditions, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
            AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Conditions, Kratos.CONSTITUTIVE_MATRIX)

            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Elements, Kratos.PRESSURE)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Elements, Kratos.DISPLACEMENT)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Elements, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Elements, Kratos.CONSTITUTIVE_MATRIX)

            AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Elements, Kratos.PRESSURE)
            AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Elements, Kratos.DISPLACEMENT)
            AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Elements, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
            AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Elements, Kratos.CONSTITUTIVE_MATRIX)

            with self.assertRaises(RuntimeError):
                AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Nodes, Kratos.PRESSURE)

            with self.assertRaises(RuntimeError):
                AddTensorAdaptor("non_hist_ghost", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().GhostMesh().Nodes, Kratos.PRESSURE)

            with self.assertRaises(RuntimeError):
                AddTensorAdaptor("non_hist_interface", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().InterfaceMesh().Nodes, Kratos.PRESSURE)

            with self.assertRaises(RuntimeError):
                AddTensorAdaptor("non_hist_ghost", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().GhostMesh().Conditions, Kratos.PRESSURE)

            with self.assertRaises(RuntimeError):
                AddTensorAdaptor("non_hist_interface", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().InterfaceMesh().Conditions, Kratos.PRESSURE)

            with self.assertRaises(RuntimeError):
                AddTensorAdaptor("non_hist_ghost", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().GhostMesh().Elements, Kratos.PRESSURE)

            with self.assertRaises(RuntimeError):
                AddTensorAdaptor("non_hist_interface", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().InterfaceMesh().Elements, Kratos.PRESSURE)

            vtu_output.PrintOutput(f"temp/vtu_output/{output_type}_output")

        # now check the files
        data_communicator: Kratos.DataCommunicator = self.model_part.GetCommunicator().GetDataCommunicator()
        proc_id = data_communicator.Rank()

        def check_pvtu(model_part: Kratos.ModelPart, suffix: str, step_id: int, point_fields, cell_fields):
            if data_communicator.Rank() == 0:
                pvtu_file_name = f"temp/vtu_output/{output_type}_output/{model_part.FullName()}_{suffix}_{step_id}.pvtu"

                tree = ET.parse(pvtu_file_name)
                root = tree.getroot()

                self.assertEqual(root.tag, "VTKFile")
                self.assertEqual(root.get("type"), "PUnstructuredGrid")
                self.assertEqual(root.get("version"), "0.1")

                unstructured_grid = root.find("PUnstructuredGrid")
                self.assertEqual(unstructured_grid.get("GhostLevel"), "0")

                ppoints = unstructured_grid.find("PPoints")
                test_vtu_output.TestVtuOutput.CheckDataArray(self, ppoints.find("PDataArray"), 3, "Position", "Float64")

                pcells = unstructured_grid.find("PCells")
                p_data_arrays = pcells.findall("PDataArray")
                test_vtu_output.TestVtuOutput.CheckDataArray(self, p_data_arrays[0], 1, "connectivity", "Int32")
                test_vtu_output.TestVtuOutput.CheckDataArray(self, p_data_arrays[1], 1, "offsets", "Int32")
                test_vtu_output.TestVtuOutput.CheckDataArray(self, p_data_arrays[2], 1, "types", "UInt8")

                ppoint_data = unstructured_grid.find("PPointData")
                for data_field_name, (number_of_components, data_type) in point_fields:
                    found_field = False
                    for p_data_array in ppoint_data.findall("PDataArray"):
                        if p_data_array.get("Name") == data_field_name:
                            found_field = True
                            test_vtu_output.TestVtuOutput.CheckDataArray(self, p_data_array, number_of_components, data_field_name, data_type)
                            break
                    self.assertTrue(found_field)

                pcell_data = unstructured_grid.find("PCellData")
                for data_field_name, (number_of_components, data_type) in cell_fields:
                    found_field = False
                    for p_data_array in pcell_data.findall("PDataArray"):
                        if p_data_array.get("Name") == data_field_name:
                            found_field = True
                            test_vtu_output.TestVtuOutput.CheckDataArray(self, p_data_array, number_of_components, data_field_name, data_type)
                            break
                    self.assertTrue(found_field)

                pieces = unstructured_grid.findall("Piece")
                self.assertEqual(len(pieces), data_communicator.Size())
                for proc_id in range(data_communicator.Size()):
                    vtu_file_name = f"temp/vtu_output/{output_type}_output/{model_part.FullName()}_{suffix}_{step_id}_{proc_id}.vtu"
                    relative_path = str(Path(vtu_file_name).absolute().relative_to(Path(pvtu_file_name).parent.absolute()))
                    self.assertEqual(pieces[proc_id].get("Source"), relative_path)

                kratos_utils.DeleteFileIfExisting(pvtu_file_name)
                return pvtu_file_name

        def check(model_part: Kratos.ModelPart, number_of_nodes: int, number_of_cells: int, step_id: int, suffix: str, point_fields, cell_fields):
            file_name = f"temp/vtu_output/{output_type}_output/{model_part.FullName()}_{suffix}_{step_id}_{proc_id}.vtu"
            test_vtu_output.TestVtuOutput.CheckVtuFile(
                self,
                file_name,
                number_of_nodes, number_of_cells, output_type, point_fields, cell_fields)
            kratos_utils.DeleteFileIfExisting(file_name)

            return check_pvtu(model_part, suffix, step_id, point_fields, cell_fields)

        def check_gauss(model_part: Kratos.ModelPart, container, step_id: int, point_fields):
            file_name = test_vtu_output.TestVtuOutput.CheckGaussVtuFile(
                self,
                model_part,
                container,
                f"temp/vtu_output/{output_type}_output",
                step_id,
                output_type,
                model_part.GetRootModelPart().GetCommunicator().GetDataCommunicator(),
                point_fields
            )

            if isinstance(container, Kratos.ConditionsArray):
                suffix = "condition"
            elif isinstance(container, Kratos.ElementsArray):
                suffix = "element"
            else:
                raise RuntimeError("Unsupported container type.")

            kratos_utils.DeleteFileIfExisting(file_name)
            return check_pvtu(model_part, suffix + "_gauss", step_id, point_fields, {})

        list_of_pvtu_file_names: 'list[str]' = []
        for step_id in range(2):
            list_of_pvtu_file_names.append(check(self.model["test"], self.model["test"].NumberOfNodes(), len(
                self.model["test"].Elements), step_id, "elements", {}, {}))
            list_of_pvtu_file_names.append(check_gauss(self.model["test"], self.model["test"].Elements, step_id, {}))

            list_of_pvtu_file_names.append(check(self.model["test"], test_vtu_output.TestVtuOutput.GetNumberOfNodes(
                self.model["test"].Conditions), len(self.model["test"].Conditions), step_id, "conditions", {}, {}))
            list_of_pvtu_file_names.append(check_gauss(self.model["test"], self.model["test"].Conditions, step_id, {}))

            list_of_pvtu_file_names.append(check(self.model["test.sub_1"], test_vtu_output.TestVtuOutput.GetNumberOfNodes(
                self.model["test.sub_1"].Elements), len(self.model["test.sub_1"].Elements), step_id, "elements", {}, {}))
            list_of_pvtu_file_names.append(check_gauss(self.model["test.sub_1"], self.model["test.sub_1"].Elements, step_id, {}))

            list_of_pvtu_file_names.append(check(self.model["test.sub_1"], test_vtu_output.TestVtuOutput.GetNumberOfNodes(
                self.model["test.sub_1"].Conditions), len(self.model["test.sub_1"].Conditions), step_id, "conditions", {}, {}))
            list_of_pvtu_file_names.append(check_gauss(self.model["test.sub_1"], self.model["test.sub_1"].Conditions, step_id, {}))

            list_of_pvtu_file_names.append(check(self.model["test.sub_2"], test_vtu_output.TestVtuOutput.GetNumberOfNodes(
                self.model["test.sub_2"].Conditions), len(self.model["test.sub_2"].Conditions), step_id, "conditions", {}, {}))
            list_of_pvtu_file_names.append(check_gauss(self.model["test.sub_2"], self.model["test.sub_2"].Conditions, step_id, {}))

            list_of_pvtu_file_names.append(check(self.model["test.sub_2.sub_1"], test_vtu_output.TestVtuOutput.GetNumberOfNodes(
                self.model["test.sub_2.sub_1"].Conditions), len(self.model["test.sub_2.sub_1"].Conditions), step_id, "conditions", {}, {}))
            list_of_pvtu_file_names.append(check_gauss(self.model["test.sub_2.sub_1"], self.model["test.sub_2.sub_1"].Conditions, step_id, {}))

        if data_communicator.Rank() == 0:
            test_vtu_output.TestVtuOutput.CheckPvdFile(self, f"temp/vtu_output/{output_type}_output.pvd", list_of_pvtu_file_names, [1.2345, 2.2345])

    def test_OutputASCII(self):
        self.__OutputTest("ascii")

    def test_OutputBinary(self):
        self.__OutputTest("binary")

    @classmethod
    def tearDownClass(cls):
        cls.model_part.GetCommunicator().GetDataCommunicator().Barrier()
        kratos_utils.DeleteDirectoryIfExistingAndEmpty("temp")

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.INFO)
    kratos_unittest.main()