import math
import KratosMultiphysics as Kratos
import KratosMultiphysics.mpi as KratosMPI
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

from pathlib import Path
class TestVtuOutput(kratos_unittest.TestCase):
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
            partition_index = TestVtuOutput.__CreateEntities(model_part, my_pid, num_proc - 1)

        model_part.SetBufferSize(1)
        TestVtuOutput.__SetAllVariables(model_part)

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
        TestVtuOutput.__SetVariable(model_part.Nodes, Kratos.DISPLACEMENT, lambda x, y, z : x.SetSolutionStepValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Nodes, Kratos.PRESSURE, lambda x, y, z : x.SetSolutionStepValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Nodes, Kratos.ACTIVATION_LEVEL, lambda x, y, z : x.SetSolutionStepValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Nodes, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, lambda x, y, z : x.SetSolutionStepValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Nodes, Kratos.CONSTITUTIVE_MATRIX, lambda x, y, z : x.SetSolutionStepValue(y, z), offset)

        TestVtuOutput.__SetVariable(model_part.Nodes, Kratos.DISPLACEMENT, lambda x, y, z : x.SetValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Nodes, Kratos.PRESSURE, lambda x, y, z : x.SetValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Nodes, Kratos.ACTIVATION_LEVEL, lambda x, y, z : x.SetValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Nodes, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, lambda x, y, z : x.SetValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Nodes, Kratos.CONSTITUTIVE_MATRIX, lambda x, y, z : x.SetValue(y, z), offset)

        TestVtuOutput.__SetVariable(model_part.Conditions, Kratos.DISPLACEMENT, lambda x, y, z : x.SetValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Conditions, Kratos.PRESSURE, lambda x, y, z : x.SetValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Conditions, Kratos.ACTIVATION_LEVEL, lambda x, y, z : x.SetValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Conditions, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, lambda x, y, z : x.SetValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Conditions, Kratos.CONSTITUTIVE_MATRIX, lambda x, y, z : x.SetValue(y, z), offset)

        TestVtuOutput.__SetVariable(model_part.Elements, Kratos.DISPLACEMENT, lambda x, y, z : x.SetValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Elements, Kratos.PRESSURE, lambda x, y, z : x.SetValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Elements, Kratos.ACTIVATION_LEVEL, lambda x, y, z : x.SetValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Elements, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, lambda x, y, z : x.SetValue(y, z), offset)
        TestVtuOutput.__SetVariable(model_part.Elements, Kratos.CONSTITUTIVE_MATRIX, lambda x, y, z : x.SetValue(y, z), offset)

    def __CheckFile(self, file_name: str):
        data_communicator: Kratos.DataCommunicator = self.model_part.GetCommunicator().GetDataCommunicator()

        if data_communicator.Size() == 2 and file_name.rfind("ascii") != -1:
            ## Settings string in json format
            params = Kratos.Parameters("""{
                "reference_file_name" : "",
                "output_file_name"    : "",
                "comparison_type"     : "deterministic"
            }""")
            params["reference_file_name"].SetString(str(Path(f"auxiliar_files_for_python_unittest/reference_files/{file_name}")))
            params["output_file_name"].SetString(file_name)
            CompareTwoFilesCheckProcess(params).Execute()
        else:
            if data_communicator.Rank() == 0:
                self.assertTrue(Path(file_name).is_file())
                kratos_utils.DeleteFileIfExisting(file_name)

    def __OutputTest(self, output_type: str):
        if output_type == "ascii":
            vtu_output = Kratos.VtuOutput(self.model_part, binary_output=Kratos.VtuOutput.ASCII, output_sub_model_parts=True, echo_level=0)
        elif output_type == "binary":
            vtu_output = Kratos.VtuOutput(self.model_part, binary_output=Kratos.VtuOutput.BINARY, output_sub_model_parts=True, echo_level=0)

        for data_location in [Kratos.Globals.DataLocation.NodeNonHistorical, Kratos.Globals.DataLocation.Condition, Kratos.Globals.DataLocation.Element]:
            vtu_output.AddVariable(Kratos.DISPLACEMENT, data_location)
            vtu_output.AddVariable(Kratos.PRESSURE, data_location)
            vtu_output.AddVariable(Kratos.ACTIVATION_LEVEL, data_location)
            vtu_output.AddVariable(Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, data_location)
            vtu_output.AddVariable(Kratos.CONSTITUTIVE_MATRIX, data_location)

        for i in range(2):
            self.model_part.ProcessInfo[Kratos.TIME] += i
            self.__SetAllVariables(self.model_part, i)
            def AddExpression(exp_name_prefix, exp_type, variable, *args):
                c_exp = exp_type(self.model_part)
                Kratos.Expression.VariableExpressionIO.Read(c_exp, variable, *args)
                if (i == 0):
                    vtu_output.AddContainerExpression(f"exp_{exp_name_prefix}_{variable.Name()}", c_exp)
                else:
                    vtu_output.UpdateContainerExpression(f"exp_{exp_name_prefix}_{variable.Name()}", c_exp)

            AddExpression("hist", Kratos.Expression.NodalExpression, Kratos.PRESSURE, True)
            AddExpression("hist", Kratos.Expression.NodalExpression, Kratos.DISPLACEMENT, True)
            AddExpression("hist", Kratos.Expression.NodalExpression, Kratos.ACTIVATION_LEVEL, True)
            AddExpression("hist", Kratos.Expression.NodalExpression, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, True)
            AddExpression("hist", Kratos.Expression.NodalExpression, Kratos.CONSTITUTIVE_MATRIX, True)
            AddExpression("non_hist", Kratos.Expression.NodalExpression, Kratos.PRESSURE, False)
            AddExpression("non_hist", Kratos.Expression.NodalExpression, Kratos.DISPLACEMENT, False)
            AddExpression("non_hist", Kratos.Expression.NodalExpression, Kratos.ACTIVATION_LEVEL, False)
            AddExpression("non_hist", Kratos.Expression.NodalExpression, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT, False)
            AddExpression("non_hist", Kratos.Expression.NodalExpression, Kratos.CONSTITUTIVE_MATRIX, False)

            AddExpression("var", Kratos.Expression.ConditionExpression, Kratos.PRESSURE)
            AddExpression("var", Kratos.Expression.ConditionExpression, Kratos.DISPLACEMENT)
            AddExpression("var", Kratos.Expression.ConditionExpression, Kratos.ACTIVATION_LEVEL)
            AddExpression("var", Kratos.Expression.ConditionExpression, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
            AddExpression("var", Kratos.Expression.ConditionExpression, Kratos.CONSTITUTIVE_MATRIX)

            AddExpression("var", Kratos.Expression.ElementExpression, Kratos.PRESSURE)
            AddExpression("var", Kratos.Expression.ElementExpression, Kratos.DISPLACEMENT)
            AddExpression("var", Kratos.Expression.ElementExpression, Kratos.ACTIVATION_LEVEL)
            AddExpression("var", Kratos.Expression.ElementExpression, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
            AddExpression("var", Kratos.Expression.ElementExpression, Kratos.CONSTITUTIVE_MATRIX)

            def AddTensorAdaptor(ta_name_prefix, ta_type, container, variable):
                ta = ta_type(container, variable)
                ta.Check()
                ta.CollectData()
                if i == 0:
                    vtu_output.AddTensorAdaptor(f"ta_{ta_name_prefix}_{variable.Name()}", ta)
                else:
                    vtu_output.UpdateTensorAdaptor(f"ta_{ta_name_prefix}_{variable.Name()}", ta)

            AddTensorAdaptor("hist", Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, self.model_part.Nodes, Kratos.PRESSURE)
            AddTensorAdaptor("hist", Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, self.model_part.Nodes, Kratos.DISPLACEMENT)
            AddTensorAdaptor("hist", Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, self.model_part.Nodes, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
            AddTensorAdaptor("hist", Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, self.model_part.Nodes, Kratos.CONSTITUTIVE_MATRIX)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Nodes, Kratos.PRESSURE)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Nodes, Kratos.DISPLACEMENT)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Nodes, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
            AddTensorAdaptor("non_hist", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.Nodes, Kratos.CONSTITUTIVE_MATRIX)

            AddTensorAdaptor("hist_local", Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Nodes, Kratos.PRESSURE)
            AddTensorAdaptor("hist_local", Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Nodes, Kratos.DISPLACEMENT)
            AddTensorAdaptor("hist_local", Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Nodes, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
            AddTensorAdaptor("hist_local", Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Nodes, Kratos.CONSTITUTIVE_MATRIX)
            AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Nodes, Kratos.PRESSURE)
            AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Nodes, Kratos.DISPLACEMENT)
            AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Nodes, Kratos.DETERMINANTS_OF_JACOBIAN_PARENT)
            AddTensorAdaptor("non_hist_local", Kratos.TensorAdaptors.VariableTensorAdaptor, self.model_part.GetCommunicator().LocalMesh().Nodes, Kratos.CONSTITUTIVE_MATRIX)

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

            vtu_output.PrintOutput(f"vtu_output/{output_type}_output")

        # now check the files
        data_communicator: Kratos.DataCommunicator = self.model_part.GetCommunicator().GetDataCommunicator()
        n_procs = data_communicator.Size()

        for step_id in range(2):
            for proc_id in range(n_procs):
                self.__CheckFile(f"vtu_output/{output_type}_output/test_conditions_{step_id}_{proc_id}.vtu")
                self.__CheckFile(f"vtu_output/{output_type}_output/test_elements_{step_id}_{proc_id}.vtu")
                self.__CheckFile(f"vtu_output/{output_type}_output/test.sub_1_conditions_{step_id}_{proc_id}.vtu")
                self.__CheckFile(f"vtu_output/{output_type}_output/test.sub_1_elements_{step_id}_{proc_id}.vtu")
                self.__CheckFile(f"vtu_output/{output_type}_output/test.sub_2_conditions_{step_id}_{proc_id}.vtu")
                self.__CheckFile(f"vtu_output/{output_type}_output/test.sub_2.sub_1_conditions_{step_id}_{proc_id}.vtu")

            self.__CheckFile(f"vtu_output/{output_type}_output/test_conditions_{step_id}.pvtu")
            self.__CheckFile(f"vtu_output/{output_type}_output/test_elements_{step_id}.pvtu")
            self.__CheckFile(f"vtu_output/{output_type}_output/test.sub_1_conditions_{step_id}.pvtu")
            self.__CheckFile(f"vtu_output/{output_type}_output/test.sub_1_elements_{step_id}.pvtu")
            self.__CheckFile(f"vtu_output/{output_type}_output/test.sub_2_conditions_{step_id}.pvtu")
            self.__CheckFile(f"vtu_output/{output_type}_output/test.sub_2.sub_1_conditions_{step_id}.pvtu")

        # checking the collection file
        self.__CheckFile(f"vtu_output/{output_type}_output.pvd")

    def test_OutputASCII(self):
        self.__OutputTest("ascii")

    def test_OutputBinary(self):
        self.__OutputTest("binary")

    @classmethod
    def tearDownClass(cls):
        cls.model_part.GetCommunicator().GetDataCommunicator().Barrier()
        kratos_utils.DeleteDirectoryIfExisting("vtu_output")

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.INFO)
    kratos_unittest.main()