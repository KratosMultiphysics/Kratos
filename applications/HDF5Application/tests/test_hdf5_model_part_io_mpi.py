import os
from KratosMultiphysics import *
import KratosMultiphysics.mpi as KratosMPI
from KratosMultiphysics.HDF5Application import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import random, math
import pathlib

class ControlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class TestCase(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def __create_entities(self, model_part: ModelPart, my_pid: int, num_proc: int) -> None:
        my_num_quad = 2 # user-defined.
        my_num_tri = 2 * my_num_quad # splits each quad into 2 triangles.
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

    def _initialize_model_part(self, model_part: ModelPart):
        # Add variables.
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT) # array_1d
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(ACCELERATION)
        model_part.AddNodalSolutionStepVariable(PRESSURE) # double
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(ACTIVATION_LEVEL) # int
        model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
        # Make a mesh out of two structured rings (inner triangles, outer quads).
        data_communicator = model_part.GetCommunicator().GetDataCommunicator()
        num_proc = data_communicator.Size()
        my_pid = data_communicator.Rank()

        # the last rank is kept empty for empty rank check
        if my_pid != num_proc - 1:
            partition_index = self.__create_entities(model_part, my_pid, num_proc - 1)

        model_part.SetBufferSize(2)
        # Write some data to the nodal solution steps variables.
        for node in model_part.Nodes:
            node.SetSolutionStepValue(PARTITION_INDEX, partition_index[node.Id])
            node.SetSolutionStepValue(DISPLACEMENT_X, random.random())
            node.SetSolutionStepValue(DISPLACEMENT_Y, random.random())
            node.SetSolutionStepValue(DISPLACEMENT_Z, random.random())
            node.SetSolutionStepValue(VELOCITY_X, random.random())
            node.SetSolutionStepValue(VELOCITY_Y, random.random())
            node.SetSolutionStepValue(VELOCITY_Z, random.random())
            node.SetSolutionStepValue(ACCELERATION_X, random.random())
            node.SetSolutionStepValue(ACCELERATION_Y, random.random())
            node.SetSolutionStepValue(ACCELERATION_Z, random.random())
            node.SetSolutionStepValue(PRESSURE, random.random())
            node.SetSolutionStepValue(VISCOSITY, random.random())
            node.SetSolutionStepValue(DENSITY, random.random())
            node.SetSolutionStepValue(ACTIVATION_LEVEL, random.randint(-100, 100))

            node.SetValue(DISPLACEMENT_X, random.random())
            node.SetValue(DISPLACEMENT_Y, random.random())
            node.SetValue(DISPLACEMENT_Z, random.random())
            node.SetValue(VELOCITY_X, random.random())
            node.SetValue(VELOCITY_Y, random.random())
            node.SetValue(VELOCITY_Z, random.random())
            node.SetValue(ACCELERATION_X, random.random())
            node.SetValue(ACCELERATION_Y, random.random())
            node.SetValue(ACCELERATION_Z, random.random())
            node.SetValue(PRESSURE, random.random())
            node.SetValue(VISCOSITY, random.random())
            node.SetValue(DENSITY, random.random())
            node.SetValue(ACTIVATION_LEVEL, random.randint(-100, 100))

            node.Set(SLIP, bool(random.randint(-100, 100) % 2))
            node.Set(ACTIVE, bool(random.randint(-100, 100) % 2))

        KratosMPI.ParallelFillCommunicator(model_part.GetRootModelPart(), data_communicator).Execute()
        model_part.GetCommunicator().SynchronizeNodalSolutionStepsData()

        model_part.GetCommunicator().SynchronizeNonHistoricalVariable(DISPLACEMENT)
        model_part.GetCommunicator().SynchronizeNonHistoricalVariable(VELOCITY)
        model_part.GetCommunicator().SynchronizeNonHistoricalVariable(ACCELERATION)
        model_part.GetCommunicator().SynchronizeNonHistoricalVariable(PRESSURE)
        model_part.GetCommunicator().SynchronizeNonHistoricalVariable(VISCOSITY)
        model_part.GetCommunicator().SynchronizeNonHistoricalVariable(DENSITY)
        model_part.GetCommunicator().SynchronizeNonHistoricalVariable(ACTIVATION_LEVEL)

        model_part.GetCommunicator().SynchronizeNodalFlags()

        # Set some process info variables.
        model_part.ProcessInfo[DOMAIN_SIZE] = 3 # int
        model_part.ProcessInfo[TIME] = 1.2345 # float
        initial_strain = Vector(6)
        for i in range(6):
            initial_strain[i] = math.cos(i)
        model_part.ProcessInfo[INITIAL_STRAIN] = initial_strain # vector
        gl_strain_tensor = Matrix(5,5)
        for i in range(5):
            for j in range(5):
                gl_strain_tensor[i,j] = math.cos(i + j)
        model_part.ProcessInfo[GREEN_LAGRANGE_STRAIN_TENSOR] = gl_strain_tensor # matrix

    def _get_file(self, communicator: DataCommunicator):
        params = Parameters("""
        {
            "file_name" : "test_hdf5_model_part_io_mpi.h5",
            "file_access_mode" : "truncate",
            "file_driver" : "mpio",
            "echo_level" : 0
        }""")
        return HDF5File(communicator, params)

    def _get_model_part_io(self, hdf5_file: HDF5File):
        return HDF5PartitionedModelPartIO(hdf5_file, "/ModelData")

    def _get_nodal_solution_step_data_io(self, hdf5_file: HDF5File):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["DISPLACEMENT", "VELOCITY", "ACCELERATION", "PRESSURE", "VISCOSITY", "DENSITY", "ACTIVATION_LEVEL"]
        }""")
        return HDF5NodalSolutionStepDataIO(params, hdf5_file)

    def _get_nodal_data_io(self, hdf5_file: HDF5File):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["DISPLACEMENT", "VELOCITY", "ACCELERATION", "PRESSURE", "VISCOSITY", "DENSITY", "ACTIVATION_LEVEL"]
        }""")
        return HDF5NodalDataValueIO(params, hdf5_file)

    def _get_nodal_flag_io(self, hdf5_file: HDF5File):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["SLIP", "ACTIVE"]
        }""")
        return HDF5NodalFlagValueIO(params, hdf5_file)

    def test_GetListOfAvailableVariables(self):
        with ControlledExecutionScope(pathlib.Path(__file__).absolute().parent):
            current_model = Model()
            model_part = current_model.CreateModelPart("write")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part, Testing.GetDefaultDataCommunicator())
            self._initialize_model_part(model_part)

            data_communicator: DataCommunicator = model_part.GetCommunicator().GetDataCommunicator()

            my_pid = data_communicator.Rank()
            element: Element
            for element in model_part.Elements:
                if my_pid % 2 == 0 and element.Id % 2 == 0:
                    element.SetValue(DENSITY, 1.0)
                elif my_pid % 2 == 1 and element.Id % 2 == 1:
                    element.SetValue(DISTANCE, 2.0)
                elif my_pid % 2 == 1 and element.Id % 2 == 0:
                    element.SetValue(RADIUS, 3.0)
                else:
                    element.SetValue(BOSSAK_ALPHA, 4.0)

            if data_communicator.Size() == 2:
                # if this test was run with 2 processes in mpi, the last rank will be
                # empty by design, hence there wont be any entities with DISTANCE and RADIUS.
                # so the following check is done. Otherwise, if we check for "DISTANCE" and RADIUS,
                # then an error saying there are no entities having those variables will be thrown.
                self.assertEqual(Utilities.GetListOfAvailableVariables(model_part.Elements, data_communicator), ['BOSSAK_ALPHA', 'DENSITY'])
            else:
                self.assertEqual(Utilities.GetListOfAvailableVariables(model_part.Elements, data_communicator), ['BOSSAK_ALPHA', 'DENSITY', 'DISTANCE', 'RADIUS'])

    def test_HDF5PropertiesIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            current_model = Model()
            write_model_part = current_model.CreateModelPart("write")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(write_model_part, Testing.GetDefaultDataCommunicator())

            number_of_properties = 100
            data_communicator: DataCommunicator = write_model_part.GetCommunicator().GetDataCommunicator()
            num_proc = data_communicator.Size()
            my_pid = data_communicator.Rank()

            # the last rank is kept empty for empty rank check
            for prop_index in range(number_of_properties):
                # last rank does not have any properties
                if (prop_index % (my_pid + 1) == 0 and my_pid != num_proc - 1):
                    props: Properties = write_model_part.CreateNewProperties(prop_index + 1)
                    if (prop_index % 2 == 0): props.SetValue(DISTANCE, prop_index)
                    if (prop_index % 3 == 0): props.SetValue(VELOCITY, Vector([prop_index, prop_index + 1, prop_index + 2]))

                    sub_props: Properties = write_model_part.CreateNewProperties(props.Id + number_of_properties)
                    if (prop_index % 4 == 0): sub_props.SetValue(DENSITY, prop_index * 2)
                    if (prop_index % 5 == 0): sub_props.SetValue(ACCELERATION, Vector([prop_index * 2, prop_index * 2 + 1, prop_index * 2 + 2]))
                    props.AddSubProperties(sub_props)

            hdf5_file = self._get_file(write_model_part.GetCommunicator().GetDataCommunicator())
            HDF5PropertiesIO.Write(hdf5_file, "/Properties", write_model_part.Properties)

            # check if the last rank does not have any properties
            self.assertEqual(len(write_model_part.Properties) == 0, my_pid == num_proc - 1)

            read_model_part = current_model.CreateModelPart("read")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(read_model_part, Testing.GetDefaultDataCommunicator())
            HDF5PropertiesIO.Read(hdf5_file, "/Properties", read_model_part.Properties)

            write_prop: Properties
            for write_prop, read_prop in zip(write_model_part.Properties, read_model_part.Properties):
                for var in [DISTANCE, VELOCITY_X, VELOCITY_Y, VELOCITY_Z, DENSITY, ACCELERATION_X, ACCELERATION_Y, ACCELERATION_Z]:
                    if (write_prop.Has(var)): self.assertEqual(write_prop[var], read_prop[var])

                if (write_prop.Id <= number_of_properties):
                    self.assertTrue(write_prop.HasSubProperties(write_prop.Id + number_of_properties))
                    sub_write_prop = write_prop.GetSubProperties(write_prop.Id + number_of_properties)
                    sub_read_prop = read_prop.GetSubProperties(write_prop.Id + number_of_properties)
                    for var in [DISTANCE, VELOCITY_X, VELOCITY_Y, VELOCITY_Z, DENSITY, ACCELERATION_X, ACCELERATION_Y, ACCELERATION_Z]:
                        if (sub_write_prop.Has(var)): self.assertEqual(sub_write_prop[var], sub_read_prop[var])

    def test_HDF5ElementDataAndFlags(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            current_model = Model()
            write_model_part = current_model.CreateModelPart("write")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(write_model_part, Testing.GetDefaultDataCommunicator())
            self._initialize_model_part(write_model_part)

            data_communicator: DataCommunicator = write_model_part.GetCommunicator().GetDataCommunicator()

            my_pid = data_communicator.Rank()
            element: Element
            for element in write_model_part.Elements:
                element.SetValue(TEMPERATURE, random.random())
                element.Set(STRUCTURE, element.Id % 2)
                if my_pid % 2 == 0 and element.Id % 2 == 0:
                    element.SetValue(DENSITY, 1.0)
                    element.Set(SLIP, True)
                elif my_pid % 2 == 1 and element.Id % 2 == 1:
                    element.SetValue(DISPLACEMENT, Vector([random.random(), random.random(), random.random()]))
                    element.SetValue(DISTANCE, random.random())
                elif my_pid % 2 == 1 and element.Id % 2 == 0:
                    element.SetValue(RADIUS, 3.0)
                else:
                    element.SetValue(BOSSAK_ALPHA, 4.0)

            params = Parameters("""
            {
                "prefix" : "/ResultsData",
                "list_of_variables" : []
            }""")
            params["list_of_variables"].SetStringArray(Utilities.GetListOfAvailableVariables(write_model_part.Elements, write_model_part.GetCommunicator().GetDataCommunicator()))
            hdf5_file = self._get_file(write_model_part.GetCommunicator().GetDataCommunicator())
            data_value_io = HDF5ElementDataValueIO(params, hdf5_file)
            data_value_io.Write(write_model_part)

            params = Parameters("""
            {
                "prefix" : "/ResultsData",
                "list_of_variables" : ["SLIP", "STRUCTURE"]
            }""")
            flag_value_io = HDF5ElementFlagValueIO(params, hdf5_file)
            flag_value_io.Write(write_model_part)

            read_model_part = current_model.CreateModelPart("read")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(read_model_part, Testing.GetDefaultDataCommunicator())
            self._initialize_model_part(read_model_part)
            data_value_io.Read(read_model_part)
            flag_value_io.Read(read_model_part)

            assert_variables_list = [DENSITY, DISTANCE, RADIUS, BOSSAK_ALPHA, TEMPERATURE]
            assert_flags_list = [SLIP, STRUCTURE]
            for read_element, write_element in zip(read_model_part.Elements, write_model_part.Elements):
                for var in assert_variables_list:
                    if write_element.Has(var):
                        self.assertTrue(read_element.Has(var))
                        self.assertEqual(read_element.GetValue(var), write_element.GetValue(var))
                    else:
                        self.assertFalse(read_element.Has(var))

                for flag in assert_flags_list:
                    if write_element.IsDefined(flag):
                        self.assertTrue(read_element.IsDefined(flag))
                        self.assertEqual(read_element.Is(flag), write_element.Is(flag))
                    else:
                        self.assertFalse(read_element.IsDefined(flag))

            kratos_utilities.DeleteFileIfExisting("test_hdf5_model_part_io_mpi.h5")

    def test_HDF5ModelPartIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            current_model = Model()
            write_model_part = current_model.CreateModelPart("write")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(write_model_part, Testing.GetDefaultDataCommunicator())
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file(write_model_part.GetCommunicator().GetDataCommunicator())
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            read_model_part = current_model.CreateModelPart("read")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(read_model_part, Testing.GetDefaultDataCommunicator())
            hdf5_model_part_io.ReadModelPart(read_model_part)
            KratosMPI.ParallelFillCommunicator(read_model_part.GetRootModelPart(), Testing.GetDefaultDataCommunicator()).Execute()
            read_model_part.GetCommunicator().SynchronizeNodalSolutionStepsData()
            # Check nodes (node order should be preserved on read/write to ensure consistency with nodal results)
            self.assertEqual(read_model_part.NumberOfNodes(), write_model_part.NumberOfNodes())
            for read_node, write_node in zip(read_model_part.Nodes, write_model_part.Nodes):
                self.assertEqual(read_node.Id, write_node.Id)
                self.assertEqual(read_node.X, write_node.X)
                self.assertEqual(read_node.Y, write_node.Y)
                self.assertEqual(read_node.Z, write_node.Z)
            # Check elements
            self.assertEqual(read_model_part.NumberOfElements(), write_model_part.NumberOfElements())
            if read_model_part.NumberOfElements() > 0:
                first_elem_id = next(iter(read_model_part.Elements)).Id
                read_model_part.GetElement(first_elem_id) # Force a sort since order is mixed by openmp. TODO: to be removed once the PointerVectorSet is guaranteed to be sorted.
                for read_elem, write_elem in zip(read_model_part.Elements, write_model_part.Elements):
                    self.assertEqual(read_elem.Id, write_elem.Id)
                    self.assertEqual(read_elem.Properties.Id, write_elem.Properties.Id)
                    self.assertEqual(len(read_elem.GetNodes()), len(write_elem.GetNodes()))
                    for read_elem_node, write_elem_node in zip(read_elem.GetNodes(), write_elem.GetNodes()):
                        self.assertEqual(read_elem_node.Id, write_elem_node.Id)
            # Check conditions
            self.assertEqual(read_model_part.NumberOfConditions(), write_model_part.NumberOfConditions())
            if read_model_part.NumberOfConditions() > 0:
                first_cond_id = next(iter(read_model_part.Conditions)).Id
                read_model_part.GetCondition(first_cond_id) # Force a sort since order is mixed by openmp.  TODO: to be removed once the PointerVectorSet is guaranteed to be sorted.
                for read_cond, write_cond in zip(read_model_part.Conditions, write_model_part.Conditions):
                    self.assertEqual(read_cond.Id, write_cond.Id)
                    self.assertEqual(read_cond.Properties.Id, write_cond.Properties.Id)
                    self.assertEqual(len(read_cond.GetNodes()), len(write_cond.GetNodes()))
                    for read_cond_node, write_cond_node in zip(read_cond.GetNodes(), write_cond.GetNodes()):
                        self.assertEqual(read_cond_node.Id, write_cond_node.Id)
            # Check process info
            self.assertEqual(read_model_part.ProcessInfo[DOMAIN_SIZE], write_model_part.ProcessInfo[DOMAIN_SIZE])
            self.assertEqual(read_model_part.ProcessInfo[TIME], write_model_part.ProcessInfo[TIME])
            read_vector = read_model_part.ProcessInfo[INITIAL_STRAIN]
            write_vector = write_model_part.ProcessInfo[INITIAL_STRAIN]
            self.assertEqual(read_vector.Size(), write_vector.Size())
            for i in range(len(read_vector)):
                self.assertEqual(read_vector[i], write_vector[i])
            read_matrix = read_model_part.ProcessInfo[GREEN_LAGRANGE_STRAIN_TENSOR]
            write_matrix = write_model_part.ProcessInfo[GREEN_LAGRANGE_STRAIN_TENSOR]
            self.assertEqual(read_matrix.Size1(), write_matrix.Size1())
            self.assertEqual(read_matrix.Size2(), write_matrix.Size2())
            for i in range(read_matrix.Size1()):
                for j in range(read_matrix.Size2()):
                    self.assertEqual(read_matrix[i,j], write_matrix[i,j])
            kratos_utilities.DeleteFileIfExisting("test_hdf5_model_part_io_mpi.h5")

    def test_HDF5NodalSolutionStepDataIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            current_model = Model()
            write_model_part = current_model.CreateModelPart("write")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(write_model_part, Testing.GetDefaultDataCommunicator())
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file(write_model_part.GetCommunicator().GetDataCommunicator())
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            read_model_part = current_model.CreateModelPart("read")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(read_model_part, Testing.GetDefaultDataCommunicator())
            hdf5_model_part_io.ReadModelPart(read_model_part)
            KratosMPI.ParallelFillCommunicator(read_model_part.GetRootModelPart(), Testing.GetDefaultDataCommunicator()).Execute()
            hdf5_nodal_solution_step_data_io = self._get_nodal_solution_step_data_io(hdf5_file)
            hdf5_nodal_solution_step_data_io.Write(write_model_part, Parameters("""{}"""), 0)
            hdf5_nodal_solution_step_data_io.Read(read_model_part, 0)

            # Check data.
            for read_node, write_node in zip(read_model_part.Nodes, write_model_part.Nodes):
                self.assertEqual(read_node.GetSolutionStepValue(DISPLACEMENT_X), write_node.GetSolutionStepValue(DISPLACEMENT_X))
                self.assertEqual(read_node.GetSolutionStepValue(DISPLACEMENT_Y), write_node.GetSolutionStepValue(DISPLACEMENT_Y))
                self.assertEqual(read_node.GetSolutionStepValue(DISPLACEMENT_Z), write_node.GetSolutionStepValue(DISPLACEMENT_Z))
                self.assertEqual(read_node.GetSolutionStepValue(VELOCITY_X), write_node.GetSolutionStepValue(VELOCITY_X))
                self.assertEqual(read_node.GetSolutionStepValue(VELOCITY_Y), write_node.GetSolutionStepValue(VELOCITY_Y))
                self.assertEqual(read_node.GetSolutionStepValue(VELOCITY_Z), write_node.GetSolutionStepValue(VELOCITY_Z))
                self.assertEqual(read_node.GetSolutionStepValue(ACCELERATION_X), write_node.GetSolutionStepValue(ACCELERATION_X))
                self.assertEqual(read_node.GetSolutionStepValue(ACCELERATION_Y), write_node.GetSolutionStepValue(ACCELERATION_Y))
                self.assertEqual(read_node.GetSolutionStepValue(ACCELERATION_Z), write_node.GetSolutionStepValue(ACCELERATION_Z))
                self.assertEqual(read_node.GetSolutionStepValue(PRESSURE), write_node.GetSolutionStepValue(PRESSURE))
                self.assertEqual(read_node.GetSolutionStepValue(VISCOSITY), write_node.GetSolutionStepValue(VISCOSITY))
                self.assertEqual(read_node.GetSolutionStepValue(DENSITY), write_node.GetSolutionStepValue(DENSITY))
                self.assertEqual(read_node.GetSolutionStepValue(ACTIVATION_LEVEL), write_node.GetSolutionStepValue(ACTIVATION_LEVEL))
                self.assertEqual(read_node.GetSolutionStepValue(PARTITION_INDEX), write_node.GetSolutionStepValue(PARTITION_INDEX))
            kratos_utilities.DeleteFileIfExisting("test_hdf5_model_part_io_mpi.h5")

    def test_HDF5NodalDataIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            current_model = Model()
            write_model_part = current_model.CreateModelPart("write")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(write_model_part, Testing.GetDefaultDataCommunicator())
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file(write_model_part.GetCommunicator().GetDataCommunicator())
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            read_model_part = current_model.CreateModelPart("read")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(read_model_part, Testing.GetDefaultDataCommunicator())
            hdf5_model_part_io.ReadModelPart(read_model_part)
            KratosMPI.ParallelFillCommunicator(read_model_part.GetRootModelPart(), Testing.GetDefaultDataCommunicator()).Execute()
            hdf5_nodal_data_io = self._get_nodal_data_io(hdf5_file)
            hdf5_nodal_data_io.Write(write_model_part)
            hdf5_nodal_data_io.Read(read_model_part)

            # # Check data.
            for read_node, write_node in zip(read_model_part.Nodes, write_model_part.Nodes):
                self.assertEqual(read_node.GetValue(DISPLACEMENT_X), write_node.GetValue(DISPLACEMENT_X))
                self.assertEqual(read_node.GetValue(DISPLACEMENT_Y), write_node.GetValue(DISPLACEMENT_Y))
                self.assertEqual(read_node.GetValue(DISPLACEMENT_Z), write_node.GetValue(DISPLACEMENT_Z))
                self.assertEqual(read_node.GetValue(VELOCITY_X), write_node.GetValue(VELOCITY_X))
                self.assertEqual(read_node.GetValue(VELOCITY_Y), write_node.GetValue(VELOCITY_Y))
                self.assertEqual(read_node.GetValue(VELOCITY_Z), write_node.GetValue(VELOCITY_Z))
                self.assertEqual(read_node.GetValue(ACCELERATION_X), write_node.GetValue(ACCELERATION_X))
                self.assertEqual(read_node.GetValue(ACCELERATION_Y), write_node.GetValue(ACCELERATION_Y))
                self.assertEqual(read_node.GetValue(ACCELERATION_Z), write_node.GetValue(ACCELERATION_Z))
                self.assertEqual(read_node.GetValue(PRESSURE), write_node.GetValue(PRESSURE))
                self.assertEqual(read_node.GetValue(VISCOSITY), write_node.GetValue(VISCOSITY))
                self.assertEqual(read_node.GetValue(DENSITY), write_node.GetValue(DENSITY))
                self.assertEqual(read_node.GetValue(ACTIVATION_LEVEL), write_node.GetValue(ACTIVATION_LEVEL))
            kratos_utilities.DeleteFileIfExisting("test_hdf5_model_part_io_mpi.h5")

    def test_HDF5NodalFlagIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            current_model = Model()
            write_model_part = current_model.CreateModelPart("write")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(write_model_part, Testing.GetDefaultDataCommunicator())
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file(write_model_part.GetCommunicator().GetDataCommunicator())
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            read_model_part = current_model.CreateModelPart("read")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(read_model_part, Testing.GetDefaultDataCommunicator())
            hdf5_model_part_io.ReadModelPart(read_model_part)
            KratosMPI.ParallelFillCommunicator(read_model_part.GetRootModelPart(), Testing.GetDefaultDataCommunicator()).Execute()
            hdf5_nodal_flag_io = self._get_nodal_flag_io(hdf5_file)
            hdf5_nodal_flag_io.Write(write_model_part)
            hdf5_nodal_flag_io.Read(read_model_part)

            # Check flag.
            for read_node, write_node in zip(read_model_part.Nodes, write_model_part.Nodes):
                self.assertEqual(read_node.Is(SLIP), write_node.Is(SLIP))
                self.assertEqual(read_node.Is(ACTIVE), write_node.Is(ACTIVE))
            kratos_utilities.DeleteFileIfExisting("test_hdf5_model_part_io_mpi.h5")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
