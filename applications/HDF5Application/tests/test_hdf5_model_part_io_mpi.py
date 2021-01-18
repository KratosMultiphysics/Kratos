import os
from KratosMultiphysics import *
import KratosMultiphysics.mpi as KratosMPI
from KratosMultiphysics.HDF5Application import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import random, math

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

    def _initialize_model_part(self, model_part):
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
        num_proc = DataCommunicator.GetDefault().Size()
        my_pid = DataCommunicator.GetDefault().Rank()
        my_num_quad = 20 # user-defined.
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

        KratosMPI.ParallelFillCommunicator(model_part.GetRootModelPart()).Execute()
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

    def _get_file(self):
        params = Parameters("""
        {
            "file_name" : "test_hdf5_model_part_io_mpi.h5",
            "file_access_mode" : "truncate",
            "file_driver" : "mpio",
            "echo_level" : 0
        }""")
        return HDF5FileParallel(params)

    def _get_model_part_io(self, hdf5_file):
        return HDF5PartitionedModelPartIO(hdf5_file, "/ModelData")

    def _get_nodal_solution_step_data_io(self, hdf5_file):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["DISPLACEMENT", "VELOCITY", "ACCELERATION", "PRESSURE", "VISCOSITY", "DENSITY", "ACTIVATION_LEVEL"]
        }""")
        return HDF5NodalSolutionStepDataIO(params, hdf5_file)

    def _get_nodal_data_io(self, hdf5_file):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["DISPLACEMENT", "VELOCITY", "ACCELERATION", "PRESSURE", "VISCOSITY", "DENSITY", "ACTIVATION_LEVEL"]
        }""")
        return HDF5NodalDataValueIO(params, hdf5_file)

    def _get_nodal_flag_io(self, hdf5_file):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["SLIP", "ACTIVE"]
        }""")
        return HDF5NodalFlagValueIO(params, hdf5_file)

    def test_HDF5ModelPartIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            current_model = Model()
            write_model_part = current_model.CreateModelPart("write")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(write_model_part)
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            read_model_part = current_model.CreateModelPart("read")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(read_model_part)
            hdf5_model_part_io.ReadModelPart(read_model_part)
            KratosMPI.ParallelFillCommunicator(read_model_part.GetRootModelPart()).Execute()
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
            first_elem_id = next(iter(read_model_part.Elements)).Id
            read_model_part.GetElement(first_elem_id) # Force a sort since order is mixed by openmp.
            for read_elem, write_elem in zip(read_model_part.Elements, write_model_part.Elements):
                self.assertEqual(read_elem.Id, write_elem.Id)
                self.assertEqual(read_elem.Properties.Id, write_elem.Properties.Id)
                self.assertEqual(len(read_elem.GetNodes()), len(write_elem.GetNodes()))
                for read_elem_node, write_elem_node in zip(read_elem.GetNodes(), write_elem.GetNodes()):
                    self.assertEqual(read_elem_node.Id, write_elem_node.Id)
            # Check conditions
            self.assertEqual(read_model_part.NumberOfConditions(), write_model_part.NumberOfConditions())
            first_cond_id = next(iter(read_model_part.Conditions)).Id
            read_model_part.GetCondition(first_cond_id) # Force a sort since order is mixed by openmp.
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
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(write_model_part)
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            read_model_part = current_model.CreateModelPart("read")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(read_model_part)
            hdf5_model_part_io.ReadModelPart(read_model_part)
            KratosMPI.ParallelFillCommunicator(read_model_part.GetRootModelPart()).Execute()
            hdf5_nodal_solution_step_data_io = self._get_nodal_solution_step_data_io(hdf5_file)
            hdf5_nodal_solution_step_data_io.WriteNodalResults(write_model_part.Nodes, 0)
            hdf5_nodal_solution_step_data_io.ReadNodalResults(read_model_part.Nodes, read_model_part.GetCommunicator(), 0)

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
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(write_model_part)
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            read_model_part = current_model.CreateModelPart("read")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(read_model_part)
            hdf5_model_part_io.ReadModelPart(read_model_part)
            KratosMPI.ParallelFillCommunicator(read_model_part.GetRootModelPart()).Execute()
            hdf5_nodal_data_io = self._get_nodal_data_io(hdf5_file)
            hdf5_nodal_data_io.WriteNodalResults(write_model_part.Nodes)
            hdf5_nodal_data_io.ReadNodalResults(read_model_part.Nodes, read_model_part.GetCommunicator())

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
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(write_model_part)
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            read_model_part = current_model.CreateModelPart("read")
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(read_model_part)
            hdf5_model_part_io.ReadModelPart(read_model_part)
            KratosMPI.ParallelFillCommunicator(read_model_part.GetRootModelPart()).Execute()
            hdf5_nodal_flag_io = self._get_nodal_flag_io(hdf5_file)
            hdf5_nodal_flag_io.WriteNodalFlags(write_model_part.Nodes)
            hdf5_nodal_flag_io.ReadNodalFlags(read_model_part.Nodes, read_model_part.GetCommunicator())

            # # Check flag.
            for read_node, write_node in zip(read_model_part.Nodes, write_model_part.Nodes):
                self.assertEqual(read_node.Is(SLIP), write_node.Is(SLIP))
                self.assertEqual(read_node.Is(ACTIVE), write_node.Is(ACTIVE))
            kratos_utilities.DeleteFileIfExisting("test_hdf5_model_part_io_mpi.h5")

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
