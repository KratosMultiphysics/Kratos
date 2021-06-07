import os
from KratosMultiphysics import *
from KratosMultiphysics.HDF5Application import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import random

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
        num_tri_elems = 10 # Use enough elements for detecting possible openmp issues.
        num_quad_elems = 15
        num_tri_conds = 5
        num_quad_conds = 10
        num_nodes = (num_tri_elems + 2) + (num_quad_elems + 2)
        # Add variables.
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT) # array_1d
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(ACCELERATION)
        model_part.AddNodalSolutionStepVariable(PRESSURE) # double
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(ACTIVATION_LEVEL) # int
        # Create nodes.
        for i in range(num_nodes):
            x = random.random()
            y = random.random()
            z = random.random()
            model_part.CreateNewNode(i + 1, x, y, z)
        # Create triangle elements.
        for i in range(num_tri_elems):
            elem_id = i + 1
            node_ids = [i + 1, i + 2, i + 3]
            prop_id = random.randint(0,4)
            prop = model_part.GetProperties()[prop_id]
            model_part.CreateNewElement("Element2D3N", elem_id, node_ids, prop)
        # Create quad elements.
        for i in range(num_quad_elems):
            elem_id = num_tri_elems + i + 1
            node_ids = [num_tri_elems + i + 1, num_tri_elems + i + 2, num_tri_elems + i + 3, num_tri_elems + i + 4]
            prop_id = random.randint(0,4)
            prop = model_part.GetProperties()[prop_id]
            model_part.CreateNewElement("Element2D4N", elem_id, node_ids, prop)
        # Create triangle conditions.
        for i in range(num_tri_conds):
            cond_id = i + 1
            node_ids = [i + 1, i + 2, i + 3]
            prop_id = random.randint(0,4)
            prop = model_part.GetProperties()[prop_id]
            model_part.CreateNewCondition("SurfaceCondition3D3N", cond_id, node_ids, prop)
        # Create quad conditions.
        for i in range(num_quad_conds):
            cond_id = num_tri_conds + i + 1
            node_ids = [num_tri_conds + i + 1, num_tri_conds + i + 2, num_tri_conds + i + 3, num_tri_conds + i + 4]
            prop_id = random.randint(0,4)
            prop = model_part.GetProperties()[prop_id]
            model_part.CreateNewCondition("SurfaceCondition3D4N", cond_id, node_ids, prop)
        model_part.SetBufferSize(2)
        for node in model_part.Nodes:
            # Write some data to the nodal solution steps variables.
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
            # Write some data to the nodal data container variables.
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
            # Write some flags
            node.Set(SLIP, bool(random.randint(-100, 100) % 2))
            node.Set(ACTIVE, bool(random.randint(-100, 100) % 2))
            node.Set(STRUCTURE, bool(random.randint(-100, 100) % 2))

        for element in model_part.Elements:
            element.SetValue(ACCELERATION, Vector([random.random(), random.random(), random.random()]))
            element.SetValue(PRESSURE, random.random())
            element.SetValue(VISCOSITY, random.random())
            element.SetValue(DENSITY, random.random())
            element.SetValue(ACTIVATION_LEVEL, random.randint(-100, 100))

            element.Set(SLIP, bool(random.randint(-100, 100) % 2))
            element.Set(ACTIVE, bool(random.randint(-100, 100) % 2))
            element.Set(STRUCTURE, bool(random.randint(-100, 100) % 2))

        for condition in model_part.Conditions:
            condition.SetValue(ACCELERATION, Vector([random.random(), random.random(), random.random()]))
            condition.SetValue(PRESSURE, random.random())
            condition.SetValue(VISCOSITY, random.random())
            condition.SetValue(DENSITY, random.random())
            condition.SetValue(ACTIVATION_LEVEL, random.randint(-100, 100))

            condition.Set(SLIP, bool(random.randint(-100, 100) % 2))
            condition.Set(ACTIVE, bool(random.randint(-100, 100) % 2))
            condition.Set(STRUCTURE, bool(random.randint(-100, 100) % 2))
        # Set some process info variables.
        model_part.ProcessInfo[DOMAIN_SIZE] = 3 # int
        model_part.ProcessInfo[TIME] = random.random() # float
        initial_strain = Vector(6)
        for i in range(6):
            initial_strain[i] = random.random()
        model_part.ProcessInfo[INITIAL_STRAIN] = initial_strain # vector
        gl_strain_tensor = Matrix(5,5)
        for i in range(5):
            for j in range(5):
                gl_strain_tensor[i,j] = random.random()
        model_part.ProcessInfo[GREEN_LAGRANGE_STRAIN_TENSOR] = gl_strain_tensor # matrix

    def _get_file(self):
        params = Parameters("""
        {
            "file_name" : "test_hdf5_model_part_io.h5",
            "file_access_mode" : "exclusive",
            "file_driver" : "core"
        }""")
        return HDF5FileSerial(params)

    def _get_model_part_io(self, hdf5_file):
        return HDF5ModelPartIO(hdf5_file, "/ModelData")

    def _get_nodal_solution_step_data_io(self, hdf5_file):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["DISPLACEMENT", "VELOCITY", "ACCELERATION", "PRESSURE", "VISCOSITY", "DENSITY", "ACTIVATION_LEVEL"]
        }""")
        return HDF5NodalSolutionStepDataIO(params, hdf5_file)

    def _get_element_data_value_io(self, hdf5_file):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["DISPLACEMENT", "VELOCITY", "ACCELERATION", "PRESSURE", "VISCOSITY", "DENSITY", "ACTIVATION_LEVEL"]
        }""")
        return HDF5ElementDataValueIO(params, hdf5_file)

    def _get_element_flag_value_io(self, hdf5_file):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["SLIP", "ACTIVE", "STRUCTURE"]
        }""")
        return HDF5ElementFlagValueIO(params, hdf5_file)

    def _get_condition_data_value_io(self, hdf5_file):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["DISPLACEMENT", "VELOCITY", "ACCELERATION", "PRESSURE", "VISCOSITY", "DENSITY", "ACTIVATION_LEVEL"]
        }""")
        return HDF5ConditionDataValueIO(params, hdf5_file)

    def _get_condition_flag_value_io(self, hdf5_file):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["SLIP", "ACTIVE", "STRUCTURE"]
        }""")
        return HDF5ConditionFlagValueIO(params, hdf5_file)

    def _get_nodal_data_value_io(self, hdf5_file):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["DISPLACEMENT", "VELOCITY", "ACCELERATION", "PRESSURE", "VISCOSITY", "DENSITY", "ACTIVATION_LEVEL"]
        }""")
        return HDF5NodalDataValueIO(params, hdf5_file)

    def _get_nodal_flag_value_io(self, hdf5_file):
        params = Parameters("""
        {
            "prefix" : "/ResultsData",
            "list_of_variables" : ["SLIP", "ACTIVE", "STRUCTURE"]
        }""")
        return HDF5NodalFlagValueIO(params, hdf5_file)

    def test_HDF5ModelPartIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model = Model()
            write_model_part = model.CreateModelPart("write", 2)
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            read_model_part = model.CreateModelPart("read", 2)
            hdf5_model_part_io.ReadModelPart(read_model_part)
            # Check nodes (node order should be preserved on read/write to ensure consistency with nodal results)
            self.assertEqual(read_model_part.NumberOfNodes(), write_model_part.NumberOfNodes())
            for read_node, write_node in zip(read_model_part.Nodes, write_model_part.Nodes):
                self.assertEqual(read_node.Id, write_node.Id)
                self.assertEqual(read_node.X, write_node.X)
                self.assertEqual(read_node.Y, write_node.Y)
                self.assertEqual(read_node.Z, write_node.Z)
            # Check elements
            self.assertEqual(read_model_part.NumberOfElements(), write_model_part.NumberOfElements())
            read_model_part.GetElement(1) # Force a sort since order is mixed by openmp.
            for read_elem, write_elem in zip(read_model_part.Elements, write_model_part.Elements):
                self.assertEqual(read_elem.Id, write_elem.Id)
                self.assertEqual(read_elem.Properties.Id, write_elem.Properties.Id)
                self.assertEqual(len(read_elem.GetNodes()), len(write_elem.GetNodes()))
                for read_elem_node, write_elem_node in zip(read_elem.GetNodes(), write_elem.GetNodes()):
                    self.assertEqual(read_elem_node.Id, write_elem_node.Id)
            # Check conditions
            self.assertEqual(read_model_part.NumberOfConditions(), write_model_part.NumberOfConditions())
            read_model_part.GetCondition(1) # Force a sort since order is mixed by openmp.
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

    def test_RecursiveSubModelParts(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model = Model()

            write_model_part = model.CreateModelPart("write_model_part")
            def create_entities(model_part, num_nodes):
                global_number_of_nodes = write_model_part.GetCommunicator().GlobalNumberOfNodes()
                for i in range(num_nodes):
                    x = random.random()
                    y = random.random()
                    z = random.random()
                    model_part.CreateNewNode(global_number_of_nodes + i + 1, x, y, z)

                prop = prop = model_part.GetProperties()[0]
                global_number_of_elements = write_model_part.GetCommunicator().GlobalNumberOfElements()
                for i in range(num_nodes - 2):
                    elem_id = global_number_of_elements + i + 1
                    node_ids = [global_number_of_nodes + i + 1, global_number_of_nodes + i + 2, global_number_of_nodes + i + 3]
                    model_part.CreateNewElement("Element3D3N", elem_id, node_ids, prop)

                global_number_of_elements = write_model_part.GetCommunicator().GlobalNumberOfElements()
                for i in range(num_nodes - 3):
                    elem_id = global_number_of_elements + i + 1
                    node_ids = [global_number_of_nodes + i + 1, global_number_of_nodes + i + 2, global_number_of_nodes + i + 3, global_number_of_nodes + i + 4]
                    model_part.CreateNewElement("Element3D4N", elem_id, node_ids, prop)

                global_number_of_conditions = write_model_part.GetCommunicator().GlobalNumberOfConditions()
                for i in range(num_nodes - 2):
                    cond_id = global_number_of_conditions + i + 1
                    node_ids = [global_number_of_nodes + i + 1, global_number_of_nodes + i + 2]
                    model_part.CreateNewCondition("LineCondition2D2N", cond_id, node_ids, prop)

                global_number_of_conditions = write_model_part.GetCommunicator().GlobalNumberOfConditions()
                for i in range(num_nodes - 3):
                    cond_id = global_number_of_conditions + i + 1
                    node_ids = [global_number_of_nodes + i + 1, global_number_of_nodes + i + 2, global_number_of_nodes + i + 3]
                    model_part.CreateNewCondition("SurfaceCondition3D3N", cond_id, node_ids, prop)

            create_entities(write_model_part, 3)
            sub_model_part = write_model_part.CreateSubModelPart("sub_model_part1")
            create_entities(sub_model_part, 5)
            sub_sub_model_part = sub_model_part.CreateSubModelPart("section_1")
            create_entities(sub_sub_model_part, 3)
            sub_sub_model_part = sub_sub_model_part.CreateSubModelPart("section_2")
            create_entities(sub_sub_model_part, 4)
            sub_sub_model_part = sub_sub_model_part.CreateSubModelPart("section_3")
            create_entities(sub_sub_model_part, 5)
            sub_sub_model_part = sub_model_part.CreateSubModelPart("section_2")
            create_entities(sub_sub_model_part, 4)

            sub_model_part = write_model_part.CreateSubModelPart("sub_model_part2")
            create_entities(sub_model_part, 3)
            sub_sub_model_part = sub_model_part.CreateSubModelPart("section_1")
            create_entities(sub_sub_model_part, 3)
            sub_sub_model_part = sub_sub_model_part.CreateSubModelPart("section_2")
            create_entities(sub_sub_model_part, 6)
            sub_sub_model_part = sub_sub_model_part.CreateSubModelPart("section_3")
            create_entities(sub_sub_model_part, 3)
            sub_sub_model_part = sub_model_part.CreateSubModelPart("section_2")
            create_entities(sub_sub_model_part, 3)

            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)

            read_model_part = model.CreateModelPart("read_model_part")
            hdf5_model_part_io.ReadModelPart(read_model_part)

            def check_model_parts(model_part1, model_part2):
                def check_nodes(nodes_1, nodes_2):
                    self.assertEqual(len(nodes_1), len(nodes_2))
                    for n1, n2 in zip(nodes_1, nodes_2):
                        self.assertEqual(n1.Id, n2.Id)
                        self.assertEqual(n1.X, n2.X)
                        self.assertEqual(n1.Y, n2.Y)
                        self.assertEqual(n1.Z, n2.Z)

                check_nodes(model_part1.Nodes, model_part2.Nodes)

                self.assertEqual(model_part1.NumberOfElements(), model_part2.NumberOfElements())
                for e1, e2 in zip(model_part1.Elements, model_part2.Elements):
                    self.assertEqual(e1.Id, e2.Id)
                    check_nodes(e1.GetNodes(), e2.GetNodes())

                self.assertEqual(model_part1.NumberOfConditions(), model_part2.NumberOfConditions())
                for c1, c2 in zip(model_part1.Conditions, model_part2.Conditions):
                    self.assertEqual(c1.Id, c2.Id)
                    check_nodes(c1.GetNodes(), c2.GetNodes())

            check_model_parts(read_model_part, write_model_part)
            def recursive_submodel_parts(sub_model_parts_1, sub_model_parts_2):
                for sub_model_part_1, sub_model_part_2 in zip(sub_model_parts_1, sub_model_parts_2):
                    self.assertEqual(sub_model_part_1.Name, sub_model_part_2.Name)
                    check_model_parts(sub_model_part_1, sub_model_part_2)
                    self.assertEqual(sub_model_part_1.NumberOfSubModelParts(), sub_model_part_2.NumberOfSubModelParts())
                    recursive_submodel_parts(sub_model_part_1.SubModelParts, sub_model_part_2.SubModelParts)

    def test_HDF5NodalSolutionStepDataIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model = Model()
            write_model_part = model.CreateModelPart("write", 2)
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_nodal_solution_step_data_io = self._get_nodal_solution_step_data_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            hdf5_nodal_solution_step_data_io.WriteNodalResults(write_model_part, 0)
            read_model_part = model.CreateModelPart("read", 2)
            hdf5_model_part_io.ReadModelPart(read_model_part)
            hdf5_nodal_solution_step_data_io.ReadNodalResults(read_model_part, 0)

            assert_variables_list = [DISPLACEMENT_X, DISPLACEMENT_Y, DISPLACEMENT_Z,
                                     VELOCITY_X, VELOCITY_Y, VELOCITY_Z,
                                     ACCELERATION_X, ACCELERATION_Y, ACCELERATION_Z,
                                     PRESSURE, VISCOSITY, DENSITY, ACTIVATION_LEVEL]
            # Check data.
            for read_node, write_node in zip(read_model_part.Nodes, write_model_part.Nodes):
                for var in assert_variables_list:
                    self.assertEqual(read_node.GetSolutionStepValue(var), write_node.GetSolutionStepValue(var))

    def test_HDF5ElementDataValueIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model = Model()
            write_model_part = model.CreateModelPart("write", 2)
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_element_data_value_io = self._get_element_data_value_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            hdf5_element_data_value_io.WriteElementResults(write_model_part.Elements)
            read_model_part = model.CreateModelPart("read", 2)
            hdf5_model_part_io.ReadModelPart(read_model_part)
            hdf5_element_data_value_io.ReadElementResults(read_model_part.Elements,
                                                          read_model_part.GetCommunicator())

            assert_variables_list = [PRESSURE, VISCOSITY, DENSITY, ACTIVATION_LEVEL]
            # Check data.
            for read_element, write_element in zip(read_model_part.Elements, write_model_part.Elements):
                for var in assert_variables_list:
                    self.assertEqual(read_element.GetValue(var), write_element.GetValue(var))

    def test_HDF5ElementFlagValueIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model = Model()
            write_model_part = model.CreateModelPart("write", 2)
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_element_flag_value_io = self._get_element_flag_value_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            hdf5_element_flag_value_io.WriteElementFlags(write_model_part.Elements)
            read_model_part = model.CreateModelPart("read", 2)
            hdf5_model_part_io.ReadModelPart(read_model_part)
            hdf5_element_flag_value_io.ReadElementFlags(read_model_part.Elements,
                                                        read_model_part.GetCommunicator())

            assert_variables_list = [SLIP, ACTIVE, STRUCTURE]
            # Check data.
            for read_element, write_element in zip(read_model_part.Elements, write_model_part.Elements):
                for var in assert_variables_list:
                    self.assertEqual(read_element.Is(var), write_element.Is(var))

    def test_HDF5ConditionDataValueIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model = Model()
            write_model_part = model.CreateModelPart("write", 2)
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_condition_data_value_io = self._get_condition_data_value_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            hdf5_condition_data_value_io.WriteConditionResults(write_model_part.Conditions)
            read_model_part = model.CreateModelPart("read", 2)
            hdf5_model_part_io.ReadModelPart(read_model_part)
            hdf5_condition_data_value_io.ReadConditionResults(read_model_part.Conditions,
                                                              read_model_part.GetCommunicator())

            assert_variables_list = [PRESSURE, VISCOSITY, DENSITY, ACTIVATION_LEVEL]
            # Check data.
            for read_condition, write_condition in zip(read_model_part.Conditions, write_model_part.Conditions):
                for var in assert_variables_list:
                    self.assertEqual(read_condition.GetValue(var), write_condition.GetValue(var))

    def test_HDF5ConditionFlagValueIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model = Model()
            write_model_part = model.CreateModelPart("write", 2)
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_condition_flag_value_io = self._get_condition_flag_value_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            hdf5_condition_flag_value_io.WriteConditionFlags(write_model_part.Conditions)
            read_model_part = model.CreateModelPart("read", 2)
            hdf5_model_part_io.ReadModelPart(read_model_part)
            hdf5_condition_flag_value_io.ReadConditionFlags(read_model_part.Conditions,
                                                            read_model_part.GetCommunicator())

            assert_variables_list = [SLIP, ACTIVE, STRUCTURE]
            # Check data.
            for read_condition, write_condition in zip(read_model_part.Conditions, write_model_part.Conditions):
                for var in assert_variables_list:
                    self.assertEqual(read_condition.Is(var), write_condition.Is(var))

    def test_HDF5NodalDataValueIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model = Model()
            write_model_part = model.CreateModelPart("write", 2)
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_nodal_data_value_io = self._get_nodal_data_value_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            hdf5_nodal_data_value_io.WriteNodalResults(write_model_part.Nodes)
            read_model_part = model.CreateModelPart("read", 2)
            hdf5_model_part_io.ReadModelPart(read_model_part)
            hdf5_nodal_data_value_io.ReadNodalResults(read_model_part.Nodes, read_model_part.GetCommunicator())

            assert_variables_list = [DISPLACEMENT_X, DISPLACEMENT_Y, DISPLACEMENT_Z,
                                     VELOCITY_X, VELOCITY_Y, VELOCITY_Z,
                                     ACCELERATION_X, ACCELERATION_Y, ACCELERATION_Z,
                                     PRESSURE, VISCOSITY, DENSITY, ACTIVATION_LEVEL]
            # Check data.
            for read_node, write_node in zip(read_model_part.Nodes, write_model_part.Nodes):
                for var in assert_variables_list:
                    self.assertEqual(read_node.GetValue(var), write_node.GetValue(var))

    def test_HDF5NodalFlagValueIO(self):
        with ControlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            model = Model()
            write_model_part = model.CreateModelPart("write", 2)
            self._initialize_model_part(write_model_part)
            hdf5_file = self._get_file()
            hdf5_model_part_io = self._get_model_part_io(hdf5_file)
            hdf5_nodal_flag_value_io = self._get_nodal_flag_value_io(hdf5_file)
            hdf5_model_part_io.WriteModelPart(write_model_part)
            hdf5_nodal_flag_value_io.WriteNodalFlags(write_model_part.Nodes)
            read_model_part = model.CreateModelPart("read", 2)
            hdf5_model_part_io.ReadModelPart(read_model_part)
            hdf5_nodal_flag_value_io.ReadNodalFlags(read_model_part.Nodes, read_model_part.GetCommunicator())

            assert_variables_list = [SLIP, ACTIVE, STRUCTURE]
            # Check flag.
            for read_node, write_node in zip(read_model_part.Nodes, write_model_part.Nodes):
                for var in assert_variables_list:
                    self.assertEqual(read_node.Is(var), write_node.Is(var))

    def tearDown(self):
        pass

if __name__ == '__main__':
    KratosUnittest.main()
