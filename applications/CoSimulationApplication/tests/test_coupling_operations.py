import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
from KratosMultiphysics.CoSimulationApplication.factories import coupling_operation_factory
from testing_utilities import DummySolverWrapper

from math import sqrt, pi

class TestScalingOperation(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("default")
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part.ProcessInfo[KM.TIME] = 0.0
        self.model_part.ProcessInfo[KM.STEP] = 0

        for i in range(5):
            new_node = self.model_part.CreateNewNode(i+1, i*0.1, 0.0, 0.0)
            new_node.SetSolutionStepValue(KM.PRESSURE, 0, i+1.3)

        data_settings = KM.Parameters("""{
            "model_part_name" : "default",
            "variable_name"   : "PRESSURE"
        }""")
        self.interface_data = CouplingInterfaceData(data_settings, self.model)

        self.solver_wrappers = {"dummy_solver" : DummySolverWrapper({"data_4_testing" : self.interface_data})}

        self.solver_process_info = KM.ProcessInfo()

    def test_constant_scaling(self):
        scaling_op_settings = KM.Parameters("""{
            "type"           : "scaling",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "scaling_factor" : 1.5,
            "echo_level"     : 0
        }""")

        scaling_op = coupling_operation_factory.CreateCouplingOperation(scaling_op_settings, self.solver_wrappers, self.solver_process_info, KM.Testing.GetDefaultDataCommunicator())

        factors = [1.5] * 3

        self.__ExecuteTest(scaling_op, factors)

    def test_constant_scaling_from_string(self):
        scaling_op_settings = KM.Parameters("""{
            "type"           : "scaling",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "scaling_factor" : "1.5",
            "echo_level"     : 0
        }""")

        scaling_op = coupling_operation_factory.CreateCouplingOperation(scaling_op_settings, self.solver_wrappers, self.solver_process_info, KM.Testing.GetDefaultDataCommunicator())

        factors = [1.5] * 3

        self.__ExecuteTest(scaling_op, factors)

    def test_variable_scaling_time(self):
        scaling_op_settings = KM.Parameters("""{
            "type"           : "scaling",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "scaling_factor" : "1.5*t",
            "echo_level"     : 0
        }""")

        scaling_op = coupling_operation_factory.CreateCouplingOperation(scaling_op_settings, self.solver_wrappers, self.solver_process_info, KM.Testing.GetDefaultDataCommunicator())

        factors = [1.5*0.25, 1.5*0.5, 1.5*0.75]

        self.__ExecuteTest(scaling_op, factors)

    def test_variable_scaling_step(self):
        scaling_op_settings = KM.Parameters("""{
            "type"           : "scaling",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "scaling_factor" : "1.5*sqrt(step)*pi",
            "echo_level"     : 0
        }""")

        scaling_op = coupling_operation_factory.CreateCouplingOperation(scaling_op_settings, self.solver_wrappers, self.solver_process_info, KM.Testing.GetDefaultDataCommunicator())

        factors = [1.5*pi*sqrt(1), 1.5*pi*sqrt(2), 1.5*pi*sqrt(3), 1.5*pi*sqrt(4), 1.5*pi*sqrt(5)]

        self.__ExecuteTest(scaling_op, factors)

    def test_scaling_in_interval(self):
        scaling_op_settings = KM.Parameters("""{
            "type"           : "scaling",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "scaling_factor" : 1.22,
            "interval"       : [0.0, 0.3]
        }""")

        scaling_op = coupling_operation_factory.CreateCouplingOperation(scaling_op_settings, self.solver_wrappers, self.solver_process_info, KM.Testing.GetDefaultDataCommunicator())

        factors = [1.0] * 5
        factors[0] = 1.22

        self.__ExecuteTest(scaling_op, factors)

    def test_scaling_in_interval_2(self):
        scaling_op_settings = KM.Parameters("""{
            "type"           : "scaling",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "scaling_factor" : 1.22,
            "interval"       : [0.8, "End"]
        }""")

        scaling_op = coupling_operation_factory.CreateCouplingOperation(scaling_op_settings, self.solver_wrappers, self.solver_process_info, KM.Testing.GetDefaultDataCommunicator())

        factors = [1.0] * 3
        factors.extend([1.22] * 3)

        self.__ExecuteTest(scaling_op, factors)

    def __ExecuteTest(self, scaling_operation, factors):
        scaling_operation.Check()

        for fac in factors:
            self.model_part.ProcessInfo[KM.TIME] += 0.25
            self.model_part.ProcessInfo[KM.STEP] += 1

            old_data = self.interface_data.GetData()

            scaling_operation.Execute()

            new_data = self.interface_data.GetData()

            self.__CompareValues(old_data, new_data, fac)

    def __CompareValues(self, old_data, new_data, factor):
        for old_val, new_val in zip(old_data, new_data):
            self.assertAlmostEqual(old_val*factor, new_val)

class TestConversionOperation(KratosUnittest.TestCase):

    def test_elemental_to_nodal_conversion(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("default")
        self.model_part.AddNodalSolutionStepVariable(KM.FORCE)
        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
        props = self.model_part.CreateNewProperties(1)

        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        self.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        new_element = self.model_part.CreateNewElement("Element2D4N", 1, [1,2,3,4], props)
        new_element.SetValue(KM.FORCE, [12.0, 8.0, 0.0])

        self.model_part.CreateNewNode(5, 2.0, 0.0, 0.0)
        self.model_part.CreateNewNode(6, 2.0, 1.0, 0.0)
        new_element = self.model_part.CreateNewElement("Element2D4N", 2, [2,5,6,3], props)
        new_element.SetValue(KM.FORCE, [16.0, 4.0, 0.0])

        elemental_data = KM.Parameters("""{
            "model_part_name" : "default",
            "location"        : "element",
            "variable_name"   : "FORCE",
            "dimension"       : 3
        }""")

        self.interface_data = CouplingInterfaceData(elemental_data, self.model)

        self.solver_wrappers = {"dummy_solver" : DummySolverWrapper({"elemental_data" : self.interface_data})}

        self.solver_process_info = KM.ProcessInfo()

        conversion_op_settings = KM.Parameters("""{
            "type"           : "elemental_data_to_nodal_data",
            "solver"         : "dummy_solver",
            "data_name"      : "elemental_data",
            "echo_level"     : 0
        }""")

        conversion_operation = coupling_operation_factory.CreateCouplingOperation(conversion_op_settings, self.solver_wrappers, self.solver_process_info, KM.Testing.GetDefaultDataCommunicator())

        conversion_operation.Check()

        conversion_operation.Execute()

        nodal_data_output_setting = KM.Parameters("""{
            "model_part_name"       : "default",
            "variable_name"         : "FORCE",
            "location"              : "node_historical",
            "dimension"             : 3
        }""")

        nodal_data_output = CouplingInterfaceData(nodal_data_output_setting,  self.model)

        expected_nodal_values = [3, 2 ,0, 7, 3, 0, 7, 3 ,0, 3, 2 ,0, 4, 1 ,0, 4, 1 ,0 ]

        self.assertVectorAlmostEqual(expected_nodal_values, nodal_data_output.GetData())

    def test_elemental_to_nodal_conversion(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("default")
        self.model_part.AddNodalSolutionStepVariable(KM.HEAT_FLUX)
        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
        props = self.model_part.CreateNewProperties(1)

        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        self.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        self.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        new_element = self.model_part.CreateNewElement("Element2D4N", 1, [1,2,3,4], props)
        new_element.SetValue(KM.HEAT_FLUX, 12.0)

        self.model_part.CreateNewNode(5, 2.0, 0.0, 0.0)
        self.model_part.CreateNewNode(6, 2.0, 1.0, 0.0)
        new_element = self.model_part.CreateNewElement("Element2D4N", 2, [2,5,6,3], props)
        new_element.SetValue(KM.HEAT_FLUX, 16.0)

        elemental_data = KM.Parameters("""{
            "model_part_name" : "default",
            "location"        : "element",
            "variable_name"   : "HEAT_FLUX"
        }""")

        self.interface_data = CouplingInterfaceData(elemental_data, self.model)

        self.solver_wrappers = {"dummy_solver" : DummySolverWrapper({"elemental_data" : self.interface_data})}

        self.solver_process_info = KM.ProcessInfo()

        conversion_op_settings = KM.Parameters("""{
            "type"           : "elemental_data_to_nodal_data",
            "solver"         : "dummy_solver",
            "data_name"      : "elemental_data",
            "echo_level"     : 0
        }""")

        conversion_operation = coupling_operation_factory.CreateCouplingOperation(conversion_op_settings, self.solver_wrappers, self.solver_process_info, KM.Testing.GetDefaultDataCommunicator())

        conversion_operation.Check()

        conversion_operation.Execute()

        nodal_data_output_setting = KM.Parameters("""{
            "model_part_name"       : "default",
            "variable_name"         : "HEAT_FLUX",
            "location"              : "node_historical"
        }""")

        nodal_data_output = CouplingInterfaceData(nodal_data_output_setting,  self.model)

        expected_nodal_values = [3, 7, 7, 3, 4, 4]

        self.assertVectorAlmostEqual(expected_nodal_values, nodal_data_output.GetData())

    def test_nodal_to_elemental_conversion_scalar(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("default")
        self.model_part.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
        props = self.model_part.CreateNewProperties(1)

        new_node = self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        new_node.SetSolutionStepValue(KM.TEMPERATURE, 5)
        new_node = self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        new_node.SetSolutionStepValue(KM.TEMPERATURE, 4)
        new_node = self.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        new_node.SetSolutionStepValue(KM.TEMPERATURE, 4)
        new_node = self.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
        new_node.SetSolutionStepValue(KM.TEMPERATURE, 5)

        new_element = self.model_part.CreateNewElement("Element2D4N", 1, [1,2,3,4], props)

        new_node = self.model_part.CreateNewNode(5, 2.0, 0.0, 0.0)
        new_node.SetSolutionStepValue(KM.TEMPERATURE, 1)
        new_node = self.model_part.CreateNewNode(6, 2.0, 1.0, 0.0)
        new_node.SetSolutionStepValue(KM.TEMPERATURE, 1)
        new_element = self.model_part.CreateNewElement("Element2D4N", 2, [2,5,6,3], props)

        nodal_data = KM.Parameters("""{
            "model_part_name" : "default",
            "location"        : "node_historical",
            "variable_name"   : "TEMPERATURE"
        }""")

        self.nodal_data = CouplingInterfaceData(nodal_data, self.model)

        self.solver_wrappers = {"dummy_solver" : DummySolverWrapper({"nodal_data" : self.nodal_data})}

        self.solver_process_info = KM.ProcessInfo()

        conversion_op_settings = KM.Parameters("""{
            "type"           : "nodal_data_to_elemental_data",
            "solver"         : "dummy_solver",
            "data_name"      : "nodal_data",
            "echo_level"     : 0
        }""")

        conversion_operation = coupling_operation_factory.CreateCouplingOperation(conversion_op_settings, self.solver_wrappers, self.solver_process_info, KM.Testing.GetDefaultDataCommunicator())

        conversion_operation.Check()

        conversion_operation.Execute()

        elemental_data_output_setting = KM.Parameters("""{
            "model_part_name"       : "default",
            "variable_name"         : "TEMPERATURE",
            "location"              : "element"
        }""")

        elemental_data_output = CouplingInterfaceData(elemental_data_output_setting,  self.model)

        expected_elemental_values = [4.5, 2.5]
        expected_nodal_values = [5.0, 4.0, 4.0, 5.0, 1.0, 1.0]

        self.assertVectorAlmostEqual(expected_elemental_values, elemental_data_output.GetData())
        self.assertVectorAlmostEqual(expected_nodal_values, self.nodal_data.GetData())

    def test_nodal_to_elemental_conversion_vector(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("default")
        self.model_part.AddNodalSolutionStepVariable(KM.FORCE)
        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
        props = self.model_part.CreateNewProperties(1)

        new_node = self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        new_node.SetSolutionStepValue(KM.FORCE, [5, 4, 2])
        new_node = self.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        new_node.SetSolutionStepValue(KM.FORCE, [5, 4, 2])
        new_node = self.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        new_node.SetSolutionStepValue(KM.FORCE, [4, 3, 5])
        new_node = self.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
        new_node.SetSolutionStepValue(KM.FORCE, [4, 3, 5])

        new_element = self.model_part.CreateNewElement("Element2D4N", 1, [1,2,3,4], props)

        new_node = self.model_part.CreateNewNode(5, 2.0, 0.0, 0.0)
        new_node.SetSolutionStepValue(KM.FORCE, [1, 2, 10])
        new_node = self.model_part.CreateNewNode(6, 2.0, 1.0, 0.0)
        new_node.SetSolutionStepValue(KM.FORCE, [1, 2, 10])
        new_element = self.model_part.CreateNewElement("Element2D4N", 2, [2,5,6,3], props)

        nodal_data = KM.Parameters("""{
            "model_part_name" : "default",
            "location"        : "node_historical",
            "variable_name"   : "FORCE",
            "dimension"       : 3
        }""")

        self.nodal_data = CouplingInterfaceData(nodal_data, self.model)

        self.solver_wrappers = {"dummy_solver" : DummySolverWrapper({"nodal_data" : self.nodal_data})}

        self.solver_process_info = KM.ProcessInfo()

        conversion_op_settings = KM.Parameters("""{
            "type"           : "nodal_data_to_elemental_data",
            "solver"         : "dummy_solver",
            "data_name"      : "nodal_data",
            "echo_level"     : 0
        }""")

        conversion_operation = coupling_operation_factory.CreateCouplingOperation(conversion_op_settings, self.solver_wrappers, self.solver_process_info, KM.Testing.GetDefaultDataCommunicator())

        conversion_operation.Check()

        conversion_operation.Execute()

        elemental_data_output_setting = KM.Parameters("""{
            "model_part_name"       : "default",
            "variable_name"         : "FORCE",
            "location"              : "element",
            "dimension"             : 3
        }""")

        elemental_data_output = CouplingInterfaceData(elemental_data_output_setting,  self.model)

        expected_elemental_values = [4.5, 3.5, 3.5, 2.75, 2.75, 6.75]
        expected_nodal_values = [5, 4, 2, 5, 4, 2, 4, 3, 5, 4, 3, 5, 1, 2, 10, 1, 2, 10]

        self.assertVectorAlmostEqual(expected_elemental_values, elemental_data_output.GetData())
        self.assertVectorAlmostEqual(expected_nodal_values, self.nodal_data.GetData())

if __name__ == '__main__':
    KratosUnittest.main()