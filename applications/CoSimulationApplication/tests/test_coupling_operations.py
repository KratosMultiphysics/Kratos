from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
from KratosMultiphysics.CoSimulationApplication.factories import coupling_operation_factory
from testing_utilities import DummySolverWrapper

from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import UsingPyKratos
using_pykratos = UsingPyKratos()

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
        self.interface_data.Initialize()

        self.solver_wrappers = {"dummy_solver" : DummySolverWrapper({"data_4_testing" : self.interface_data})}

    def test_constant_scaling(self):
        scaling_op_settings = KM.Parameters("""{
            "type"           : "scaling",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "scaling_factor" : 1.5,
            "echo_level"     : 0
        }""")

        scaling_op = coupling_operation_factory.CreateCouplingOperation(scaling_op_settings, self.solver_wrappers)

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

        scaling_op = coupling_operation_factory.CreateCouplingOperation(scaling_op_settings, self.solver_wrappers)

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

        scaling_op = coupling_operation_factory.CreateCouplingOperation(scaling_op_settings, self.solver_wrappers)

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

        scaling_op = coupling_operation_factory.CreateCouplingOperation(scaling_op_settings, self.solver_wrappers)

        factors = [1.5*pi*sqrt(1), 1.5*pi*sqrt(2), 1.5*pi*sqrt(3), 1.5*pi*sqrt(4), 1.5*pi*sqrt(5)]

        self.__ExecuteTest(scaling_op, factors)

    def test_scaling_in_interval(self):
        if using_pykratos:
            self.skipTest("This test can only be run with pyKratos after the IntervalUtility is implemented!")

        scaling_op_settings = KM.Parameters("""{
            "type"           : "scaling",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "scaling_factor" : 1.22,
            "interval"       : [0.0, 0.3]
        }""")

        scaling_op = coupling_operation_factory.CreateCouplingOperation(scaling_op_settings, self.solver_wrappers)

        factors = [1.0] * 5
        factors[0] = 1.22

        self.__ExecuteTest(scaling_op, factors)

    def test_scaling_in_interval_2(self):
        if using_pykratos:
            self.skipTest("This test can only be run with pyKratos after the IntervalUtility is implemented!")

        scaling_op_settings = KM.Parameters("""{
            "type"           : "scaling",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "scaling_factor" : 1.22,
            "interval"       : [0.8, "End"]
        }""")

        scaling_op = coupling_operation_factory.CreateCouplingOperation(scaling_op_settings, self.solver_wrappers)

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


if __name__ == '__main__':
    KratosUnittest.main()
