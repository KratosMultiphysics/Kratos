from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
from testing_utilities import DummySolverWrapper

from KratosMultiphysics.CoSimulationApplication.factories.convergence_criterion_factory import CreateConvergenceCriterion
from KratosMultiphysics.CoSimulationApplication.convergence_criteria.convergence_criteria_wrapper import ConvergenceCriteriaWrapper

from unittest.mock import Mock
import numpy as np
from random import uniform

if KM.IsDistributedRun():
    import KratosMultiphysics.mpi as KratosMPI

class TestConvergenceCriteriaWrapper(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("default")
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)
        self.dimension = 3
        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = self.dimension

        self.my_pid = KM.DataCommunicator.GetDefault().Rank()
        self.num_nodes = self.my_pid % 5 + 3 # num_nodes in range (3 ... 7)
        if self.my_pid == 4:
            self.num_nodes = 0 # in order to emulate one partition not having local nodes

        for i in range(self.num_nodes):
            node = self.model_part.CreateNewNode(i, 0.1*i, 0.0, 0.0) # this creates the same coords in different ranks, which does not matter for this test

            node.SetSolutionStepValue(KM.PARTITION_INDEX, self.my_pid)
            node.SetSolutionStepValue(KM.PRESSURE, uniform(-10, 50))

        if KM.IsDistributedRun():
            KratosMPI.ParallelFillCommunicator(self.model_part).Execute()

        data_settings = KM.Parameters("""{
            "model_part_name" : "default",
            "variable_name"   : "PRESSURE"
        }""")
        self.interface_data = CouplingInterfaceData(data_settings, self.model)
        self.interface_data.Initialize()

        self.dummy_solver_wrapper = DummySolverWrapper({"data_4_testing" : self.interface_data})

    def test_wrapper(self):
        conv_acc_settings = KM.Parameters("""{
            "type"      : "relative_norm_previous_residual",
            "data_name" : "data_4_testing"
        }""")
        conv_crit_wrapper = ConvergenceCriteriaWrapper(conv_acc_settings, self.dummy_solver_wrapper)

        data_init = self.interface_data.GetData()

        conv_crit_mock = Mock()

        is_converged = True
        attrs = {
            'IsConverged.return_value' : is_converged
        }
        conv_crit_mock.configure_mock(**attrs)

        conv_crit_wrapper.conv_crit = conv_crit_mock

        conv_crit_wrapper.InitializeSolutionStep()

        conv_crit_wrapper.InitializeNonLinearIteration()

        # setting new solution for computing the residual
        rand_data = [uniform(-10, 50) for _ in range(self.num_nodes)]
        self.interface_data.SetData(rand_data)
        exp_res = rand_data - data_init

        self.assertEqual(conv_crit_wrapper.IsConverged(), is_converged)

        self.assertEqual(conv_crit_mock.IsConverged.call_count, int(self.my_pid == 0)) # only one rank calls "IsConverged"

        global_exp_res = np.array(np.concatenate(KM.DataCommunicator.GetDefault().GathervDoubles(exp_res, 0)))
        global_rand_data_inp = np.array(np.concatenate(KM.DataCommunicator.GetDefault().GathervDoubles(rand_data, 0)))
        if self.my_pid == 0:
            # numpy arrays cannot be compared using the mock-functions, hence using the numpy functions
            np.testing.assert_array_equal(global_exp_res, conv_crit_mock.IsConverged.call_args[0][0])
            np.testing.assert_array_equal(global_rand_data_inp, conv_crit_mock.IsConverged.call_args[0][1])


class TestConvergenceCriteria(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("default")
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)

        for i in range(10): # using 10 nodes gives suitable values in the tests
            node = self.model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0) # position of nodes does not matter for this test
            node.SetSolutionStepValue(KM.PRESSURE, 1.0)

        data_settings = KM.Parameters("""{
            "model_part_name" : "default",
            "variable_name"   : "PRESSURE"
        }""")
        self.interface_data = CouplingInterfaceData(data_settings, self.model)
        self.interface_data.Initialize()

        self.dummy_solver_wrapper = DummySolverWrapper({"data_4_testing" : self.interface_data})

    def test_RelativeNormInitialResidual_abs_tol(self):
        conv_crit_settings = KM.Parameters("""{
            "type"           : "relative_norm_initial_residual",
            "data_name"      : "data_4_testing",
            "abs_tolerance"  : 1e-5,
            "rel_tolerance"  : 1e-12,
            "echo_level"     : 0
        }""")
        conv_crit = ConvergenceCriteriaWrapper(conv_crit_settings, self.dummy_solver_wrapper)

        sol_values = [
            (2e-1, False),
            (2e-4, False),
            (2e-6, False),
            (2e-7, True),
            (2e-8, True),
            (2e-9, True),
            (2e-4, False),
            (2e-6, False),
            (2e-7, True)
        ]

        self.__ExecuteTest(conv_crit, sol_values)

    def test_RelativeNormInitialResidual_rel_tol(self):
        conv_crit_settings = KM.Parameters("""{
            "type"           : "relative_norm_initial_residual",
            "data_name"      : "data_4_testing",
            "abs_tolerance"  : 1e-12,
            "rel_tolerance"  : 1e-5,
            "echo_level"     : 0
        }""")
        conv_crit = ConvergenceCriteriaWrapper(conv_crit_settings, self.dummy_solver_wrapper)

        sol_values = [
            (2e-1, False),
            (2e-2, False),
            (2e-3, False),
            (2e-4, False),
            (2e-5, False),
            (2e-6, False),
            (1e-6, True),
            (1e-7, True),
            (2e-5, False),
            (2e-6, False),
            (1e-6, True),
        ]

        self.__ExecuteTest(conv_crit, sol_values)

    def test_RelativeNormPreviousResidual_abs_tol(self):
        conv_crit_settings = KM.Parameters("""{
            "type"           : "relative_norm_previous_residual",
            "data_name"      : "data_4_testing",
            "abs_tolerance"  : 1e-5,
            "rel_tolerance"  : 1e-12,
            "echo_level"     : 0
        }""")
        conv_crit = ConvergenceCriteriaWrapper(conv_crit_settings, self.dummy_solver_wrapper)

        sol_values = [
            (2e-1, False),
            (2e-4, False),
            (2e-6, False),
            (2e-7, True),
            (2e-8, True),
            (2e-9, True),
            (2e-4, False),
            (2e-6, False),
            (2e-7, True)
        ]

        self.__ExecuteTest(conv_crit, sol_values)

    def test_RelativeNormPreviousResidual_rel_tol(self):
        conv_crit_settings = KM.Parameters("""{
            "type"           : "relative_norm_previous_residual",
            "data_name"      : "data_4_testing",
            "abs_tolerance"  : 1e-12,
            "rel_tolerance"  : 1e-5,
            "echo_level"     : 0
        }""")
        conv_crit = ConvergenceCriteriaWrapper(conv_crit_settings, self.dummy_solver_wrapper)

        sol_values = [
            (2e-1, False),
            (2e-2, False),
            (2.00001e-2, True), # small change in the values leads to convergence
            (2e-4, False),
            (8e-6, False),
            (7e-6, False),
            (7.00001e-6, True) # small change in the values leads to convergence
        ]

        self.__ExecuteTest(conv_crit, sol_values)


    def __ExecuteTest(self, conv_crit, solution_values):
        conv_crit.Initialize()
        conv_crit.Check() # not yet implemented

        conv_crit.InitializeSolutionStep()

        for vals_tuple in solution_values:

            conv_crit.InitializeNonLinearIteration()

            KM.VariableUtils().SetScalarVar(KM.PRESSURE, vals_tuple[0], self.model_part.Nodes)

            self.assertEqual(vals_tuple[1], conv_crit.IsConverged())

            conv_crit.FinalizeNonLinearIteration()

        conv_crit.FinalizeSolutionStep()

        conv_crit.Finalize()


if __name__ == '__main__':
    KratosUnittest.main()
