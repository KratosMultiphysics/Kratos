from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
from KratosMultiphysics.CoSimulationApplication.convergence_accelerators.convergence_accelerator_wrapper import ConvergenceAcceleratorWrapper
from testing_utilities import DummySolverWrapper

from unittest.mock import Mock
import numpy as np
from random import uniform

if KM.IsDistributedRun():
    import KratosMultiphysics.mpi as KratosMPI

class TestConvergenceAcceleratorWrapper(KratosUnittest.TestCase):

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

    def test_accelerator_without_support_for_distributed_data(self):
        conv_acc_settings = KM.Parameters("""{
            "type"      : "constant_relaxation",
            "data_name" : "data_4_testing"
        }""")
        conv_acc_wrapper = ConvergenceAcceleratorWrapper(conv_acc_settings, self.dummy_solver_wrapper)

        exp_inp = self.interface_data.GetData()
        update_solution_return_value = [uniform(-10, 50) for _ in range(self.num_nodes)]

        global_update_solution_return_value = np.array(np.concatenate(KM.DataCommunicator.GetDefault().GathervDoubles(update_solution_return_value, 0)))

        conv_acc_mock = Mock()

        attrs = {
            'SupportsDistributedData.return_value': False,
            'UpdateSolution.return_value' : global_update_solution_return_value
        }
        conv_acc_mock.configure_mock(**attrs)

        conv_acc_wrapper.conv_acc = conv_acc_mock

        conv_acc_wrapper.InitializeSolutionStep()

        self.assertEqual(conv_acc_mock.SupportsDistributedData.call_count, 1)
        self.assertEqual(conv_acc_wrapper.gather_scatter_required, self.interface_data.IsDistributed()) # gather-scatter is only required in case of distributed data
        self.assertEqual(conv_acc_wrapper.executing_rank, self.my_pid == 0)

        conv_acc_wrapper.InitializeNonLinearIteration()

        # setting new solution for computing the residual
        rand_data = [uniform(-10, 50) for _ in range(self.num_nodes)]
        self.interface_data.SetData(rand_data)
        exp_res = rand_data - exp_inp

        conv_acc_wrapper.ComputeAndApplyUpdate()

        self.assertEqual(conv_acc_mock.UpdateSolution.call_count, int(self.my_pid == 0)) # only one rank calls "UpdateSolution"
        global_exp_res = np.array(np.concatenate(KM.DataCommunicator.GetDefault().GathervDoubles(exp_res, 0)))
        global_exp_inp = np.array(np.concatenate(KM.DataCommunicator.GetDefault().GathervDoubles(exp_inp, 0)))
        if self.my_pid == 0:
            # numpy arrays cannot be compared using the mock-functions, hence using the numpy functions
            np.testing.assert_array_equal(global_exp_res, conv_acc_mock.UpdateSolution.call_args[0][0])
            np.testing.assert_array_equal(global_exp_inp, conv_acc_mock.UpdateSolution.call_args[0][1])

        np.testing.assert_array_equal(exp_inp + update_solution_return_value, self.interface_data.GetData())


    def test_accelerator_with_support_for_distributed_data(self):
        conv_acc_settings = KM.Parameters("""{
            "type"      : "constant_relaxation",
            "data_name" : "data_4_testing"
        }""")
        conv_acc_wrapper = ConvergenceAcceleratorWrapper(conv_acc_settings, self.dummy_solver_wrapper)

        exp_inp = self.interface_data.GetData()
        update_solution_return_value = [uniform(-10, 50) for _ in range(self.num_nodes)]

        conv_acc_mock = Mock()

        attrs = {
            'SupportsDistributedData.return_value': True,
            'UpdateSolution.return_value' : update_solution_return_value
        }
        conv_acc_mock.configure_mock(**attrs)

        conv_acc_wrapper.conv_acc = conv_acc_mock

        conv_acc_wrapper.InitializeSolutionStep()

        self.assertEqual(conv_acc_mock.SupportsDistributedData.call_count, 1)
        self.assertFalse(conv_acc_wrapper.gather_scatter_required)
        self.assertTrue(conv_acc_wrapper.executing_rank)

        conv_acc_wrapper.InitializeNonLinearIteration()

        # setting new solution for computing the residual
        rand_data = [uniform(-10, 50) for _ in range(self.num_nodes)]
        self.interface_data.SetData(rand_data)
        exp_res = rand_data - exp_inp

        conv_acc_wrapper.ComputeAndApplyUpdate()

        self.assertEqual(conv_acc_mock.UpdateSolution.call_count, 1)

        # numpy arrays cannot be compared using the mock-functions, hence using the numpy functions
        np.testing.assert_array_equal(exp_res, conv_acc_mock.UpdateSolution.call_args[0][0])
        np.testing.assert_array_equal(exp_inp, conv_acc_mock.UpdateSolution.call_args[0][1])

        np.testing.assert_array_equal(exp_inp + update_solution_return_value, self.interface_data.GetData())


if __name__ == '__main__':
    KratosUnittest.main()
