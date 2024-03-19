import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
from KratosMultiphysics.CoSimulationApplication.convergence_accelerators.convergence_accelerator_wrapper import BlockConvergenceAcceleratorWrapper
from testing_utilities import DummySolverWrapper

from unittest.mock import Mock, patch
import numpy as np
from random import uniform

if KM.IsDistributedRun():
    import KratosMultiphysics.mpi as KratosMPI

class TestBlockConvergenceAcceleratorWrapper(KratosUnittest.TestCase):

    def setUp(self):
        self.dimension = 3

        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("default")
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)
        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = self.dimension

        self.model2 = KM.Model()
        self.model_part2 = self.model2.CreateModelPart("default2")
        self.model_part2.AddNodalSolutionStepVariable(KM.MESH_DISPLACEMENT_X)
        self.model_part2.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)
        self.model_part2.ProcessInfo[KM.DOMAIN_SIZE] = self.dimension

        self.my_pid = KM.Testing.GetDefaultDataCommunicator().Rank()
        self.num_nodes = self.my_pid % 5 + 3 # num_nodes in range (3 ... 7)
        self.num_nodes2 = self.my_pid % 5 + 3 # num_nodes in range (3 ... 7)

        # in order to emulate one partition not having local nodes
        if self.my_pid == 4:
            self.num_nodes = 0
        if self.my_pid == 3:
            self.num_nodes2 = 0

        for i in range(self.num_nodes):
            node = self.model_part.CreateNewNode(i, 0.1*i, 0.0, 0.0) # this creates the same coords in different ranks, which does not matter for this test
            node.SetSolutionStepValue(KM.PARTITION_INDEX, self.my_pid)
            node.SetSolutionStepValue(KM.PRESSURE, uniform(-10, 50))

        for i in range(self.num_nodes2):
            node = self.model_part2.CreateNewNode(i, 0.1*i, 0.0, 0.0) # this creates the same coords in different ranks, which does not matter for this test
            node.SetSolutionStepValue(KM.PARTITION_INDEX, self.my_pid)
            node.SetSolutionStepValue(KM.MESH_DISPLACEMENT_X, uniform(-10, 50))

        if KM.IsDistributedRun():
            KratosMPI.ParallelFillCommunicator(self.model_part, KM.Testing.GetDefaultDataCommunicator()).Execute()
            KratosMPI.ParallelFillCommunicator(self.model_part2, KM.Testing.GetDefaultDataCommunicator()).Execute()

        data_settings = KM.Parameters("""{
            "model_part_name" : "default",
            "variable_name"   : "PRESSURE"
        }""")
        self.interface_data = CouplingInterfaceData(data_settings, self.model)
        self.dummy_solver_wrapper = DummySolverWrapper({"data_4_testing" : self.interface_data})

        data_settings2 = KM.Parameters("""{
            "model_part_name" : "default2",
            "variable_name"   : "MESH_DISPLACEMENT_X"
        }""")
        self.interface_data2 = CouplingInterfaceData(data_settings2, self.model2)
        self.dummy_solver_wrapper2 = DummySolverWrapper({"data_4_testing2" : self.interface_data2})

        self.interface_data_dict = {}
        self.interface_data_dict["data_4_testing"] = self.dummy_solver_wrapper.interface_data_dict["data_4_testing"]
        self.interface_data_dict["data_4_testing2"] = self.dummy_solver_wrapper2.interface_data_dict["data_4_testing2"]

    def test_accelerator_without_support_for_distributed_data(self):
        conv_acc_settings = KM.Parameters("""{
            "type"      : "patched_mock_testing",
                "solver_sequence" :
                [
                    {
                        "solver": "fluid",
                        "data_name"  : "data_4_testing",
                        "coupled_data_name" : "data_4_testing2"
                    },
                    {
                        "solver": "structure",
                        "data_name"  : "data_4_testing2",
                        "coupled_data_name" : "data_4_testing"
                    }
                ]
            }""")

        exp_inp = {}
        exp_res = {}

        exp_inp["data_4_testing"] = self.interface_data.GetData()
        exp_inp["data_4_testing2"] = self.interface_data2.GetData()

        update_solution_return_value = [uniform(-10, 50) for _ in range(self.num_nodes)]
        global_update_solution_return_value = np.array(np.concatenate(KM.Testing.GetDefaultDataCommunicator().GathervDoubles(update_solution_return_value, 0)))
        conv_acc_mock = Mock()
        attrs = {
            'SupportsDistributedData.return_value': False,
            'UpdateSolution.return_value' : global_update_solution_return_value
        }
        conv_acc_mock.configure_mock(**attrs)

        with patch('KratosMultiphysics.CoSimulationApplication.convergence_accelerators.convergence_accelerator_wrapper.CreateConvergenceAccelerator') as p:
            p.return_value = conv_acc_mock
            conv_acc_wrapper = BlockConvergenceAcceleratorWrapper(conv_acc_settings, self.interface_data_dict, KM.Testing.GetDefaultDataCommunicator())

        conv_acc_wrapper.Initialize()
        conv_acc_wrapper.InitializeSolutionStep()

        self.assertEqual(conv_acc_mock.SupportsDistributedData.call_count, 2)
        self.assertEqual(conv_acc_wrapper.gather_scatter_required["data_4_testing"], self.interface_data_dict["data_4_testing"].IsDistributed()) # gather-scatter is only required in case of distributed data
        self.assertEqual(conv_acc_wrapper.gather_scatter_required["data_4_testing2"], self.interface_data_dict["data_4_testing2"].IsDistributed()) # gather-scatter is only required in case of distributed data
        self.assertEqual(conv_acc_wrapper.executing_rank["data_4_testing"], self.my_pid == 0)
        self.assertEqual(conv_acc_wrapper.executing_rank["data_4_testing2"], self.my_pid == 0)

        conv_acc_wrapper.InitializeNonLinearIteration()

        # setting new solution for computing the residual
        rand_data1 = [uniform(-10, 50) for _ in range(self.num_nodes)]
        self.interface_data_dict["data_4_testing"].SetData(rand_data1)
        exp_res["data_4_testing"] = rand_data1 - exp_inp["data_4_testing"]
        exp_res_y = self.num_nodes2 * [0.]
        conv_acc_wrapper.ComputeAndApplyUpdate()


        global_exp_res = np.array(np.concatenate(KM.Testing.GetDefaultDataCommunicator().GathervDoubles(exp_res["data_4_testing"], 0)))
        global_exp_inp = np.array(np.concatenate(KM.Testing.GetDefaultDataCommunicator().GathervDoubles(exp_inp["data_4_testing"], 0)))
        global_exp_inp_y = np.array(np.concatenate(KM.Testing.GetDefaultDataCommunicator().GathervDoubles(exp_inp["data_4_testing2"], 0)))
        global_exp_res_y = np.array(np.concatenate(KM.Testing.GetDefaultDataCommunicator().GathervDoubles(exp_res_y, 0)))
        if self.my_pid == 0:
            np.testing.assert_array_equal(global_exp_res, conv_acc_mock.UpdateSolution.call_args[0][0])
            np.testing.assert_array_equal(global_exp_inp, conv_acc_mock.UpdateSolution.call_args[0][1])
            np.testing.assert_array_equal(global_exp_inp_y, conv_acc_mock.UpdateSolution.call_args[0][2])
            np.testing.assert_array_equal(global_exp_res_y, conv_acc_mock.UpdateSolution.call_args[0][4])

        solver_1_output = exp_inp["data_4_testing"] + update_solution_return_value
        np.testing.assert_array_equal(solver_1_output, self.interface_data_dict["data_4_testing"].GetData())


        update_solution_return_value2 = [uniform(-10, 50) for _ in range(self.num_nodes2)]
        global_update_solution_return_value2 = np.array(np.concatenate(KM.Testing.GetDefaultDataCommunicator().GathervDoubles(update_solution_return_value2, 0)))
        conv_acc_mock.UpdateSolution.return_value = global_update_solution_return_value2.copy()

        # setting new solution (other solver) for computing the residual
        rand_data2 = [uniform(-10, 50) for _ in range(self.num_nodes2)]
        self.interface_data_dict["data_4_testing2"].SetData(rand_data2)
        exp_res["data_4_testing2"] = rand_data2 - exp_inp["data_4_testing2"]
        exp_res_y =  solver_1_output - rand_data1
        conv_acc_wrapper.ComputeAndApplyUpdate()

        self.assertEqual(conv_acc_mock.UpdateSolution.call_count, 2*int(self.my_pid == 0)) # only one rank calls "UpdateSolution"
        global_exp_res2 = np.array(np.concatenate(KM.Testing.GetDefaultDataCommunicator().GathervDoubles(exp_res["data_4_testing2"], 0)))
        global_exp_inp2 = np.array(np.concatenate(KM.Testing.GetDefaultDataCommunicator().GathervDoubles(exp_inp["data_4_testing2"], 0)))
        global_exp_inp2_y = np.array(np.concatenate(KM.Testing.GetDefaultDataCommunicator().GathervDoubles(solver_1_output, 0)))
        global_exp_res2_y = np.array(np.concatenate(KM.Testing.GetDefaultDataCommunicator().GathervDoubles(exp_res_y, 0)))
        if self.my_pid == 0:
            np.testing.assert_array_equal(global_exp_res2, conv_acc_mock.UpdateSolution.call_args[0][0])
            np.testing.assert_array_equal(global_exp_inp2, conv_acc_mock.UpdateSolution.call_args[0][1])
            np.testing.assert_array_equal(global_exp_inp2_y, conv_acc_mock.UpdateSolution.call_args[0][2])
            np.testing.assert_array_equal(global_exp_res2_y, conv_acc_mock.UpdateSolution.call_args[0][4])

        np.testing.assert_array_equal(exp_inp["data_4_testing2"] + update_solution_return_value2, self.interface_data_dict["data_4_testing2"].GetData())

    def test_accelerator_with_support_for_distributed_data(self):
        conv_acc_settings = KM.Parameters("""{
            "type"      : "patched_mock_testing",
                "solver_sequence" :
                [
                    {
                        "solver": "fluid",
                        "data_name"  : "data_4_testing",
                        "coupled_data_name" : "data_4_testing2"
                    },
                    {
                        "solver": "structure",
                        "data_name"  : "data_4_testing2",
                        "coupled_data_name" : "data_4_testing"
                    }
                ]
            }""")

        exp_inp = {}
        exp_res = {}

        exp_inp["data_4_testing"] = self.interface_data.GetData()
        exp_inp["data_4_testing2"] = self.interface_data2.GetData()

        update_solution_return_value = [uniform(-10, 50) for _ in range(self.num_nodes)]
        conv_acc_mock = Mock()
        attrs = {
            'SupportsDistributedData.return_value': True,
            'UpdateSolution.return_value' : update_solution_return_value
        }
        conv_acc_mock.configure_mock(**attrs)

        with patch('KratosMultiphysics.CoSimulationApplication.convergence_accelerators.convergence_accelerator_wrapper.CreateConvergenceAccelerator') as p:
            p.return_value = conv_acc_mock
            conv_acc_wrapper = BlockConvergenceAcceleratorWrapper(conv_acc_settings, self.interface_data_dict, KM.Testing.GetDefaultDataCommunicator())

        conv_acc_wrapper.Initialize()
        conv_acc_wrapper.InitializeSolutionStep()

        self.assertEqual(conv_acc_mock.SupportsDistributedData.call_count, 2)
        self.assertFalse(conv_acc_wrapper.gather_scatter_required["data_4_testing"])
        self.assertFalse(conv_acc_wrapper.gather_scatter_required["data_4_testing2"])
        self.assertTrue(conv_acc_wrapper.executing_rank["data_4_testing"])
        self.assertTrue(conv_acc_wrapper.executing_rank["data_4_testing2"])

        conv_acc_wrapper.InitializeNonLinearIteration()

        # setting new solution for computing the residual
        rand_data1 = [uniform(-10, 50) for _ in range(self.num_nodes)]
        self.interface_data_dict["data_4_testing"].SetData(rand_data1)
        exp_res["data_4_testing"] = rand_data1 - exp_inp["data_4_testing"]
        exp_res_y = self.num_nodes2 * [0.]
        conv_acc_wrapper.ComputeAndApplyUpdate()


        np.testing.assert_array_equal(exp_res["data_4_testing"], conv_acc_mock.UpdateSolution.call_args[0][0])
        np.testing.assert_array_equal(exp_inp["data_4_testing"], conv_acc_mock.UpdateSolution.call_args[0][1])
        np.testing.assert_array_equal(exp_inp["data_4_testing2"], conv_acc_mock.UpdateSolution.call_args[0][2])
        np.testing.assert_array_equal(exp_res_y, conv_acc_mock.UpdateSolution.call_args[0][4])

        solver_1_output = exp_inp["data_4_testing"] + update_solution_return_value
        np.testing.assert_array_equal(solver_1_output, self.interface_data_dict["data_4_testing"].GetData())


        update_solution_return_value2 = [uniform(-10, 50) for _ in range(self.num_nodes2)]
        conv_acc_mock.UpdateSolution.return_value = update_solution_return_value2.copy()

        # setting new solution (other solver) for computing the residual
        rand_data2 = [uniform(-10, 50) for _ in range(self.num_nodes2)]
        self.interface_data_dict["data_4_testing2"].SetData(rand_data2)
        exp_res["data_4_testing2"] = rand_data2 - exp_inp["data_4_testing2"]
        exp_res_y =  solver_1_output - rand_data1
        conv_acc_wrapper.ComputeAndApplyUpdate()

        self.assertEqual(conv_acc_mock.UpdateSolution.call_count, 2) # only one rank calls "UpdateSolution"
        np.testing.assert_array_equal(exp_res["data_4_testing2"], conv_acc_mock.UpdateSolution.call_args[0][0])
        np.testing.assert_array_equal(exp_inp["data_4_testing2"], conv_acc_mock.UpdateSolution.call_args[0][1])
        np.testing.assert_array_equal(solver_1_output, conv_acc_mock.UpdateSolution.call_args[0][2])
        np.testing.assert_array_equal(exp_res_y, conv_acc_mock.UpdateSolution.call_args[0][4])

        np.testing.assert_array_equal(exp_inp["data_4_testing2"] + update_solution_return_value2, self.interface_data_dict["data_4_testing2"].GetData())


if __name__ == '__main__':
    KratosUnittest.main()
