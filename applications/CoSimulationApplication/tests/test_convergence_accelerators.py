from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
from KratosMultiphysics.CoSimulationApplication.convergence_accelerators.convergence_accelerator_wrapper import ConvergenceAcceleratorWrapper
from testing_utilities import DummySolverWrapper

if KM.IsDistributedRun():
    import KratosMultiphysics.mpi as KratosMPI

class TestConvergenceAccelerators(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("default")
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)
        self.dimension = 3
        self.model_part.ProcessInfo[KM.DOMAIN_SIZE] = self.dimension

        num_proc = KM.DataCommunicator.GetDefault().Size()
        my_pid = KM.DataCommunicator.GetDefault().Rank()
        self.num_nodes = my_pid % 5 + 3 # num_nodes in range (3 ... 7)
        if my_pid == 1:
            self.num_nodes = 0 # in order to emulate one partition not having local nodes

        for i in range(self.num_nodes):
            node = self.model_part.CreateNewNode(i, 0.1*i, 0.0, 0.0) # this creates the same coords in different ranks, which does not matter for this test

            node.SetSolutionStepValue(KM.PARTITION_INDEX, my_pid)

        if KM.IsDistributedRun():
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(self.model_part)
            KratosMPI.ParallelFillCommunicator(self.model_part).Execute()

        print(my_pid, self.model_part)

        data_settings_vector = KM.Parameters("""{
            "model_part_name" : "default",
            "variable_name"   : "PRESSURE"
        }""")
        self.interface_data = CouplingInterfaceData(data_settings_vector, self.model)
        self.interface_data.Initialize()


        self.dummy_solver_wrapper = DummySolverWrapper({"data_4_testing" : self.interface_data})

    def test_convergence_accelerator_wrapper(self):
        conv_acc_settings = KM.Parameters("""{
            "type"      : "constant_relaxation",
            "data_name" : "data_4_testing"
        }""")
        conv_acc_wrapper = ConvergenceAcceleratorWrapper(conv_acc_settings, self.dummy_solver_wrapper)

        # TODO this is the order as is done in the FSI-App, we have to update this accordingly in CoSim
        conv_acc_wrapper.Initialize()
        conv_acc_wrapper.InitializeSolutionStep()
        conv_acc_wrapper.InitializeNonLinearIteration()
        conv_acc_wrapper.ComputeAndApplyUpdate()
        conv_acc_wrapper.FinalizeNonLinearIteration()
        conv_acc_wrapper.FinalizeSolutionStep()
        conv_acc_wrapper.Finalize()


if __name__ == '__main__':
    KratosUnittest.main()
