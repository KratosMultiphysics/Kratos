from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
from testing_utilities import DummySolverWrapper

from KratosMultiphysics.CoSimulationApplication.factories.convergence_criterion_factory import CreateConvergenceCriterion


class TestConvergenceCriteria(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("default")
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)

        for i in range(5):
            new_node = self.model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0) # position of nodes does not matter for this test

        data_settings = KM.Parameters("""{
            "model_part_name" : "default",
            "variable_name"   : "PRESSURE"
        }""")
        self.interface_data = CouplingInterfaceData(data_settings, self.model)
        self.interface_data.Initialize()

        self.solver_wrappers = {"dummy_solver" : DummySolverWrapper({"data_4_testing" : self.interface_data})}

    def test_RelativeNormInitialResidual_abs_tol(self):
        conv_crit_settings = KM.Parameters("""{
            "type"           : "relative_norm_initial_residual",
            "abs_tolerance"  : 1e-5,
            "rel_tolerance"  : 1e-5,
            "echo_level"     : 0
        }""")
        conv_crit = CreateConvergenceCriterion(conv_crit_settings)
        self.__ExecuteTest(conv_crit)

    def test_RelativeNormInitialResidual_rel_tol(self):
        conv_crit_settings = KM.Parameters("""{
            "type"           : "relative_norm_initial_residual",
            "abs_tolerance"  : 1e-5,
            "rel_tolerance"  : 1e-5,
            "echo_level"     : 0
        }""")
        conv_crit = CreateConvergenceCriterion(conv_crit_settings)
        self.__ExecuteTest(conv_crit)

    def test_RelativeNormPreviousResidual_abs_tol(self):
        conv_crit_settings = KM.Parameters("""{
            "type"           : "relative_norm_previous_residual",
            "abs_tolerance"  : 1e-5,
            "rel_tolerance"  : 1e-5,
            "echo_level"     : 0
        }""")
        conv_crit = CreateConvergenceCriterion(conv_crit_settings)
        self.__ExecuteTest(conv_crit)

    def test_RelativeNormPreviousResidual_rel_tol(self):
        conv_crit_settings = KM.Parameters("""{
            "type"           : "relative_norm_previous_residual",
            "abs_tolerance"  : 1e-5,
            "rel_tolerance"  : 1e-5,
            "echo_level"     : 0
        }""")
        conv_crit = CreateConvergenceCriterion(conv_crit_settings)
        self.__ExecuteTest(conv_crit)


    def __ExecuteTest(self, conv_crit):
        self.__SetInitialValues()
        conv_crit.Initialize()
        conv_crit.Check() # not yet implemented

        '''
        TS 1: no convergence
        TS 2: res decreases, abs
        TS 3: res decreases, rel
        TS 4: res does not change, abs
        TS 5: res does not change, rel
        '''
        time_steps = 5
        max_inner_iter = 5

        for i in range(time_steps):
            conv_crit.InitializeSolutionStep()

            for j in range(max_inner_iter):

                conv_crit.InitializeNonLinearIteration()

                self.__SetValues(1.0)

                # self.assertEqual(exp_conv_status, conv_crit.IsConverged())

                conv_crit.FinalizeNonLinearIteration()

            conv_crit.FinalizeSolutionStep()

        conv_crit.Finalize()

    def __SetInitialValues(self):
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(KM.PRESSURE, 1.0)

    def __SetValues(self, factor):
        pass


if __name__ == '__main__':
    KratosUnittest.main()
