from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
from KratosMultiphysics.CoSimulationApplication.factories import convergence_criterion_factory
from testing_utilities import DummySolverWrapper


class TestConvergenceCriteria(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("default")
        self.model_part.AddNodalSolutionStepVariable(KM.PRESSURE)

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

    def test_RelativeNormInitialResidual_abs_tol(self):
        conv_crit_settings = KM.Parameters("""{
            "type"           : "relative_norm_initial_residual",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "abs_tolerance"  : 1e-5,
            "rel_tolerance"  : 1e-5,
            "echo_level"     : 0
        }""")
        conv_crit = convergence_criterion_factory.CreateConvergenceCriterion(conv_crit_settings, self.solver_wrappers["dummy_solver"])
        self.__ExecuteTest(conv_crit)

    def test_RelativeNormInitialResidual_rel_tol(self):
        conv_crit_settings = KM.Parameters("""{
            "type"           : "relative_norm_initial_residual",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "abs_tolerance"  : 1e-5,
            "rel_tolerance"  : 1e-5,
            "echo_level"     : 0
        }""")
        conv_crit = convergence_criterion_factory.CreateConvergenceCriterion(conv_crit_settings, self.solver_wrappers["dummy_solver"])
        self.__ExecuteTest(conv_crit)

    def test_RelativeNormPreviousResidual_abs_tol(self):
        conv_crit_settings = KM.Parameters("""{
            "type"           : "relative_norm_previous_residual",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "abs_tolerance"  : 1e-5,
            "rel_tolerance"  : 1e-5,
            "echo_level"     : 0
        }""")
        conv_crit = convergence_criterion_factory.CreateConvergenceCriterion(conv_crit_settings, self.solver_wrappers["dummy_solver"])
        self.__ExecuteTest(conv_crit)

    def test_RelativeNormPreviousResidual_rel_tol(self):
        conv_crit_settings = KM.Parameters("""{
            "type"           : "relative_norm_previous_residual",
            "solver"         : "dummy_solver",
            "data_name"      : "data_4_testing",
            "abs_tolerance"  : 1e-5,
            "rel_tolerance"  : 1e-5,
            "echo_level"     : 0
        }""")
        conv_crit = convergence_criterion_factory.CreateConvergenceCriterion(conv_crit_settings, self.solver_wrappers["dummy_solver"])
        self.__ExecuteTest(conv_crit)

    def __ExecuteTest(self, conv_crit):
        conv_crit.Initialize()
        conv_crit.Check()

        time_steps = 3
        max_inner_iter = 5

        for i in range(time_steps):
            conv_crit.InitializeSolutionStep()

            for j in range(max_inner_iter):

                conv_crit.InitializeNonLinearIteration()

                SetValues(self.model_part)

                # self.assertEqual(exp_conv_status, conv_crit.IsConverged())

                conv_crit.FinalizeNonLinearIteration()

            conv_crit.FinalizeSolutionStep()

        conv_crit.Finalize()

def SetValues(model_part):
    pass


if __name__ == '__main__':
    KratosUnittest.main()
