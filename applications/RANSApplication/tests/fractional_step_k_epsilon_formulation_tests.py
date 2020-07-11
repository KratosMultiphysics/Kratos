import KratosMultiphysics.KratosUnittest as UnitTest

import turbulence_modelling_test_case


class FractionalStepKEpsilonTest(
        turbulence_modelling_test_case.TurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(FractionalStepKEpsilonTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_fractional_step_k_epsilon_parameters.json",
            False)


if __name__ == '__main__':
    UnitTest.main()