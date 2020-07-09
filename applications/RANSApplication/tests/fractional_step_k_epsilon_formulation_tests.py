import KratosMultiphysics.KratosUnittest as UnitTest

import evm_turbulence_modelling_test_case


class FractionalStepKEpsilonTest(
        evm_turbulence_modelling_test_case.EvmTurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(FractionalStepKEpsilonTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_fractional_step_k_epsilon_parameters.json",
            False)


if __name__ == '__main__':
    UnitTest.main()