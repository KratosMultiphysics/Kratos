import KratosMultiphysics.KratosUnittest as UnitTest

import evm_turbulence_modelling_test_case


class MonolithicKEpsilonHighReTest(
        evm_turbulence_modelling_test_case.EvmTurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(MonolithicKEpsilonHighReTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_monolithic_k_epsilon_high_re_parameters.json",
            False)


if __name__ == '__main__':
    UnitTest.main()