import KratosMultiphysics.KratosUnittest as UnitTest

import evm_turbulence_modelling_test_case


class FractionalStepKOmegaSSTTest(
        evm_turbulence_modelling_test_case.EvmTurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(FractionalStepKOmegaSSTTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_fractional_step_k_omega_sst_parameters.json",
            False)


if __name__ == '__main__':
    UnitTest.main()