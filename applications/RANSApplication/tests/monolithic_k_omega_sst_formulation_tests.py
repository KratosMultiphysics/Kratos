import KratosMultiphysics.KratosUnittest as UnitTest

import evm_turbulence_modelling_test_case


class MonolithicKOmegaSSTTest(
        evm_turbulence_modelling_test_case.EvmTurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(MonolithicKOmegaSSTTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_monolithic_k_omega_sst_parameters.json",
            False)


if __name__ == '__main__':
    UnitTest.main()
