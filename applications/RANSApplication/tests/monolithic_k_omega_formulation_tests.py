import KratosMultiphysics.KratosUnittest as UnitTest

import turbulence_modelling_test_case
import periodic_turbulence_modelling_test_case

class MonolithicKOmegaTest(turbulence_modelling_test_case.TurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(MonolithicKOmegaTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_monolithic_k_omega_parameters.json",
            False)

        cls.transient_scheme_type = "bossak"

class MonolithicKOmegaPeriodicTest(periodic_turbulence_modelling_test_case.PeriodicTurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(MonolithicKOmegaPeriodicTest, cls).setUpCase(
            "ChannelFlowTest",
            "channel_flow_monolithic_k_omega_parameters.json",
            False)

if __name__ == '__main__':
    UnitTest.main()