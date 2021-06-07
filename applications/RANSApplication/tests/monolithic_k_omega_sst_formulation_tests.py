import KratosMultiphysics.KratosUnittest as UnitTest

import turbulence_modelling_test_case
import periodic_turbulence_modelling_test_case

class MonolithicKOmegaSSTTest(turbulence_modelling_test_case.TurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(MonolithicKOmegaSSTTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_mon_kwsst_parameters.json",
            False)

        cls.transient_scheme_type = "bossak"

class MonolithicKOmegaSSTPeriodicTest(periodic_turbulence_modelling_test_case.PeriodicTurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(MonolithicKOmegaSSTPeriodicTest, cls).setUpCase(
            "ChannelFlowTest",
            "channel_flow_mon_kwsst_parameters.json",
            False)


if __name__ == '__main__':
    UnitTest.main()
