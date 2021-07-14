import KratosMultiphysics.KratosUnittest as UnitTest

import turbulence_modelling_test_case
import periodic_turbulence_modelling_test_case

class MonolithicKOmegaSSTTest(turbulence_modelling_test_case.TurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(MonolithicKOmegaSSTTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_mon_kwsst_parameters.json",
            "backward_facing_step_material_properties.json",
            False)

        cls.transient_scheme_type = "bossak"
        cls.parameters["<CONSTITUTIVE_LAW>"] = "RansKOmegaSSTNewtonian2DLaw"

class MonolithicKOmegaSSTPeriodicTest(periodic_turbulence_modelling_test_case.PeriodicTurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(MonolithicKOmegaSSTPeriodicTest, cls).setUpCase(
            "ChannelFlowTest",
            "channel_flow_mon_kwsst_parameters.json",
            "channel_flow_material_properties.json",
            False)

        cls.parameters["<CONSTITUTIVE_LAW>"] = "RansKOmegaSSTNewtonian2DLaw"


if __name__ == '__main__':
    UnitTest.main()
