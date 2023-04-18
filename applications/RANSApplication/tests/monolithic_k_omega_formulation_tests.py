import KratosMultiphysics.KratosUnittest as UnitTest

import turbulence_modelling_test_case
import periodic_turbulence_modelling_test_case

class MonolithicKOmegaTest(turbulence_modelling_test_case.TurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(MonolithicKOmegaTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_mon_kw_parameters.json",
            "backward_facing_step_material_properties.json",
            False)

        cls.transient_scheme_type = "bossak"
        cls.parameters["<CONSTITUTIVE_LAW>"] = "RansKOmegaNewtonian2DLaw"

class MonolithicKOmegaPeriodicTest(periodic_turbulence_modelling_test_case.PeriodicTurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(MonolithicKOmegaPeriodicTest, cls).setUpCase(
            "ChannelFlowTest",
            "channel_flow_mon_kw_parameters.json",
            "channel_flow_material_properties.json",
            False)

        cls.parameters["<CONSTITUTIVE_LAW>"] = "RansKOmegaNewtonian2DLaw"

if __name__ == '__main__':
    UnitTest.main()