import KratosMultiphysics.KratosUnittest as UnitTest

import turbulence_modelling_test_case


class FractionalStepKOmegaSSTTest(turbulence_modelling_test_case.TurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(FractionalStepKOmegaSSTTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_fractional_step_k_omega_sst_parameters.json",
            False)

        cls.transient_scheme_type = "bdf2"

if __name__ == '__main__':
    UnitTest.main()