import KratosMultiphysics.KratosUnittest as UnitTest

import turbulence_modelling_test_case


class FractionalStepKOmegaTest(turbulence_modelling_test_case.TurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(FractionalStepKOmegaTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_fs_kw_parameters.json",
            False)

        cls.transient_scheme_type = "bdf2"


if __name__ == '__main__':
    UnitTest.main()