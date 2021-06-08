import KratosMultiphysics.KratosUnittest as UnitTest

import turbulence_modelling_test_case


class FractionalStepKEpsilonTest(turbulence_modelling_test_case.TurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(FractionalStepKEpsilonTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_fs_ke_parameters.json",
            False)

        cls.transient_scheme_type = "bdf2"


if __name__ == '__main__':
    UnitTest.main()