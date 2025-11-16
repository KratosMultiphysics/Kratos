import KratosMultiphysics.KratosUnittest as UnitTest

import turbulence_modelling_test_case


class FractionalStepKOmegaSSTTest(turbulence_modelling_test_case.TurbulenceModellingTestCase):
    @classmethod
    def setUpClass(cls):
        super(FractionalStepKOmegaSSTTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_fs_kwsst_parameters.json",
            "backward_facing_step_material_properties.json",
            False)

        cls.transient_scheme_type = "bdf2"
        cls.parameters["<CONSTITUTIVE_LAW>"] = "RansKOmegaSSTNewtonian2DLaw"

if __name__ == '__main__':
    UnitTest.main()