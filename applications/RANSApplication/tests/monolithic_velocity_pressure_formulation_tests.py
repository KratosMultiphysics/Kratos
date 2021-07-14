import KratosMultiphysics.KratosUnittest as UnitTest

import flow_solver_test_case


class MonolithicVelocityPressureFormulationTest(flow_solver_test_case.FlowSolverTestCase):
    @classmethod
    def setUpClass(cls):
        super(MonolithicVelocityPressureFormulationTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_mon_up_parameters.json",
            "backward_facing_step_material_properties.json",
            False)

        cls.transient_scheme_type = "bossak"
        cls.parameters["<CONSTITUTIVE_LAW>"] = "Newtonian2DLaw"


if __name__ == '__main__':
    UnitTest.main()
