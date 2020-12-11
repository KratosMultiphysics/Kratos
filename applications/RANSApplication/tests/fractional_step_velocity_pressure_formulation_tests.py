import KratosMultiphysics.KratosUnittest as UnitTest

import flow_solver_test_case


class FractionalStepVelocityPressureFormulationTest(flow_solver_test_case.FlowSolverTestCase):
    @classmethod
    def setUpClass(cls):
        super(FractionalStepVelocityPressureFormulationTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_fractional_step_velocity_pressure_parameters.json",
            False)
        cls.transient_scheme_type = "bdf2"

if __name__ == '__main__':
    UnitTest.main()