import KratosMultiphysics.KratosUnittest as UnitTest

import flow_solver_test_case


class MonolithicVelocityPressureFormulationTest(
        flow_solver_test_case.FlowSolverTestCase):
    @classmethod
    def setUpClass(cls):
        super(MonolithicVelocityPressureFormulationTest, cls).setUpCase(
            "BackwardFacingStepTest",
            "backward_facing_step_monolithic_velocity_pressure_parameters.json",
            False)


if __name__ == '__main__':
    UnitTest.main()
