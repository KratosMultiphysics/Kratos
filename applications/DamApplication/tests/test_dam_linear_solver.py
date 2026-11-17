import os

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis


class _DisplacementControlledInterfaceAnalysis(DamAnalysis):
    displacement_history = (
        2.0e-4, 5.0e-4, 2.0e-5, -2.0e-5,
        -2.0e-4, -2.0e-5, 2.0e-5, 2.0e-4
    )
    boundary_node_ids = (7, 11, 12, 16, 20, 24)

    def Initialize(self):
        super().Initialize()
        for node_id in self.boundary_node_ids:
            self.main_model_part.GetNode(node_id).Fix(KM.DISPLACEMENT_X)

    def RunMainTemporalLoop(self):
        for step, prescribed_displacement in enumerate(self.displacement_history, 1):
            self.time += self.delta_time
            self.main_model_part.CloneTimeStep(self.time)
            self.main_model_part.ProcessInfo[KM.STEP] = step
            for node_id in self.boundary_node_ids:
                self.main_model_part.GetNode(node_id).SetSolutionStepValue(
                    KM.DISPLACEMENT_X, prescribed_displacement)
            for process in self.list_of_processes:
                process.ExecuteInitializeSolutionStep()
            self.solver.Solve()
            for process in self.list_of_processes:
                process.ExecuteFinalizeSolutionStep()

    def CaptureFinalState(self):
        process_info = self.main_model_part.ProcessInfo
        interface = self.main_model_part.GetElement(21)
        return {
            "displacements": {
                node.Id: list(node.GetSolutionStepValue(KM.DISPLACEMENT))
                for node in self.main_model_part.Nodes
            },
            "jump": list(interface.CalculateOnIntegrationPoints(
                KratosPoro.LOCAL_RELATIVE_DISPLACEMENT_VECTOR, process_info)[0]),
            "traction": list(interface.CalculateOnIntegrationPoints(
                KratosPoro.LOCAL_STRESS_VECTOR, process_info)[0]),
            "reaction": sum(
                self.main_model_part.GetNode(node_id).GetSolutionStepValue(KM.REACTION_X)
                for node_id in self.boundary_node_ids)
        }


class TestDamLinearSolver(KratosUnittest.TestCase):
    def _GetParameters(self, use_direct_solver):
        test_directory = os.path.dirname(os.path.realpath(__file__))
        input_directory = os.path.join(
            test_directory, "joint_elastic_cohesive_2d_normal")
        with open(os.path.join(
                input_directory,
                "joint_elastic_cohesive_2d_normal_parameters.json"), "r") as parameter_file:
            parameters = KM.Parameters(parameter_file.read())

        parameters["problem_data"]["end_time"].SetDouble(
            float(len(_DisplacementControlledInterfaceAnalysis.displacement_history)))
        parameters["problem_data"]["time_step"].SetDouble(1.0)
        parameters["problem_data"]["number_of_threads"].SetInt(1)
        parameters["solver_settings"]["echo_level"].SetInt(0)
        parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(
            os.path.join(input_directory, "joint_elastic_cohesive_2d_normal"))
        parameters["loads_process_list"][0]["Parameters"]["modulus"].SetString("0.0")
        if use_direct_solver:
            linear_settings = parameters["solver_settings"][
                "mechanical_solver_settings"]["linear_solver_settings"]
            linear_settings.RemoveValue("preconditioner_type")
            linear_settings.RemoveValue("max_iteration")
            linear_settings.RemoveValue("tolerance")
            linear_settings.RemoveValue("scaling")
            linear_settings["solver_type"].SetString(
                "skyline_lu_factorization")
        return parameters

    def _RunCase(self, use_direct_solver):
        analysis = _DisplacementControlledInterfaceAnalysis(
            KM.Model(), self._GetParameters(use_direct_solver))
        analysis.Run()
        return analysis.CaptureFinalState()

    def _AssertFinalStatesAlmostEqual(self, expected, actual):
        for node_id, expected_displacement in expected["displacements"].items():
            for expected_value, actual_value in zip(
                    expected_displacement, actual["displacements"][node_id]):
                self.assertAlmostEqual(actual_value, expected_value, delta=1.0e-12)
        for variable_name in ("jump", "traction"):
            for expected_value, actual_value in zip(
                    expected[variable_name], actual[variable_name]):
                self.assertAlmostEqual(actual_value, expected_value, delta=1.0e-10)
        self.assertAlmostEqual(actual["reaction"], expected["reaction"], delta=1.0e-8)

    def test_iterative_and_direct_solvers_produce_equivalent_solution(self):
        with KratosUnittest.WorkFolderScope(".", __file__):
            iterative_state = self._RunCase(False)
            direct_state = self._RunCase(True)
        self._AssertFinalStatesAlmostEqual(iterative_state, direct_state)


if __name__ == "__main__":
    KratosUnittest.main()
