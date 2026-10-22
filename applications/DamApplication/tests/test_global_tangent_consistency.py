import math
import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.DamApplication as DamApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.PoromechanicsApplication as PoromechanicsApplication

from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis


class _DisplacementControlledAnalysis(DamAnalysis):
    boundary_node_ids = (7, 11, 12, 16, 20, 24)

    def Initialize(self):
        super().Initialize()
        for node_id in self.boundary_node_ids:
            self.main_model_part.GetNode(node_id).Fix(Kratos.DISPLACEMENT_X)


class _AssembledTangentCheck:
    interface_element_ids = range(21, 26)

    def __init__(self, model_part):
        self.model_part = model_part
        self.scheme = DamApplication.IncrementalUpdateStaticSmoothingScheme()
        self.space = Kratos.SparseSpace()
        self.builder = Kratos.ResidualBasedBlockBuilderAndSolver(
            Kratos.SkylineLUFactorizationSolver())
        self.A = self.space.CreateEmptyMatrixPointer()
        self.Dx = self.space.CreateEmptyVectorPointer()
        self.b = self.space.CreateEmptyVectorPointer()

        self.scheme.Initialize(model_part)
        self.scheme.InitializeElements(model_part)
        self.scheme.InitializeConditions(model_part)
        self.builder.SetUpDofSet(self.scheme, model_part)
        self.builder.SetUpSystem(model_part)
        self.builder.ResizeAndInitializeVectors(
            self.scheme, self.A, self.Dx, self.b, model_part)
        self.dofs = self.builder.GetDofSet()
        self.free_equation_ids = [
            dof.EquationId for dof in self.dofs if not dof.IsFixed()]

    @staticmethod
    def _norm(values):
        return math.sqrt(sum(value * value for value in values))

    def _move_mesh(self):
        for node in self.model_part.Nodes:
            displacement = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            node.X = node.X0 + displacement[0]
            node.Y = node.Y0 + displacement[1]
            node.Z = node.Z0 + displacement[2]

    def _jumps(self):
        jumps = []
        for element_id in self.interface_element_ids:
            values = self.model_part.GetElement(element_id).CalculateOnIntegrationPoints(
                PoromechanicsApplication.LOCAL_RELATIVE_DISPLACEMENT_VECTOR,
                self.model_part.ProcessInfo)
            jumps.extend(value[1] for value in values)
        return jumps

    def _snapshot(self):
        return (
            [dof.GetSolutionStepValue() for dof in self.dofs],
            [(node.X, node.Y, node.Z) for node in self.model_part.Nodes])

    def _restore(self, snapshot):
        dof_values, coordinates = snapshot
        for dof, value in zip(self.dofs, dof_values):
            dof.SetSolutionStepValue(value)
        for node, position in zip(self.model_part.Nodes, coordinates):
            node.X, node.Y, node.Z = position

    def _residual(self):
        residual = Kratos.Vector(len(self.b))
        self.space.SetToZeroVector(residual)
        self.builder.BuildRHS(self.scheme, self.model_part, residual)
        return residual

    def _perturb(self, snapshot, direction, amount):
        for dof, value in zip(self.dofs, snapshot[0]):
            if not dof.IsFixed():
                dof.SetSolutionStepValue(
                    value + amount * direction[dof.EquationId])
        self._move_mesh()

    def evaluate(self):
        self.space.SetToZeroMatrix(self.A)
        self.space.SetToZeroVector(self.b)
        self.builder.Build(self.scheme, self.model_part, self.A, self.b)

        direction = [0.0] * len(self.dofs)
        for dof in self.dofs:
            if not dof.IsFixed():
                equation_id = dof.EquationId
                direction[equation_id] = float((equation_id * 37) % 17 - 8)
        direction_norm = self._norm(direction)
        direction = [value / direction_norm for value in direction]

        tangent_product = Kratos.Vector(len(direction))
        self.space.Mult(self.A, Kratos.Vector(direction), tangent_product)
        tangent_product = [
            tangent_product[i] for i in self.free_equation_ids]

        snapshot = self._snapshot()
        h = 2.0e-11
        self._perturb(snapshot, direction, h)
        plus_jumps = self._jumps()
        plus_residual = self._residual()
        self._restore(snapshot)
        self._perturb(snapshot, direction, -h)
        minus_jumps = self._jumps()
        minus_residual = self._residual()
        self._restore(snapshot)

        switched_count = sum(
            (plus < 0.0) != (minus < 0.0)
            for plus, minus in zip(plus_jumps, minus_jumps))
        all_jumps = plus_jumps + minus_jumps
        if min(all_jumps) >= 0.0:
            classification = "OPENING_ONLY"
        elif max(all_jumps) < 0.0:
            classification = "COMPRESSION_ONLY"
        else:
            classification = "MIXED_NO_CROSS"

        derivative = [
            (plus_residual[i] - minus_residual[i]) / (2.0 * h)
            for i in self.free_equation_ids]
        error = [
            tangent_i + derivative_i
            for tangent_i, derivative_i in zip(tangent_product, derivative)]
        relative_error = self._norm(error) / max(
            self._norm(tangent_product), self._norm(derivative), 1.0e-30)
        return classification, switched_count, relative_error


class TestGlobalTangentConsistency(KratosUnittest.TestCase):
    displacement_history = (2.0e-4, 5.0e-4, 2.0e-5, -2.0e-5, -2.0e-4)

    @staticmethod
    def _parameters():
        case_directory = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "joint_elastic_cohesive_2d_normal")
        with open(os.path.join(
                case_directory,
                "joint_elastic_cohesive_2d_normal_parameters.json")) as parameter_file:
            parameters = Kratos.Parameters(parameter_file.read())
        parameters["problem_data"]["end_time"].SetDouble(5.0)
        parameters["problem_data"]["time_step"].SetDouble(1.0)
        parameters["problem_data"]["number_of_threads"].SetInt(1)
        parameters["solver_settings"]["echo_level"].SetInt(0)
        parameters["solver_settings"]["mechanical_solver_settings"][
            "move_mesh_flag"].SetBool(True)
        parameters["solver_settings"]["model_import_settings"][
            "input_filename"].SetString(os.path.join(
                case_directory, "joint_elastic_cohesive_2d_normal"))
        parameters["loads_process_list"][0]["Parameters"][
            "modulus"].SetString("0.0")
        return parameters

    def test_opening_and_mixed_state_assembled_tangents(self):
        analysis = _DisplacementControlledAnalysis(
            Kratos.Model(), self._parameters())
        previous_severity = Kratos.Logger.GetDefaultOutput().GetSeverity()
        Kratos.Logger.GetDefaultOutput().SetSeverity(Kratos.Logger.Severity.WARNING)
        try:
            analysis.Initialize()
            tangent_check = _AssembledTangentCheck(analysis.main_model_part)
            results = {}
            for step, prescribed_displacement in enumerate(
                    self.displacement_history, 1):
                analysis.time += analysis.delta_time
                analysis.main_model_part.CloneTimeStep(analysis.time)
                analysis.main_model_part.ProcessInfo[Kratos.STEP] = step
                for node_id in analysis.boundary_node_ids:
                    analysis.main_model_part.GetNode(node_id).SetSolutionStepValue(
                        Kratos.DISPLACEMENT_X, prescribed_displacement)
                for process in analysis.list_of_processes:
                    process.ExecuteInitializeSolutionStep()
                analysis.solver.Solve()
                if step in (2, 5):
                    results[step] = tangent_check.evaluate()
                for process in analysis.list_of_processes:
                    process.ExecuteFinalizeSolutionStep()

            self.assertEqual(results[2][0], "OPENING_ONLY")
            self.assertEqual(results[5][0], "MIXED_NO_CROSS")
            for classification, switched_count, relative_error in results.values():
                self.assertEqual(switched_count, 0, classification)
                self.assertLess(relative_error, 1.0e-6)
        finally:
            analysis.Finalize()
            Kratos.Logger.GetDefaultOutput().SetSeverity(previous_severity)


if __name__ == "__main__":
    KratosUnittest.main()
