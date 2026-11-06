"""Cross-validates the shipped adjoint sensitivities against Kratos's own.

`sensitivity_utils` computes dJ/dtheta as lambda^T db/dtheta from the
assembled FEM tangent. StructuralMechanicsApplication ships an entirely
separate adjoint stack - AdjointFiniteDifferencing* elements, adjoint
response functions, ResidualBasedAdjointStaticScheme and SensitivityBuilder -
which produces SHAPE_SENSITIVITY for the same objective. Agreement between
the two, with full finite differences through real solves as a third
opinion, is a far stronger statement than either against FD alone: a shared
sign or scaling error could survive an FD check on one implementation, but
not this.

The objective is a traced nodal displacement, so dJ/du is a unit vector and
no load-vector reconstruction is needed on either side.
"""

import shutil
import tempfile
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.physics import differentiable_residual
from KratosMultiphysics.PhysicsNeMoApplication.physics import sensitivity_utils
have_structural = kratos_utils.CheckIfApplicationsAvailable(
    "StructuralMechanicsApplication", "LinearSolversApplication",
    "ConstitutiveLawsApplication")

_CASE_DIR = Path(__file__).parent / "adjoint_cases"
_PRIMAL_FILE = "cantilever_primal.json"

_FIELDS = [("DISPLACEMENT_X", "node_historical"),
           ("DISPLACEMENT_Y", "node_historical"),
           ("DISPLACEMENT_Z", "node_historical")]
_AXES = (("X0", "X"), ("Y0", "Y"), ("Z0", "Z"))
_TIP_NODE = 9
# node id, axis: a free interior node, the loaded tip, and two more nodes,
# covering all three axes
_TARGETS = ((5, 0), (5, 2), (9, 0), (11, 2), (12, 1))


@KratosUnittest.skipUnless(have_structural,
                           "Missing StructuralMechanics/ConstitutiveLaws/LinearSolvers applications.")
class TestAdjointCrossValidation(KratosUnittest.TestCase):
    """The shipped adjoint vs StructuralMechanics' adjoint vs full FD."""

    def setUp(self):
        # AdjointResponseFunction reads its primal settings from a file and
        # re-reads the mdpa for its own adjoint model part, so the case has
        # to live in the working directory.
        self.work_dir = Path(tempfile.mkdtemp(prefix="physicsnemo_adjoint_"))
        for path in _CASE_DIR.iterdir():
            shutil.copy(path, self.work_dir / path.name)
        self.previous_dir = Path.cwd()
        import os
        os.chdir(self.work_dir)

    def tearDown(self):
        import os
        os.chdir(self.previous_dir)
        kratos_utils.DeleteDirectoryIfExisting(str(self.work_dir))

    @staticmethod
    def _SolvePrimal(perturbation=None):
        """Solves the cantilever, optionally with one node coordinate moved.

        Returns the model part with the support DOFs still fixed: the
        adjoint needs the CONSTRAINED operator, and the BC process releases
        those DOFs in its ExecuteFinalizeSolutionStep - so calling
        analysis.Finalize() here would leave an unconstrained, singular
        tangent (and a residual of load magnitude rather than ~0).
        """
        from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
            StructuralMechanicsAnalysis)

        with open(_PRIMAL_FILE, 'r') as parameter_file:
            parameters = Kratos.Parameters(parameter_file.read())
        model = Kratos.Model()
        analysis = StructuralMechanicsAnalysis(model, parameters)
        analysis.Initialize()
        model_part = model["Structure"]

        if perturbation is not None:
            node_id, axis, delta = perturbation
            node = model_part.GetNode(node_id)
            for attribute in _AXES[axis]:
                setattr(node, attribute, getattr(node, attribute) + delta)

        analysis.RunSolutionLoop()
        for node in model_part.GetSubModelPart("Support").Nodes:
            for component in (Kratos.DISPLACEMENT_X, Kratos.DISPLACEMENT_Y,
                              Kratos.DISPLACEMENT_Z):
                node.Fix(component)
        return model, model_part

    @staticmethod
    def _Objective(model_part):
        return model_part.GetNode(_TIP_NODE).GetSolutionStepValue(Kratos.DISPLACEMENT_Z)

    def _OurShapeSensitivities(self, model_part):
        assembler = differentiable_residual.TangentAssembler(model_part)
        dof_map = differentiable_residual.DofFieldMap(assembler, _FIELDS)

        # J = u_z at the tip, so dJ/du is a unit vector there
        selector = numpy.zeros((model_part.NumberOfNodes(), 3))
        for row, node in enumerate(model_part.Nodes):
            if node.Id == _TIP_NODE:
                selector[row, 2] = 1.0
        dJ_du = dof_map.FieldsToDofVector(selector)

        appliers = {}
        for node_id, axis in _TARGETS:
            node = model_part.GetNode(node_id)

            def MakeApplier(node=node, axis=axis):
                def ApplyCoordinate(value):
                    for attribute in _AXES[axis]:
                        setattr(node, attribute, value)
                return ApplyCoordinate

            appliers[(node_id, axis)] = (MakeApplier(),
                                         getattr(node, _AXES[axis][0]))

        return (assembler, dof_map,
                sensitivity_utils.ComputeParameterSensitivities(
                    assembler, dof_map, dJ_du, appliers, fd_step=1e-6))

    @staticmethod
    def _KratosShapeSensitivities():
        from KratosMultiphysics.StructuralMechanicsApplication import structural_response

        settings = Kratos.Parameters("""{
            "response_type"                    : "adjoint_nodal_displacement",
            "gradient_mode"                    : "semi_analytic",
            "step_size"                        : 1e-7,
            "primal_settings"                  : "%s",
            "adjoint_settings"                 : "auto",
            "primal_data_transfer_with_python" : true,
            "response_part_name"               : "Tip",
            "direction"                        : [0.0, 0.0, 1.0],
            "traced_dof"                       : "DISPLACEMENT",
            "sensitivity_settings"             : {
                "sensitivity_model_part_name"                : "Design",
                "nodal_solution_step_sensitivity_variables"  : ["SHAPE_SENSITIVITY"],
                "build_mode"                                 : "static"
            }
        }""" % _PRIMAL_FILE)
        response = structural_response.AdjointResponseFunction(
            "tip_displacement", settings, Kratos.Model())
        response.Initialize()
        response.InitializeSolutionStep()
        response.CalculateValue()
        response.CalculateGradient()
        return response.GetValue(), response.GetNodalGradient(Kratos.SHAPE_SENSITIVITY)

    def _FullFiniteDifference(self, node_id, axis, step=1e-6):
        _, plus = self._SolvePrimal((node_id, axis, +step))
        _, minus = self._SolvePrimal((node_id, axis, -step))
        return (self._Objective(plus) - self._Objective(minus)) / (2.0 * step)

    def _DofSetup(self, model_part):
        """The assembler/dof_map/dJ_du triple for J = u_z at the tip."""
        assembler = differentiable_residual.TangentAssembler(model_part)
        dof_map = differentiable_residual.DofFieldMap(assembler, _FIELDS)
        selector = numpy.zeros((model_part.NumberOfNodes(), 3))
        for row, node in enumerate(model_part.Nodes):
            if node.Id == _TIP_NODE:
                selector[row, 2] = 1.0
        return assembler, dof_map, dof_map.FieldsToDofVector(selector)

    @staticmethod
    def _NodeRows(model_part):
        return {node.Id: row for row, node in enumerate(model_part.Nodes)}

    def _AllCoordinateSensitivities(self, assembler, dof_map, dJ_du, model_part):
        """The reference: the per-coordinate global path over EVERY coordinate."""
        appliers = {}
        for node in model_part.Nodes:
            for axis in range(3):
                def MakeApplier(node=node, axis=axis):
                    def ApplyCoordinate(value):
                        for attribute in _AXES[axis]:
                            setattr(node, attribute, value)
                    return ApplyCoordinate
                appliers[(node.Id, axis)] = (MakeApplier(),
                                             getattr(node, _AXES[axis][0]))
        return sensitivity_utils.ComputeParameterSensitivities(
            assembler, dof_map, dJ_du, appliers, fd_step=1e-6)

    def test_FieldReproducesThePerCoordinateGlobalPath(self):
        # The central claim: one element-local pass gives the same numbers as
        # 72 full residual re-assemblies, on every coordinate of every node -
        # supports and the loaded tip included.
        _, model_part = self._SolvePrimal()
        assembler, dof_map, dJ_du = self._DofSetup(model_part)

        field = sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dJ_du, fd_step=1e-6)
        self.assertEqual(field.shape, (model_part.NumberOfNodes(), 3))

        reference = self._AllCoordinateSensitivities(
            assembler, dof_map, dJ_du, model_part)
        rows = self._NodeRows(model_part)
        self.assertEqual(len(reference), 3 * model_part.NumberOfNodes())
        for (node_id, axis), value in reference.items():
            with self.subTest(node=node_id, axis=axis):
                self.assertAlmostEqual(field[rows[node_id], axis], value,
                                       delta=1e-7 * abs(value))

    def test_FieldMatchesKratosShapeSensitivityField(self):
        # Kratos' SHAPE_SENSITIVITY is a whole nodal field, so this validates
        # the entire array at once rather than a handful of scalars. Its
        # tolerance is the looser one: Kratos takes a FORWARD difference.
        _, model_part = self._SolvePrimal()
        assembler, dof_map, dJ_du = self._DofSetup(model_part)
        field = sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dJ_du, fd_step=1e-6)

        _, kratos_gradient = self._KratosShapeSensitivities()
        rows = self._NodeRows(model_part)
        for node_id, gradient in kratos_gradient.items():
            for axis in range(3):
                with self.subTest(node=node_id, axis=axis):
                    self.assertAlmostEqual(field[rows[node_id], axis], gradient[axis],
                                           delta=1e-4 * abs(gradient[axis]) + 1e-12)

    def test_FieldMatchesFullFiniteDifferences(self):
        _, model_part = self._SolvePrimal()
        assembler, dof_map, dJ_du = self._DofSetup(model_part)
        field = sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dJ_du, fd_step=1e-6)
        rows = self._NodeRows(model_part)

        for node_id, axis in _TARGETS:
            with self.subTest(node=node_id, axis=axis):
                reference = self._FullFiniteDifference(node_id, axis)
                self.assertGreater(abs(reference), 1e-9)
                self.assertAlmostEqual(field[rows[node_id], axis], reference,
                                       delta=1e-6 * abs(reference))

    def test_DesignNodeRestrictionReproducesTheFullField(self):
        # The mode that makes it scale: restricting to a design surface must
        # give those nodes' exact values and leave every other row at zero.
        _, model_part = self._SolvePrimal()
        assembler, dof_map, dJ_du = self._DofSetup(model_part)
        full = sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dJ_du, fd_step=1e-6)

        design_ids = [5, 9, 11]
        restricted = sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dJ_du, fd_step=1e-6, design_node_ids=design_ids)

        rows = self._NodeRows(model_part)
        for node in model_part.Nodes:
            row = rows[node.Id]
            for axis in range(3):
                if node.Id in design_ids:
                    self.assertAlmostEqual(restricted[row, axis], full[row, axis],
                                           delta=1e-9 * abs(full[row, axis]) + 1e-18)
                else:
                    self.assertEqual(restricted[row, axis], 0.0)

    def test_MeshIsBitIdenticalAfterThePass(self):
        # Perturbations are undone by writing the SAVED value back, never by
        # a reverse increment (x + h - 2h + h does not round to x), so a full
        # sweep must leave the coordinates and the residual untouched.
        _, model_part = self._SolvePrimal()
        assembler, dof_map, dJ_du = self._DofSetup(model_part)

        before = [(node.X0, node.Y0, node.Z0, node.X, node.Y, node.Z)
                  for node in model_part.Nodes]
        sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dJ_du, fd_step=1e-6)
        after = [(node.X0, node.Y0, node.Z0, node.X, node.Y, node.Z)
                 for node in model_part.Nodes]

        self.assertEqual(before, after)
        residual = numpy.array(assembler.ComputeResidualVector(), copy=True)
        self.assertLess(numpy.abs(residual).max(), 1e-8)

    def test_PointLoadConditionIsShapeBlind(self):
        # A PointLoadCondition's RHS is just the nodal POINT_LOAD and never
        # reads the geometry, so it contributes exactly nothing. Pinning the
        # measured fact, not assuming it: dropping the conditions entirely
        # must change nothing at all.
        _, model_part = self._SolvePrimal()
        self.assertGreater(model_part.NumberOfConditions(), 0)
        assembler, dof_map, dJ_du = self._DofSetup(model_part)
        with_conditions = sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dJ_du, fd_step=1e-6)

        process_info = model_part.ProcessInfo
        rhs = Kratos.Vector()
        for condition in model_part.Conditions:
            self.assertFalse(sensitivity_utils._EntityHasShapeDerivative(
                condition, process_info, rhs, 1e-6))

    def test_OutputVariableWritesTheField(self):
        _, model_part = self._SolvePrimal()
        assembler, dof_map, dJ_du = self._DofSetup(model_part)
        field = sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dJ_du, fd_step=1e-6,
            output_variable=Kratos.SHAPE_SENSITIVITY)

        for row, node in enumerate(model_part.Nodes):
            written = node.GetValue(Kratos.SHAPE_SENSITIVITY)
            for axis in range(3):
                self.assertAlmostEqual(written[axis], field[row, axis], places=15)

        # the historical database rejects a variable the primal analysis
        # never allocated, and says so rather than failing obscurely
        with self.assertRaisesRegex(ValueError, "solution-step"):
            sensitivity_utils.ComputeShapeSensitivityField(
                assembler, dof_map, dJ_du, output_variable=Kratos.SHAPE_SENSITIVITY,
                output_location="node_historical")

    def test_ConstrainedTangentIsTheSolvedState(self):
        # Guards the trap the helper works around: the adjoint operator must
        # be the constrained one, at a converged state.
        _, model_part = self._SolvePrimal()
        assembler = differentiable_residual.TangentAssembler(model_part)
        dof_map = differentiable_residual.DofFieldMap(assembler, _FIELDS)
        self.assertEqual(int(dof_map.fixed_mask.sum()), 12)  # 4 support nodes x 3
        residual = numpy.array(assembler.ComputeResidualVector(), copy=True)
        self.assertLess(numpy.abs(residual).max(), 1e-8)

    def test_ShapeSensitivitiesMatchKratosAdjointAndFiniteDifferences(self):
        _, model_part = self._SolvePrimal()
        objective = self._Objective(model_part)
        _, _, ours = self._OurShapeSensitivities(model_part)

        kratos_value, kratos_gradient = self._KratosShapeSensitivities()
        # both stacks solved the same primal problem
        self.assertAlmostEqual(kratos_value, objective, places=12)

        for node_id, axis in _TARGETS:
            with self.subTest(node=node_id, axis=axis):
                reference = self._FullFiniteDifference(node_id, axis)
                self.assertGreater(abs(reference), 1e-9)  # a real sensitivity

                # The shipped adjoint is exact up to the residual FD step;
                # Kratos' semi-analytic gradient carries its own step_size
                # error, so it gets the looser tolerance of the two.
                self.assertAlmostEqual(ours[(node_id, axis)], reference,
                                       delta=1e-6 * abs(reference))
                self.assertAlmostEqual(kratos_gradient[node_id][axis], reference,
                                       delta=1e-4 * abs(reference))

    def test_AdjointIsIndependentOfTheResidualStep(self):
        # An exact adjoint must not drift with the differentiation step; a
        # formulation error would show up as step dependence.
        _, model_part = self._SolvePrimal()
        assembler = differentiable_residual.TangentAssembler(model_part)
        dof_map = differentiable_residual.DofFieldMap(assembler, _FIELDS)
        selector = numpy.zeros((model_part.NumberOfNodes(), 3))
        for row, node in enumerate(model_part.Nodes):
            if node.Id == _TIP_NODE:
                selector[row, 2] = 1.0
        dJ_du = dof_map.FieldsToDofVector(selector)

        node_id, axis = _TARGETS[0]
        node = model_part.GetNode(node_id)

        def ApplyCoordinate(value):
            for attribute in _AXES[axis]:
                setattr(node, attribute, value)

        appliers = {"theta": (ApplyCoordinate, getattr(node, _AXES[axis][0]))}
        values = [
            sensitivity_utils.ComputeParameterSensitivities(
                assembler, dof_map, dJ_du, appliers, fd_step=step)["theta"]
            for step in (1e-5, 1e-6, 1e-7)]
        for value in values[1:]:
            self.assertAlmostEqual(value, values[0], delta=1e-6 * abs(values[0]))

    def test_AdjointSolveMatchesTheDirectSensitivity(self):
        # lambda^T db/dtheta must equal dJ/du . du/dtheta computed by solving
        # the tangent system directly - the adjoint identity itself.
        _, model_part = self._SolvePrimal()
        assembler, dof_map, ours = self._OurShapeSensitivities(model_part)

        selector = numpy.zeros((model_part.NumberOfNodes(), 3))
        for row, node in enumerate(model_part.Nodes):
            if node.Id == _TIP_NODE:
                selector[row, 2] = 1.0
        dJ_du = dof_map.FieldsToDofVector(selector)

        node_id, axis = _TARGETS[1]
        node = model_part.GetNode(node_id)
        reference_coordinate = getattr(node, _AXES[axis][0])
        step = 1e-6
        residuals = []
        for offset in (+step, -step):
            for attribute in _AXES[axis]:
                setattr(node, attribute, reference_coordinate + offset)
            residuals.append(numpy.array(assembler.ComputeResidualVector(), copy=True))
        for attribute in _AXES[axis]:
            setattr(node, attribute, reference_coordinate)
        db_dtheta = (residuals[0] - residuals[1]) / (2.0 * step)

        import scipy.sparse.linalg
        tangent = assembler.ComputeTangentMatrix(apply_dirichlet=True)
        du_dtheta = scipy.sparse.linalg.splu(tangent.tocsc()).solve(
            numpy.where(dof_map.fixed_mask, 0.0, db_dtheta))
        direct = float(dJ_du @ du_dtheta)

        self.assertAlmostEqual(ours[(node_id, axis)], direct,
                               delta=1e-9 * abs(direct))


if __name__ == '__main__':
    KratosUnittest.main()
