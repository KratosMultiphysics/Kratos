import numpy

import KratosMultiphysics
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics import scipy_conversion_tools


class NumpyRomProjector:
    """Reference implementation in NumPy/SciPy of the C++ RomProjector (same interface)."""

    def __init__(self, scheme, strategy_data):
        self.scheme = scheme
        self.strategy_data = strategy_data

    def BuildEffectiveSystem(self):
        # Set the system arrays to zero
        linear_system = self.strategy_data.GetLinearSystem()
        lhs = linear_system.GetMatrix(KratosMultiphysics.Future.SparseMatrixTag.LHS)
        rhs = linear_system.GetVector(KratosMultiphysics.Future.DenseVectorTag.RHS)
        linear_system.GetVector(KratosMultiphysics.Future.DenseVectorTag.Dx).SetValue(0.0)
        lhs.SetValue(0.0)
        rhs.SetValue(0.0)

        # Build the system and apply the constraints
        self.scheme.Build(lhs, rhs)
        self.scheme.BuildLinearSystemConstraints(self.strategy_data)
        self.scheme.ApplyLinearSystemConstraints(self.strategy_data)

        # The rows of the basis of the fixed DOFs are set to zero in the projection
        self.free_dofs = numpy.array([not dof.IsFixed() for dof in self.strategy_data.GetEffectiveDofSet()])

    def Project(self, phi):
        effective_linear_system = self.strategy_data.GetEffectiveLinearSystem()
        lhs = scipy_conversion_tools.to_csr(effective_linear_system.GetMatrix(KratosMultiphysics.Future.SparseMatrixTag.LHS))
        rhs = numpy.array(effective_linear_system.GetVector(KratosMultiphysics.Future.DenseVectorTag.RHS), copy=False)

        phi = phi * self.free_dofs[:, None]
        self.reduced_lhs = phi.T @ (lhs @ phi)
        self.reduced_rhs = phi.T @ rhs
        return self.reduced_lhs, self.reduced_rhs

    def SolveReduced(self):
        return numpy.linalg.solve(self.reduced_lhs, self.reduced_rhs)

    def SetSolution(self, u):
        dof_set = self.strategy_data.GetEffectiveDofSet()
        effective_dx = self.strategy_data.GetEffectiveLinearSystem().GetVector(KratosMultiphysics.Future.DenseVectorTag.Dx)
        numpy.array(effective_dx, copy=False)[:] = u - numpy.array(dof_set.GetValues(), copy=False)
        self.scheme.Update(self.strategy_data)


class FutureRomSolver:
    """Projection-based ROM on top of the Kratos Future schemes.

    The decoder u = D(q) and its Jacobian are evaluated in Python. The projector builds the
    full order system through the Future scheme, projects it and solves the reduced system.
    """

    def __init__(self, model_part, settings=None):
        if not hasattr(KratosMultiphysics, "Future"):
            raise Exception("FutureRomSolver requires Kratos to be compiled with 'KRATOS_USE_FUTURE=ON'.")

        if settings is None:
            settings = KratosMultiphysics.Parameters("{}")
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model_part = model_part
        self.decoder = None
        self.max_iterations = settings["max_iterations"].GetInt()
        self.relative_tolerance = settings["relative_tolerance"].GetDouble()
        self.absolute_tolerance = settings["absolute_tolerance"].GetDouble()
        self.residual_relative_tolerance = settings["residual_relative_tolerance"].GetDouble()
        self.residual_absolute_tolerance = settings["residual_absolute_tolerance"].GetDouble()
        self.echo_level = settings["echo_level"].GetInt()
        self.convergence_history = []

        # Convergence criterion: "displacement_criterion" (increment of the free DOFs, "relative_tolerance" and "absolute_tolerance")
        # or "reduced_residual_criterion" (RHS of the projected system, "residual_relative_tolerance" and "residual_absolute_tolerance")
        self.convergence_criterion = settings["convergence_criterion"].GetString()
        if self.convergence_criterion not in ("displacement_criterion", "reduced_residual_criterion"):
            err_msg = "Unknown value \'{}\' for \'convergence_criterion\'. Available options are \'displacement_criterion\' and \'reduced_residual_criterion\'.".format(self.convergence_criterion)
            raise Exception(err_msg)

        # Values of the fixed DOFs given by the decoder: "basis" (u0 + Phi q, as the free ones) or "prescribed" (the values imposed
        # by the conditions, see FixedDofsDecoder). The projection does not use the rows of the fixed DOFs in either case
        self.fixed_dofs_decoding = settings["fixed_dofs_decoding"].GetString()
        if self.fixed_dofs_decoding not in ("basis", "prescribed"):
            err_msg = "Unknown value \'{}\' for \'fixed_dofs_decoding\'. Available options are \'basis\' and \'prescribed\'.".format(self.fixed_dofs_decoding)
            raise Exception(err_msg)

        self.scheme = KratosMultiphysics.Future.StaticScheme(self.model_part, settings["scheme_settings"])
        self.strategy_data = KratosMultiphysics.Future.ImplicitStrategyData()

        projection_backend = settings["projection_backend"].GetString()
        if projection_backend == "numpy":
            self.projector = NumpyRomProjector(self.scheme, self.strategy_data)
        elif projection_backend == "cpp":
            self.projector = KratosROM.Future.RomProjector(self.scheme, self.strategy_data)
        else:
            err_msg = "Unknown value \'{}\' for \'projection_backend\'. Available options are \'numpy\' and \'cpp\'.".format(projection_backend)
            raise Exception(err_msg)

    @classmethod
    def GetDefaultParameters(cls):
        return KratosMultiphysics.Parameters("""{
            "scheme_settings" : {
                "name" : "static_scheme",
                "build_settings" : {
                    "name" : "block_builder"
                },
                "echo_level" : 0,
                "move_mesh" : false,
                "reform_dofs_at_each_step" : false
            },
            "projection_backend" : "numpy",
            "fixed_dofs_decoding" : "basis",
            "max_iterations" : 10,
            "convergence_criterion" : "displacement_criterion",
            "relative_tolerance" : 1.0e-9,
            "absolute_tolerance" : 1.0e-9,
            "residual_relative_tolerance" : 1.0e-9,
            "residual_absolute_tolerance" : 1.0e-9,
            "echo_level" : 0
        }""")

    def Initialize(self):
        self.scheme.Initialize(self.strategy_data)

    def SetDecoder(self, decoder):
        """Sets the decoder. Its rows must follow the effective DOF set, which is available after Initialize."""
        self.decoder = decoder
        self.q = numpy.zeros(self.decoder.NumberOfRomDofs())

    def InitializeSolutionStep(self):
        self.scheme.InitializeSolutionStep(self.strategy_data)
        # The conditions of the step are already applied, so the decoder can read the prescribed values of the fixed DOFs
        self.decoder.InitializeSolutionStep()

    def SolveSolutionStep(self):
        """Newton-Raphson iteration on the reduced coordinates. Returns True if converged.
        With the "displacement_criterion" the convergence is checked after each solve on the increment of the free DOFs (J dq), as
        the legacy displacement criterion does. With the "reduced_residual_criterion" it is checked before each solve on the reduced
        residual (the RHS of the projected system), relative to the one of the first iteration of the step or absolute (divided by
        the square root of the number of ROM DOFs), as the legacy residual criterion does on the full residual. The step then ends
        without a further update, and the system is built once more after the last iteration to check it.
        The history of the iterations of the step is kept (see GetConvergenceHistory)."""
        self.convergence_history = []
        is_converged = False
        free_dofs = numpy.array([not dof.IsFixed() for dof in self.GetEffectiveDofSet()])
        check_residual = self.convergence_criterion == "reduced_residual_criterion"
        initial_residual_norm = None
        for iteration in range(1, self.max_iterations + (2 if check_residual else 1)):
            # u = D(q)
            self.projector.SetSolution(self.decoder.Decode(self.q))

            # Build and project the system evaluated at u = D(q)
            self.scheme.InitializeNonLinIteration(self.strategy_data)
            self.projector.BuildEffectiveSystem()
            phi = self.decoder.Jacobian(self.q)
            _, reduced_rhs = self.projector.Project(phi)
            residual_norm = float(numpy.linalg.norm(reduced_rhs))
            if initial_residual_norm is None:
                initial_residual_norm = residual_norm
            history = {"iteration": iteration, "reduced_residual_norm": residual_norm}
            self.convergence_history.append(history)

            if check_residual:
                ratio, absolute_norm = self._ResidualCriterionNorms(residual_norm, initial_residual_norm, len(reduced_rhs))
                history.update({"residual_ratio": ratio, "residual_absolute_norm": absolute_norm})
                is_converged = ratio <= self.residual_relative_tolerance or absolute_norm < self.residual_absolute_tolerance
                self._PrintIteration(history)
                if is_converged or iteration > self.max_iterations:
                    self.scheme.FinalizeNonLinIteration(self.strategy_data)
                    break

            # Solve and update the reduced coordinates
            dq = self.projector.SolveReduced()
            self.q = self.q + dq
            self.scheme.FinalizeNonLinIteration(self.strategy_data)

            ratio, absolute_norm = self._DisplacementCriterionNorms((phi @ dq)[free_dofs], self.decoder.Decode(self.q)[free_dofs])
            history.update({"reduced_increment_norm": float(numpy.linalg.norm(dq)), "displacement_ratio": ratio, "displacement_absolute_norm": absolute_norm})
            if not check_residual:
                is_converged = ratio <= self.relative_tolerance or absolute_norm < self.absolute_tolerance
                self._PrintIteration(history)
                if is_converged:
                    break

        # Leave the database at the last reduced coordinates
        self.projector.SetSolution(self.decoder.Decode(self.q))

        return is_converged

    def GetConvergenceHistory(self):
        """History of the iterations of the last step: one dict per build of the system with the "iteration", the norm of the reduced
        residual before the solve ("reduced_residual_norm") and, if a solve followed, the norm of the increment of the reduced
        coordinates ("reduced_increment_norm") and the displacement criterion values. With the "reduced_residual_criterion", also
        its values; its last entry, where it is checked without a further solve, has no increment."""
        return self.convergence_history

    def _PrintIteration(self, history):
        if self.echo_level > 0:
            KratosMultiphysics.Logger.PrintInfo("FutureRomSolver", ", ".join(f"{key} = {value}" for key, value in history.items()))

    @staticmethod
    def _ResidualCriterionNorms(residual_norm, initial_residual_norm, number_of_rom_dofs):
        """Ratio to the reduced residual of the first iteration of the step and absolute norm, as the legacy ResidualCriteria on the full residual"""
        epsilon = numpy.finfo(float).eps
        ratio = 0.0 if residual_norm < epsilon else residual_norm / max(initial_residual_norm, epsilon)
        return ratio, residual_norm / numpy.sqrt(max(number_of_rom_dofs, 1))

    @staticmethod
    def _DisplacementCriterionNorms(du, u):
        """Ratio and absolute norm of the correction du of the free DOFs, as computed by the legacy DisplacementCriteria"""
        correction_norm = numpy.linalg.norm(du)
        reference_norm = numpy.linalg.norm(u)
        epsilon = numpy.finfo(float).eps
        if reference_norm < epsilon:
            reference_norm = correction_norm
        ratio = 0.0 if correction_norm < epsilon else correction_norm / reference_norm
        return ratio, correction_norm / numpy.sqrt(max(len(du), 1))

    def FinalizeSolutionStep(self):
        self.scheme.FinalizeSolutionStep(self.strategy_data)

    def GetRomSolution(self):
        return self.q

    def SetRomSolution(self, q):
        self.q = numpy.array(q, dtype=numpy.float64)

    def GetEffectiveDofSet(self):
        return self.strategy_data.GetEffectiveDofSet()
