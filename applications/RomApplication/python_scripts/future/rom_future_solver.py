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
        self.echo_level = settings["echo_level"].GetInt()

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
            "relative_tolerance" : 1.0e-9,
            "absolute_tolerance" : 1.0e-9,
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
        The convergence criterion is the one of the legacy displacement criterion, on the increment of the free DOFs (Phi dq)."""
        is_converged = False
        iteration = 0
        free_dofs = numpy.array([not dof.IsFixed() for dof in self.GetEffectiveDofSet()])
        while iteration < self.max_iterations and not is_converged:
            iteration += 1

            # u = D(q)
            self.projector.SetSolution(self.decoder.Decode(self.q))

            # Build and project the system evaluated at u = D(q)
            self.scheme.InitializeNonLinIteration(self.strategy_data)
            self.projector.BuildEffectiveSystem()
            phi = self.decoder.Jacobian(self.q)
            self.projector.Project(phi)

            # Solve and update the reduced coordinates
            dq = self.projector.SolveReduced()
            self.q = self.q + dq
            self.scheme.FinalizeNonLinIteration(self.strategy_data)

            # Check convergence
            ratio, absolute_norm = self._DisplacementCriterionNorms((phi @ dq)[free_dofs], self.decoder.Decode(self.q)[free_dofs])
            is_converged = ratio <= self.relative_tolerance or absolute_norm < self.absolute_tolerance
            if self.echo_level > 0:
                KratosMultiphysics.Logger.PrintInfo("FutureRomSolver", "Iteration {}: ratio = {}, absolute norm = {}".format(iteration, ratio, absolute_norm))

        # Leave the database at the last reduced coordinates
        self.projector.SetSolution(self.decoder.Decode(self.q))

        return is_converged

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
