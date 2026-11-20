"""PDE residual evaluation of (ML-predicted) fields through the real solver.

Writes nothing itself: the caller puts a candidate field (e.g. an ML
prediction) into the solution-step database of the unknown variable, then
asks this module for the assembled residual of the discretized equations at
that state. The residual is assembled by the exact same elements, conditions
and builder-and-solver machinery a solve would use, so it is the physics'
own verdict on the field - no surrogate of the residual, the residual.

Uses: physics-informed *monitoring* during training (via TrainModel epoch
callbacks), residual-based active-learning query strategies, and validation.

**This evaluator is NOT differentiable with respect to the field**:
assembly happens in C++ outside any autodiff graph, so it ranks, scores
and logs. For the gradient-carrying counterpart - the SAME assembled
residual wrapped in a torch.autograd.Function whose backward is the
consistent tangent's transpose - see differentiable_residual (built on
this class), which lifts the former scores-only restriction.

Pure Kratos + numpy: this module never imports torch or physicsnemo.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics import python_linear_solver_factory


class ResidualEvaluator:
    """Assembles the free-DOF residual vector of a model part's current state.

    The block builder-and-solver zeroes the rows of fixed DOFs, so every norm
    and per-node value reported here is over free DOFs only - exactly the
    equations a solve would drive to zero.
    """

    def __init__(self, model_part: Kratos.ModelPart, linear_solver=None, scheme=None) -> None:
        self.model_part = model_part
        if linear_solver is None:
            # BuildRHS never solves; the builder constructor merely requires one.
            linear_solver = python_linear_solver_factory.ConstructSolver(
                Kratos.Parameters("""{"solver_type": "skyline_lu_factorization"}"""))
        self._space = Kratos.UblasSparseSpace()
        self.is_static_scheme = scheme is None
        self._scheme = scheme if scheme is not None else Kratos.ResidualBasedIncrementalUpdateStaticScheme()
        self._builder_and_solver = Kratos.ResidualBasedBlockBuilderAndSolver(linear_solver)
        self._initialized = False

    def _Initialize(self) -> None:
        if self._initialized:
            return
        self._builder_and_solver.SetUpDofSet(self._scheme, self.model_part)
        self._builder_and_solver.SetUpSystem(self.model_part)
        if not self.is_static_scheme:
            # Time-integration schemes need their element/condition state set
            # up; the static scheme's versions are no-ops, so this is only
            # done for an explicitly supplied scheme.
            self._scheme.Initialize(self.model_part)
            self._scheme.InitializeElements(self.model_part)
            self._scheme.InitializeConditions(self.model_part)
        self._initialized = True

    def ComputeResidualVector(self) -> Kratos.Vector:
        """Assembles and returns the RHS (residual) vector b of the current
        model part state (free-DOF rows only; fixed-DOF rows are zero)."""
        self._Initialize()
        b = Kratos.Vector(self._builder_and_solver.GetEquationSystemSize())
        self._space.SetToZeroVector(b)
        self._builder_and_solver.BuildRHS(self._scheme, self.model_part, b)
        return b

    def ComputeResidualNorm(self) -> float:
        """Euclidean norm of the assembled residual vector."""
        return float(self._space.TwoNorm(self.ComputeResidualVector()))

    def ComputeNodalResiduals(self) -> dict:
        """Absolute residual per DOF, keyed by (node_id, variable_name).

        Fixed DOFs report 0.0 (their rows are zeroed by the block builder).
        """
        residual = numpy.array(self.ComputeResidualVector(), copy=False)
        values = {}
        for dof in self._builder_and_solver.GetDofSet():
            values[(dof.Id(), dof.GetVariable().Name())] = float(abs(residual[dof.EquationId]))
        return values


def BuildResidualEvaluator(model_part: Kratos.ModelPart, linear_solver=None,
                           scheme=None) -> ResidualEvaluator:
    """Creates a ResidualEvaluator for a model part whose elements already
    carry their DOFs (i.e. after the solver's Initialize, or after a solve).

    Args:
        model_part: The computing model part (the one the solver assembles).
        linear_solver: Optional Kratos linear solver; only needed because the
            builder-and-solver constructor requires one - it is never used.
        scheme: Optional Kratos scheme. The default (None) uses the static
            scheme, which is also the right choice for solvers whose
            ELEMENTS integrate in time (ConvectionDiffusion's transient
            solver does exactly that). Pass a displacement time-integration
            scheme - e.g. ResidualBasedBossakDisplacementScheme - to assemble
            the dynamic residual and effective tangent; the caller must then
            drive InitializeSolutionStep per time step (see
            differentiable_residual.TangentAssembler).
    """
    return ResidualEvaluator(model_part, linear_solver, scheme)
