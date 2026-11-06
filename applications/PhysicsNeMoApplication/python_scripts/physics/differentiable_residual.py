"""The real assembled FEM residual as a differentiable torch operation.

This module lifts the long-standing restriction that solver-assembled
residuals are scores, never losses. solver_residuals.ResidualEvaluator
remains the plain (non-differentiable) evaluator; here a
torch.autograd.Function wraps the SAME assembly so the residual carries
gradients:

    forward(u)  : write u into the DOFs, BuildRHS -> b(u)
    backward(g) : assemble the consistent tangent K at u,
                  return (d b/d u)^T g = -(masked K)^T g

The vector-Jacobian product is a single TRANSPOSE MATVEC - no linear solve
(the adjoint solve A^T lambda = dJ/du belongs to sensitivity_utils). The
sign and masking follow the block builder-and-solver's conventions, pinned
by tests rather than trusted on paper:

- Kratos statics assembles b = f - K u for linear problems, and Build's LHS
  is the consistent tangent for nonlinear ones, so d b/d u = -K on the free
  rows.
- BuildRHS ZEROES the fixed-DOF rows (Build does not), so backward masks
  the incoming cotangent with the fixed-DOF mask before the matvec.

The assembly path is the documented core Python pattern (PFEM2's
strategy_python, StructuralMechanics' scipy base solver):
space.CreateEmptyMatrixPointer() -> builder.ResizeAndInitializeVectors ->
SetToZero -> builder.Build -> KratosMultiphysics.scipy_conversion_tools
.to_csr (which reads the CompressedMatrix CSR triple
value_data()/index2_data()/index1_data() - the value view is zero-copy).

Scope: statics by default, and one step of a TRANSIENT problem at fixed
step history. Two transient flavours, both giving d b/d u = -(masked
effective tangent) at the current step:

- Element-integrated time stepping (ConvectionDiffusion's transient
  solver): nothing changes here - the solver itself uses the static
  scheme, and the time dependence enters through ProcessInfo (DELTA_TIME,
  THETA) and the solution-step buffer.
- Displacement schemes (Bossak/Newmark/BDF): construct the assembler with
  scheme=..., call InitializeSolutionStep() once per time step (that is
  where the schemes compute their coefficients), and the assembler
  refreshes the derived VELOCITY/ACCELERATION from the written DOFs before
  every assembly. builder.Build then yields K_eff = K + M(1-alpha)c0 + D c1
  and BuildRHS the dynamic residual, which is why the same wrapper works.

Both forward and backward WRITE the model part's solution-step database -
a documented side effect (backward re-writes the forward state first, so
interleaved calls stay consistent).

torch is imported lazily; scipy imports stay function-local (it is a hard
Kratos-core dependency via scipy_conversion_tools, but this application
keeps it out of module scope like every optional import).
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.physics import solver_residuals
def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.differentiable_residual requires torch, which could "
            "not be imported. Install it with e.g. 'pip install torch'.") from e


class TangentAssembler(solver_residuals.ResidualEvaluator):
    """ResidualEvaluator plus consistent-tangent assembly to scipy CSR.

    The system pointers are created once (ResizeAndInitializeVectors) and
    reused: subsequent calls only zero and re-assemble into the existing
    sparsity graph.
    """

    def __init__(self, model_part: Kratos.ModelPart, linear_solver=None, scheme=None) -> None:
        super().__init__(model_part, linear_solver, scheme)
        self._pA = None
        self._pDx = None
        self._pb = None

    def _EnsureSystemPointers(self) -> None:
        self._Initialize()
        if self._pA is None:
            self._pA = self._space.CreateEmptyMatrixPointer()
            self._pDx = self._space.CreateEmptyVectorPointer()
            self._pb = self._space.CreateEmptyVectorPointer()
            self._builder_and_solver.ResizeAndInitializeVectors(
                self._scheme, self._pA, self._pDx, self._pb, self.model_part)

    def InitializeSolutionStep(self) -> None:
        """Prepares the assembler for the model part's CURRENT time step.

        Call once per time step, after the solver has advanced time and
        cloned the buffer, before assembling. Time-integration schemes
        compute their coefficients here (Bossak/Newmark's c0..c5 come from
        DELTA_TIME in InitializeSolutionStep), so assembling a dynamic
        system without this call uses stale coefficients. Harmless (and
        cheap) for the default static scheme.
        """
        self._EnsureSystemPointers()
        self._builder_and_solver.InitializeSolutionStep(
            self.model_part, self._pA, self._pDx, self._pb)
        self._scheme.InitializeSolutionStep(
            self.model_part, self._pA, self._pDx, self._pb)

    def RefreshTimeDerivatives(self) -> None:
        """Recomputes the scheme's derived fields from the current DOFs.

        Displacement time-integration schemes express velocity and
        acceleration as affine functions of the current displacement at
        fixed step history; after writing a new u into the DOFs those
        derived fields are stale, and the assembled dynamic residual would
        mix states. scheme.Update with a zero increment recomputes them.
        No-op for the static scheme (nothing is derived).
        """
        if self.is_static_scheme:
            return
        self._EnsureSystemPointers()
        self._space.SetToZeroVector(self._pDx)
        self._scheme.Update(self.model_part, self._builder_and_solver.GetDofSet(),
                            self._pA, self._pDx, self._pb)

    def ComputeResidualVector(self) -> Kratos.Vector:
        """Assembled residual at the current state (dynamic residual when a
        time-integration scheme was supplied; the derived velocity and
        acceleration are refreshed first)."""
        self.RefreshTimeDerivatives()
        return super().ComputeResidualVector()

    def ComputeSystem(self, apply_dirichlet: bool = False):
        """Assembles the tangent matrix and RHS at the current state.

        Args:
            apply_dirichlet: False (default) returns the raw assembled
                system - the Jacobian consistent with BuildRHS's zeroed
                fixed rows (for vector-Jacobian products). True applies the
                builder's Dirichlet treatment (unit diagonal on fixed
                rows/columns) - the well-posed operator for adjoint solves.

        Returns:
            (K, b): scipy.sparse.csr_matrix and (n_eq,) float64 numpy array.
        """
        self._EnsureSystemPointers()
        self.RefreshTimeDerivatives()
        self._space.SetToZeroMatrix(self._pA)
        self._space.SetToZeroVector(self._pb)
        self._space.SetToZeroVector(self._pDx)
        self._builder_and_solver.Build(self._scheme, self.model_part, self._pA, self._pb)
        if apply_dirichlet:
            self._builder_and_solver.ApplyDirichletConditions(
                self._scheme, self.model_part, self._pA, self._pDx, self._pb)

        import KratosMultiphysics.scipy_conversion_tools as scipy_conversion_tools
        # to_csr copies the index arrays but wraps the zero-copy value view;
        # copy() detaches the result from the reused assembly buffers.
        return (scipy_conversion_tools.to_csr(self._pA).copy(),
                numpy.array(self._pb, copy=True))

    def ComputeTangentMatrix(self, apply_dirichlet: bool = False):
        """The consistent tangent K at the current state (scipy CSR)."""
        return self.ComputeSystem(apply_dirichlet)[0]


class DofFieldMap:
    """Maps between the app's (N, total_width) nodal-field layout and the
    equation-id-ordered DOF vector of a block builder-and-solver.

    Built once from the assembler's DOF set and field specs
    [(variable_name, data_location)] (node_historical only - DOFs live in
    the historical database). Node rows follow model_part.Nodes iteration
    order - the same order every gather in this application uses. Column
    layout concatenates the specs' components in order (scalar variables
    contribute one column; array variables three, X/Y/Z).
    """

    def __init__(self, assembler: solver_residuals.ResidualEvaluator, field_specs) -> None:
        assembler._Initialize()
        self.model_part = assembler.model_part
        self.n_equations = assembler._builder_and_solver.GetEquationSystemSize()

        # column of each DOF variable name in the concatenated field layout
        column_of = {}
        offset = 0
        for variable_name, data_location in field_specs:
            if data_location != "node_historical":
                raise ValueError(
                    f"DOF field \"{variable_name}\" must live in \"node_historical\" "
                    f"(DOFs are historical); got \"{data_location}\".")
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            if isinstance(variable, Kratos.DoubleVariable):
                column_of[variable_name] = offset
                offset += 1
            else:  # array variable: the DOFs are the components
                for axis in "XYZ":
                    column_of[f"{variable_name}_{axis}"] = offset
                    offset += 1
        self.total_width = offset
        # exposed so a caller can place a per-variable quantity (an
        # objective's weights, say) into the right column block
        self.column_of = column_of

        node_row = {node.Id: row for row, node in enumerate(self.model_part.Nodes)}
        self.n_nodes = len(node_row)

        dof_set = assembler._builder_and_solver.GetDofSet()
        self._dof_set = dof_set
        equation_ids = numpy.array(dof_set.GetEquationIds(), dtype=numpy.int64)
        flat_index = numpy.empty(len(equation_ids), dtype=numpy.int64)
        fixed = numpy.empty(len(equation_ids), dtype=bool)
        for i, dof in enumerate(dof_set):
            name = dof.GetVariable().Name()
            if name not in column_of:
                raise ValueError(
                    f"DOF variable \"{name}\" is not covered by the field specs "
                    f"{[spec[0] for spec in field_specs]}.")
            flat_index[i] = node_row[dof.Id()] * self.total_width + column_of[name]
            fixed[i] = dof.IsFixed()

        # dof-set order -> equation order
        self._equation_ids = equation_ids
        self.gather_index = numpy.empty(self.n_equations, dtype=numpy.int64)
        self.gather_index[equation_ids] = flat_index
        self.fixed_mask = numpy.zeros(self.n_equations, dtype=bool)
        self.fixed_mask[equation_ids] = fixed

    def FieldsToDofVector(self, values):
        """(N, total_width) nodal values -> (n_eq,) equation-ordered vector."""
        return numpy.asarray(values, dtype=numpy.float64).reshape(-1)[self.gather_index]

    def DofVectorToFields(self, u):
        """(n_eq,) equation-ordered vector -> (N, total_width) nodal values
        (non-DOF entries stay zero)."""
        out = numpy.zeros(self.n_nodes * self.total_width)
        out[self.gather_index] = numpy.asarray(u, dtype=numpy.float64)
        return out.reshape(self.n_nodes, self.total_width)

    def WriteDofVector(self, u) -> None:
        """Writes an equation-ordered vector into the DOFs' solution-step values."""
        u = numpy.asarray(u, dtype=numpy.float64)
        self._dof_set.SetValues(Kratos.Vector(numpy.ascontiguousarray(u[self._equation_ids])))

    def ReadDofVector(self):
        """Reads the DOFs' solution-step values as an equation-ordered vector."""
        values = numpy.array(self._dof_set.GetValues(), copy=False)
        out = numpy.empty(self.n_equations)
        out[self._equation_ids] = values
        return out

    def TorchGatherIndex(self):
        """gather_index as a LongTensor, for differentiable
        prediction.reshape(-1)[index] gathers."""
        torch = _TryImportTorch()
        return torch.from_numpy(self.gather_index)


class KratosResidualFunction:
    """Namespace holder; the actual autograd.Function is created lazily so
    importing this module never requires torch. Use Apply(...)."""

    _function = None

    @classmethod
    def _Get(cls):
        if cls._function is None:
            torch = _TryImportTorch()

            class _KratosResidualFunction(torch.autograd.Function):
                @staticmethod
                def forward(ctx, u_dofs, assembler, dof_map):
                    u = u_dofs.detach().cpu().to(torch.float64).numpy()
                    dof_map.WriteDofVector(u)  # state <- u
                    b = numpy.array(assembler.ComputeResidualVector(), copy=True)
                    ctx.assembler = assembler
                    ctx.dof_map = dof_map
                    ctx.u = u
                    return torch.from_numpy(b)  # float64; fixed rows exactly 0

                @staticmethod
                def backward(ctx, grad_b):
                    # restore the linearization state before assembling
                    ctx.dof_map.WriteDofVector(ctx.u)
                    K = ctx.assembler.ComputeTangentMatrix(apply_dirichlet=False)
                    g = grad_b.detach().cpu().to(torch.float64).numpy().copy()
                    g[ctx.dof_map.fixed_mask] = 0.0  # BuildRHS zeroed those rows
                    vjp = -(K.T @ g)                 # d b/d u = -(masked K)
                    return torch.from_numpy(vjp), None, None

            cls._function = _KratosResidualFunction
        return cls._function

    @classmethod
    def Apply(cls, u_dofs, assembler, dof_map):
        """b(u): the assembled residual as a differentiable torch tensor.

        Args:
            u_dofs: (n_eq,) float64 torch tensor, equation-id ordered.
            assembler: A TangentAssembler on the computing model part.
            dof_map: The matching DofFieldMap.
        """
        return cls._Get().apply(u_dofs, assembler, dof_map)


def MakeExactResidualLossTerm(settings: Kratos.Parameters, model_part: Kratos.ModelPart,
                              full_inputs_provider, linear_solver=None):
    """Builds a gradient-carrying EXACT-residual loss term for
    training_utils.TrainModel(..., extra_loss_terms=[term]).

    Unlike physics_informed.MakePhysicsLossTerm (analytic strong-form
    residuals), this term evaluates the DISCRETE residual through the real
    element/builder assembly - the physics' own verdict, now with
    gradients. Because TrainModel batches over nodes, the term re-runs the
    model on the full case inputs (the same re-forward idiom as the
    autodiff physics term):

        loss = weight * mean(b(u_pred)^2)

    Settings:
        fields: [{"variable_name": "TEMPERATURE",
                  "data_location": "node_historical"}] - the unknown
            fields, concatenated in the model's output-channel order.
        weight: 1.0.
        use_stored_fixed_values: true - evaluate the residual with the
            model part's stored (Dirichlet) values on fixed DOFs instead
            of the model's predictions there; gradients then flow through
            free DOFs only.

    Args:
        model_part: The computing model part (elements carry their DOFs -
            i.e. after the solver's Initialize / a solve).
        full_inputs_provider: callable() -> (N, C_in) torch tensor of the
            full case inputs the model maps to the unknown fields.
        linear_solver: Optional; the builder constructor merely requires one.

    The assembler and map are built lazily on the first call. Forward and
    backward write the model part's solution-step database (documented
    side effect of the wrapped Function).
    """
    torch = _TryImportTorch()

    defaults = Kratos.Parameters("""{
        "fields" : [
            {
                "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                "data_location" : "node_historical"
            }
        ],
        "weight" : 1.0,
        "use_stored_fixed_values" : true
    }""")
    settings = settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)
    for i in range(settings["fields"].size()):
        settings["fields"][i].ValidateAndAssignDefaults(defaults["fields"][0])
    field_specs = [
        (settings["fields"][i]["variable_name"].GetString(),
         settings["fields"][i]["data_location"].GetString())
        for i in range(settings["fields"].size())
    ]
    weight = settings["weight"].GetDouble()
    use_stored_fixed_values = settings["use_stored_fixed_values"].GetBool()

    state = {}

    def _Setup():
        assembler = TangentAssembler(model_part, linear_solver)
        dof_map = DofFieldMap(assembler, field_specs)
        state["assembler"] = assembler
        state["dof_map"] = dof_map
        state["gather_index"] = dof_map.TorchGatherIndex()
        state["fixed_mask"] = torch.from_numpy(dof_map.fixed_mask)
        state["fixed_values"] = torch.from_numpy(dof_map.ReadDofVector()).to(torch.float64)

    def term(model, inputs, prediction):
        if not state:
            _Setup()
        full_prediction = model(full_inputs_provider()).to(torch.float64)
        u_dofs = full_prediction.reshape(-1)[state["gather_index"]]
        if use_stored_fixed_values:
            u_dofs = torch.where(state["fixed_mask"], state["fixed_values"], u_dofs)
        b = KratosResidualFunction.Apply(u_dofs, state["assembler"], state["dof_map"])
        return weight * b.square().mean()

    return term
