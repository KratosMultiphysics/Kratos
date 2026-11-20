"""The exact discrete FEM residual as a diffusion-guidance operator.

Diffusion posterior sampling (DPS) steers an already trained denoiser at
SAMPLING time: at every step the sampler's current estimate of the clean
field is scored by a differentiable operator, and the gradient of that
score is added to the score of the diffusion model. physicsnemo ships the
guidance terms; what a solver-coupled surrogate wants as the operator is
the physics itself.

This module supplies that operator: a callable A(x0) mapping a batch of
generated GRIDS to the assembled residual b(u) of the real Kratos system,
built on physics.differentiable_residual's autograd Function (forward =
BuildRHS, backward = the consistent tangent's transpose). Sampling with

    guidance   = {"type" : "model_consistency"}
    observation = 0

then pushes the ensemble toward fields whose PDE residual vanishes -
without retraining anything, and with the model's own prior still doing
the work everywhere the residual does not constrain.

The chain A runs, in order, is exactly the chain the inference process
runs on the way out, so a guided sample is scored in the same physical
units it will eventually be written in:

    normalized grid -> de-normalize (the model card's
    "output_normalization") -> un-squeeze the thin axis -> trilinear
    interpolation at the nodes (the differentiable
    grid_bridge.InterpolateGridAtPointsTorch) -> gather DOFs -> overwrite
    the fixed DOFs with the model part's own Dirichlet values -> assemble.

torch is imported lazily.
"""

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.bridges import grid_bridge
from KratosMultiphysics.PhysicsNeMoApplication.physics import differentiable_residual


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.diffusion_residual_operator requires torch, which "
            "could not be imported. Install it with e.g. 'pip install torch'.") from e


class KratosResidualObservationOperator:
    """A(x0) = the assembled residual of the generated field.

    Callable with a (B, C, *spatial) batch of normalized grids, returning a
    (B, n_equations) float64 tensor that carries gradients back to the
    grid. The matching observation is zero - a field that solves the PDE
    has no residual - and is exposed as .observation.

    IMPORTANT: assembling writes the trial field into the model part's own
    solution-step database (that is how Kratos evaluates a residual), so
    the operator saves the incoming state on construction and .Restore()
    puts it back. A sampler run must call Restore() when it is done, or the
    model part is left holding the last trial field instead of the solve's.
    """

    def __init__(self, model_part: Kratos.ModelPart, field_specs, grid_shape,
                 bounding_box, squeeze_axis=None, output_normalization=None,
                 linear_solver=None, scheme=None,
                 use_stored_fixed_values: bool = True) -> None:
        torch = _TryImportTorch()
        self._torch = torch
        self.model_part = model_part
        self.field_specs = list(field_specs)
        self.grid_shape = tuple(int(n) for n in grid_shape)
        self.bounding_box = bounding_box
        self.squeeze_axis = squeeze_axis
        self.output_normalization = output_normalization
        self.use_stored_fixed_values = use_stored_fixed_values

        self._assembler = differentiable_residual.TangentAssembler(
            model_part, linear_solver=linear_solver, scheme=scheme)
        self._dof_map = differentiable_residual.DofFieldMap(self._assembler, self.field_specs)
        self._gather_index = self._dof_map.TorchGatherIndex()
        self._fixed_mask = torch.from_numpy(self._dof_map.fixed_mask)
        # the solve's own state: both the Dirichlet values to re-impose and
        # what Restore() puts back
        self._saved_state = self._dof_map.ReadDofVector()
        self._fixed_values = torch.from_numpy(self._saved_state.copy())
        self._node_coordinates = numpy.array(
            [[node.X, node.Y, node.Z] for node in model_part.Nodes], dtype=float)

    @property
    def observation(self):
        """The residual a field satisfying the PDE would have: zero."""
        return self._torch.zeros(
            (1, self._dof_map.n_equations), dtype=self._torch.float64)

    def Restore(self) -> None:
        """Puts the model part's DOFs back as they were before sampling."""
        self._dof_map.WriteDofVector(self._saved_state)

    def _Unsqueeze(self, grid):
        """(C, H, W) -> (C, D, H, W), duplicating across the thin axis."""
        if self.squeeze_axis is None:
            return grid
        thin_size = self.grid_shape[self.squeeze_axis]
        axis = 1 + self.squeeze_axis
        return grid.unsqueeze(axis).expand(
            *grid.shape[:axis], thin_size, *grid.shape[axis:])

    def __call__(self, x0):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry

        torch = self._torch
        if x0.ndim < 3:
            raise ValueError(
                "The residual operator consumes a (B, C, *spatial) batch of grids, got "
                f"shape {tuple(x0.shape)}.")
        residuals = []
        for index in range(x0.shape[0]):
            grid = model_registry.ApplyOutputNormalization(
                x0[index], self.output_normalization, channel_axis=0)
            grid = self._Unsqueeze(grid).to(torch.float64)
            values = grid_bridge.InterpolateGridAtPointsTorch(
                grid, self.bounding_box, self._node_coordinates)  # (N, total_width)
            if values.shape[1] != self._dof_map.total_width:
                raise ValueError(
                    f"The generated grid has {values.shape[1]} channels but the DOF "
                    f"fields need {self._dof_map.total_width}; the guidance field specs "
                    "must match the model's output channels.")
            u = values.reshape(-1)[self._gather_index]
            if self.use_stored_fixed_values:
                # the Dirichlet rows are the solve's, not the model's guess
                u = torch.where(self._fixed_mask, self._fixed_values, u)
            residuals.append(differentiable_residual.KratosResidualFunction.Apply(
                u, self._assembler, self._dof_map))
        return torch.stack(residuals)


def MakeKratosResidualObservationOperator(model_part: Kratos.ModelPart,
                                          settings: Kratos.Parameters,
                                          output_normalization=None,
                                          linear_solver=None, scheme=None):
    """Builds the residual operator from a "guidance" block's settings.

    Args:
        model_part: The model part whose system the residual is assembled
            from. Its DOFs must be the fields the model generates.
        settings: Kratos Parameters; defaults:
            residual_fields ([] - the (variable_name, data_location) specs
            of the generated channels, in channel order; node_historical
            only, since DOFs are historical), grid_shape ([8, 8, 2]),
            bounding_box ([]), squeeze_axis (-1),
            use_stored_fixed_values (true).
        output_normalization: The model card's "output_normalization"
            entry, so the residual sees physical values.
        linear_solver, scheme: Forwarded to the TangentAssembler (a
            transient case needs its scheme).

    Returns:
        A KratosResidualObservationOperator.
    """
    default_settings = Kratos.Parameters("""{
        "residual_fields"          : [],
        "grid_shape"               : [8, 8, 2],
        "bounding_box"             : [],
        "squeeze_axis"             : -1,
        "use_stored_fixed_values"  : true
    }""")
    settings.ValidateAndAssignDefaults(default_settings)

    field_specs = [
        (settings["residual_fields"][i]["variable_name"].GetString(),
         settings["residual_fields"][i]["data_location"].GetString()
         if settings["residual_fields"][i].Has("data_location") else "node_historical")
        for i in range(settings["residual_fields"].size())
    ]
    if not field_specs:
        raise ValueError(
            "\"residual_fields\" is empty: the residual operator needs the DOF fields "
            "the generated channels correspond to.")

    grid_shape = tuple(int(n) for n in settings["grid_shape"].GetVector())
    box = settings["bounding_box"].GetVector()
    if len(box) == 6:
        bounding_box = (numpy.array(box[:3]), numpy.array(box[3:]))
    elif len(box) == 0:
        bounding_box = grid_bridge.ComputeBoundingBox(model_part)
    else:
        raise ValueError("\"bounding_box\" must be empty or [x0,y0,z0,x1,y1,z1].")

    squeeze_axis = settings["squeeze_axis"].GetInt()
    return KratosResidualObservationOperator(
        model_part, field_specs, grid_shape, bounding_box,
        squeeze_axis=None if squeeze_axis == -1 else squeeze_axis,
        output_normalization=output_normalization,
        linear_solver=linear_solver, scheme=scheme,
        use_stored_fixed_values=settings["use_stored_fixed_values"].GetBool())
