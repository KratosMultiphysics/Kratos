"""Differentiable physics-informed residuals via physicsnemo.sym.

Promotes the deliberately non-differentiable solver_residuals score to a
real training term: physicsnemo.sym.eq.phy_informer.PhysicsInformer
evaluates SymPy-defined PDE residuals differentiably (grad_method
"autodiff", "least_squares", "finite_difference", "spectral", ...), so the
residual can enter a loss and backpropagate into the surrogate.
physicsnemo.sym ships BUNDLED inside nvidia-physicsnemo 2.1.x - no extra
install. Note the solver-assembled ResidualEvaluator remains the
ground-truth *score*; this module's residuals are the analytic (strong
form) counterparts.

Tensor layout (verified against the installed 2.1.1): everything is
channels-first - fields (B, C, N) point clouds or (B, C, *spatial) grids,
coordinates (B, 3, N) with requires_grad for "autodiff".

Built-in PDEs (equations must be defined inline - physicsnemo.sym no
longer ships pre-built PDE classes):
- "builtin:diffusion": -D lap(u) (+ optional source); D=None makes the
  coefficient a spatial input named "D" - the inverse-problem mechanism
  (feed a trainable tensor, detach "u" via detach_names).
- "builtin:convection_diffusion": c . grad(u) - D lap(u).
- "builtin:linear_elasticity": Navier-Cauchy displacement form,
  -((lmbda+mu) grad(div u) + mu lap(u)) - f, components u_x/u_y/u_z
  (StructuralMechanicsApplication small-strain statics; pass E/nu or the
  Lame parameters directly).
- "builtin:incompressible_navier_stokes": steady convective-form momentum
  rho (v.grad)v_i - mu lap(v_i) + dp/dxi - rho f_i plus continuity div v,
  components velocity_x/_y/_z and pressure (the strong form of
  FluidDynamicsApplication's incompressible solvers).
Anything else: a dotted path "my_module.MyPde" to a PDE subclass,
constructed with the "pde_arguments" block.

Vector fields: physicsnemo.sym derives the informer's required input names
from the sympy Function names, so a width-3 Kratos field must be fed as
three width-1 inputs (u_x/u_y/u_z). The "fields" specs support this via an
optional "components" list; widths 2/3 auto-generate _x/_y/_z suffixes.

torch/physicsnemo are imported lazily; module import stays ML-free.
"""

import numpy

import KratosMultiphysics as Kratos

_GRAD_METHODS = ("autodiff", "least_squares", "finite_difference",
                 "spectral", "meshless_finite_difference")


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.physics_informed requires torch, which could not be "
            "imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportPhysicsNemoSym():
    try:
        from physicsnemo.sym.eq.pde import PDE
        from physicsnemo.sym.eq.phy_informer import PhysicsInformer
        return PDE, PhysicsInformer
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.physics_informed requires physicsnemo.sym, which could "
            "not be imported. It ships bundled with e.g. 'pip install nvidia-physicsnemo' "
            "(>= 2.1).") from e


def MakeDiffusionPde(D=0.1, source=0.0, dim=3):
    """-D lap(u) - source. D=None turns the coefficient into a spatial
    input named "D" (inverse problems).

    dim=2 drops the z-Laplacian and makes u a function of (x, y) only.
    That matters for planar collocation clouds (all nodes at z = 0): with
    the 3D operator, u_zz is completely unconstrained there, and a PINN
    happily dumps curvature into z to cancel the source - the loss
    converges while the in-plane solution has the wrong amplitude.
    """
    PDE, _ = _TryImportPhysicsNemoSym()
    from sympy import Function, Number, Symbol

    if dim not in (2, 3):
        raise ValueError(f"dim must be 2 or 3, got {dim}.")

    class Diffusion(PDE):
        def __init__(self):
            self.dim = dim
            x, y, z = Symbol("x"), Symbol("y"), Symbol("z")
            variables = (x, y) if dim == 2 else (x, y, z)
            u = Function("u")(*variables)
            coefficient = (Function("D")(*variables) if D is None else Number(D))
            laplacian = sum(u.diff(variable, 2) for variable in variables)
            self.equations = {
                "diffusion": -coefficient * laplacian - Number(source),
            }

    return Diffusion()


def MakeConvectionDiffusionPde(D=0.1, cx=0.0, cy=0.0, cz=0.0, source=0.0):
    """c . grad(u) - D lap(u) - source (the stationary
    ConvectionDiffusionApplication residual)."""
    PDE, _ = _TryImportPhysicsNemoSym()
    from sympy import Function, Number, Symbol

    class ConvectionDiffusion(PDE):
        def __init__(self):
            self.dim = 3
            x, y, z = Symbol("x"), Symbol("y"), Symbol("z")
            u = Function("u")(x, y, z)
            self.equations = {
                "convection_diffusion":
                    Number(cx) * u.diff(x) + Number(cy) * u.diff(y) + Number(cz) * u.diff(z)
                    - Number(D) * (u.diff(x, 2) + u.diff(y, 2) + u.diff(z, 2))
                    - Number(source),
            }

    return ConvectionDiffusion()


def MakeLinearElasticityPde(E=None, nu=None, lmbda=1.0, mu=1.0,
                            fx=0.0, fy=0.0, fz=0.0, dim=3):
    """Navier-Cauchy small-strain elastostatics in displacement form:

        elasticity_i = -((lmbda + mu) * d(div u)/dx_i + mu * lap(u_i)) - f_i

    Components are named u_x/u_y/u_z (feed them via a fields spec
    {"name": "u", "width": 3} - the components auto-generate). Pass either
    the Lame parameters (lmbda, mu) or Young's modulus / Poisson's ratio
    (E, nu), which convert via lmbda = E nu / ((1+nu)(1-2nu)),
    mu = E / (2(1+nu)).
    """
    PDE, _ = _TryImportPhysicsNemoSym()
    from sympy import Function, Number, Symbol

    if E is not None or nu is not None:
        if E is None or nu is None:
            raise ValueError("Pass BOTH E and nu (or neither, using lmbda/mu).")
        lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        mu = E / (2.0 * (1.0 + nu))
    if dim not in (2, 3):
        raise ValueError(f"dim must be 2 or 3 [ dim = {dim} ].")

    class LinearElasticity(PDE):
        def __init__(self):
            self.dim = dim
            axes = (Symbol("x"), Symbol("y"), Symbol("z"))[:dim]
            components = [Function(f"u_{name}")(*axes) for name in "xyz"[:dim]]
            body_force = (Number(fx), Number(fy), Number(fz))[:dim]
            divergence = sum(u_j.diff(axis_j) for u_j, axis_j in zip(components, axes))
            self.equations = {}
            for name, u_i, axis_i, f_i in zip("xyz", components, axes, body_force):
                laplacian = sum(u_i.diff(axis_j, 2) for axis_j in axes)
                self.equations[f"elasticity_{name}"] = (
                    -((Number(lmbda) + Number(mu)) * divergence.diff(axis_i)
                      + Number(mu) * laplacian) - f_i)

    return LinearElasticity()


def MakeIncompressibleNavierStokesPde(rho=1.0, mu=1.0, fx=0.0, fy=0.0, fz=0.0, dim=3):
    """Steady incompressible Navier-Stokes, convective form:

        momentum_i = rho (v . grad) v_i - mu lap(v_i) + dp/dx_i - rho f_i
        continuity = div v

    Components are named velocity_x/_y/_z and pressure (feed them via
    fields specs {"name": "velocity", "width": 3} and
    {"name": "pressure", "width": 1}). The strong form of
    FluidDynamicsApplication's incompressible solvers.
    """
    PDE, _ = _TryImportPhysicsNemoSym()
    from sympy import Function, Number, Symbol

    if dim not in (2, 3):
        raise ValueError(f"dim must be 2 or 3 [ dim = {dim} ].")

    class IncompressibleNavierStokes(PDE):
        def __init__(self):
            self.dim = dim
            axes = (Symbol("x"), Symbol("y"), Symbol("z"))[:dim]
            velocity = [Function(f"velocity_{name}")(*axes) for name in "xyz"[:dim]]
            pressure = Function("pressure")(*axes)
            body_force = (Number(fx), Number(fy), Number(fz))[:dim]
            self.equations = {}
            for name, v_i, axis_i, f_i in zip("xyz", velocity, axes, body_force):
                convection = sum(v_j * v_i.diff(axis_j) for v_j, axis_j in zip(velocity, axes))
                laplacian = sum(v_i.diff(axis_j, 2) for axis_j in axes)
                self.equations[f"momentum_{name}"] = (
                    Number(rho) * convection - Number(mu) * laplacian
                    + pressure.diff(axis_i) - Number(rho) * f_i)
            self.equations["continuity"] = sum(
                v_j.diff(axis_j) for v_j, axis_j in zip(velocity, axes))

    return IncompressibleNavierStokes()


def _KwargsFromParameters(settings: Kratos.Parameters) -> dict:
    kwargs = {}
    for key in settings.keys():
        value = settings[key]
        if value.IsBool():
            kwargs[key] = value.GetBool()
        elif value.IsInt():
            kwargs[key] = value.GetInt()
        elif value.IsDouble():
            kwargs[key] = value.GetDouble()
        elif value.IsString():
            kwargs[key] = value.GetString()
        elif value.IsNull():
            kwargs[key] = None
        else:
            raise ValueError(f"Unsupported pde_argument type for \"{key}\".")
    return kwargs


def CreatePde(pde_name: str, pde_arguments: Kratos.Parameters):
    """Resolves "builtin:<name>" or a dotted "module.Class" path."""
    kwargs = _KwargsFromParameters(pde_arguments)
    if pde_name == "builtin:diffusion":
        return MakeDiffusionPde(**kwargs)
    if pde_name == "builtin:convection_diffusion":
        return MakeConvectionDiffusionPde(**kwargs)
    if pde_name == "builtin:linear_elasticity":
        return MakeLinearElasticityPde(**kwargs)
    if pde_name == "builtin:incompressible_navier_stokes":
        return MakeIncompressibleNavierStokesPde(**kwargs)
    if pde_name.startswith("builtin:"):
        raise ValueError(
            f"Unknown builtin PDE \"{pde_name}\". Builtins: \"builtin:diffusion\", "
            "\"builtin:convection_diffusion\", \"builtin:linear_elasticity\", "
            "\"builtin:incompressible_navier_stokes\"; anything else is a dotted "
            "\"module.Class\" path.")
    import importlib
    module_path, _, class_name = pde_name.rpartition(".")
    if not module_path:
        raise ValueError(
            f"\"{pde_name}\" is neither a builtin PDE nor a dotted \"module.Class\" path.")
    return getattr(importlib.import_module(module_path), class_name)(**kwargs)


def _ValidateSettings(settings: Kratos.Parameters) -> Kratos.Parameters:
    defaults = Kratos.Parameters("""{
        "pde"            : "builtin:diffusion",
        "pde_arguments"  : {},
        "residual_names" : [],
        "grad_method"    : "autodiff",
        "fd_dx"          : 0.001,
        "bounds"         : [],
        "detach_names"   : [],
        "weight"         : 1.0,
        "boundary_trim"  : 0,
        "coordinate_channels" : 3,
        "grid_shape"     : []
    }""")
    settings = settings.Clone()
    settings.ValidateAndAssignDefaults(defaults)
    return settings


def CreatePhysicsInformer(settings: Kratos.Parameters):
    """Builds a PhysicsInformer from a "physics_loss"-style settings block.

    Settings (defaults): pde ("builtin:diffusion"), pde_arguments ({}),
    residual_names ([] = every equation of the PDE), grad_method
    ("autodiff"), fd_dx (0.001), bounds ([]), detach_names ([]).

    Returns:
        (informer, residual_names, grad_method)
    """
    _, PhysicsInformer = _TryImportPhysicsNemoSym()
    settings = _ValidateSettings(settings)

    grad_method = settings["grad_method"].GetString()
    if grad_method not in _GRAD_METHODS:
        raise ValueError(
            f"Unsupported grad_method \"{grad_method}\". Use one of {_GRAD_METHODS}.")

    pde = CreatePde(settings["pde"].GetString(), settings["pde_arguments"])
    residual_names = settings["residual_names"].GetStringArray()
    if not residual_names:
        residual_names = list(pde.equations.keys())

    kwargs = {}
    if grad_method in ("finite_difference", "meshless_finite_difference"):
        kwargs["fd_dx"] = settings["fd_dx"].GetDouble()
    if grad_method == "spectral":
        bounds = list(settings["bounds"].GetVector())
        if bounds:
            kwargs["bounds"] = bounds
    detach_names = settings["detach_names"].GetStringArray()
    if detach_names:
        kwargs["detach_names"] = detach_names

    informer = PhysicsInformer(
        required_outputs=residual_names, equations=pde, grad_method=grad_method, **kwargs)
    return informer, residual_names, grad_method


def _SplitFields(prediction, field_specs):
    """(N, total_width) prediction -> {component_name: (1, 1, N)} channels-first.

    Every component becomes its OWN width-1 informer input, keyed by its
    component name - physicsnemo.sym resolves inputs by the sympy Function
    names, so a width-3 field named "u" must arrive as u_x/u_y/u_z. Width-1
    fields keep their single component named after the field itself.
    """
    fields = {}
    offset = 0
    for name, width, components in field_specs:
        for column, component in enumerate(components):
            fields[component] = prediction[:, offset + column:offset + column + 1].T[None]
        offset += width
    if offset != prediction.shape[1]:
        raise ValueError(
            f"The physics fields consume {offset} channels but the prediction has "
            f"{prediction.shape[1]}.")
    return fields


_COMPONENT_SUFFIXES = ("x", "y", "z")


def _ReadFieldSpecs(settings: Kratos.Parameters):
    """Parses the "fields" block into (name, width, components) triples.

    "components" names the per-channel sympy functions; it defaults to
    [name] for width 1 and to name_x/name_y/name_z for widths 2/3 (matching
    the builtin vector PDEs). Width > 3 requires an explicit list. (Before
    the components extension, a width > 1 field was passed as one
    multi-channel input under the field name, which could never match any
    informer required_inputs - the auto-split is a fix, not a break.)
    """
    if settings["fields"].size() == 0:
        return [("u", 1, ["u"])]
    specs = []
    for i in range(settings["fields"].size()):
        entry = settings["fields"][i]
        components = None
        if entry.Has("components"):
            components = entry["components"].GetStringArray()
            entry.RemoveValue("components")
        entry.ValidateAndAssignDefaults(Kratos.Parameters('{"name": "u", "width": 1}'))
        name = entry["name"].GetString()
        width = entry["width"].GetInt()
        if components is None:
            if width == 1:
                components = [name]
            elif width <= 3:
                components = [f"{name}_{suffix}" for suffix in _COMPONENT_SUFFIXES[:width]]
            else:
                raise ValueError(
                    f"Field \"{name}\" has width {width} > 3; give an explicit "
                    "\"components\" list naming its per-channel sympy functions.")
        if len(components) != width:
            raise ValueError(
                f"Field \"{name}\" has width {width} but {len(components)} components "
                f"{components}.")
        specs.append((name, width, components))
    return specs


def MakePhysicsLossTerm(settings: Kratos.Parameters, connectivity_provider=None):
    """Builds an extra loss term for training_utils.TrainModel.

    The returned callable term(model, inputs, prediction) -> scalar adds
    weight * mean(residual^2) over the configured residuals. Per
    grad_method:

    - "autodiff": the first "coordinate_channels" input channels are the
      point coordinates (the point-cloud "generic" training layout); the
      model is RE-RUN with those channels requiring grad so du/dx exists,
      and the residual is evaluated at the points.
    - "least_squares": connectivity_provider() -> (nodes (N, 3),
      edges (2, E)) (e.g. graph_bridge/particle_bridge output); the
      residual differentiates the prediction itself - no re-forward.
    - "finite_difference"/"spectral": the prediction rows are a flattened
      grid; "grid_shape" (e.g. [16, 16, 16]) reshapes them.

    Settings: CreatePhysicsInformer's keys plus weight (1.0),
    coordinate_channels (3), fields ([{"name": "u", "width": 1,
    "components": [...]}] - the PDE unknowns, matched to the prediction
    channels in order; see _ReadFieldSpecs for the component naming),
    grid_shape, and boundary_trim (0): for the grid methods, drop that
    many cells per side of the RESIDUAL grid before averaging. The
    upstream finite-difference stencils are wrong on the outermost shell
    of a non-periodic field (a field with an FD-exact interior residual of
    zero still averages O(1) over the full grid), so a non-zero trim is
    what makes the term usable as a physics loss on bounded domains.
    """
    torch = _TryImportTorch()

    settings = settings.Clone()
    if not settings.Has("fields"):
        settings.AddEmptyArray("fields")
    field_specs = _ReadFieldSpecs(settings)
    settings.RemoveValue("fields")

    settings = _ValidateSettings(settings)
    informer, residual_names, grad_method = CreatePhysicsInformer(settings)
    weight = settings["weight"].GetDouble()
    coordinate_channels = settings["coordinate_channels"].GetInt()
    grid_shape = tuple(int(n) for n in settings["grid_shape"].GetVector())
    boundary_trim = settings["boundary_trim"].GetInt()
    if boundary_trim < 0:
        raise ValueError(f"\"boundary_trim\" must be >= 0, got {boundary_trim}.")
    if boundary_trim and min(grid_shape or (0,)) <= 2 * boundary_trim:
        raise ValueError(
            f"\"boundary_trim\" = {boundary_trim} leaves no interior for "
            f"grid_shape {list(grid_shape)}.")

    def term(model, inputs, prediction):
        if grad_method == "autodiff":
            # the (1, 3, N) coordinates tensor the informer differentiates
            # against must be the leaf the model input derives from, or the
            # autodiff graph breaks at the transpose
            coordinates = inputs[:, :coordinate_channels].detach().clone().T[None]
            coordinates.requires_grad_(True)
            rebuilt = torch.cat([coordinates[0].T, inputs[:, coordinate_channels:]], dim=1)
            prediction = model(rebuilt)
            informer_inputs = _SplitFields(prediction, field_specs)
            informer_inputs["coordinates"] = coordinates
        elif grad_method == "least_squares":
            # point layout (N, C): fields (N, width), coordinates/nodes (N, 3),
            # edges (E, 2)
            if connectivity_provider is None:
                raise ValueError(
                    "grad_method \"least_squares\" needs a connectivity_provider returning "
                    "(nodes (N, 3), edges (2, E)).")
            coordinates, edges = connectivity_provider()
            coordinates = torch.as_tensor(numpy.asarray(coordinates), dtype=prediction.dtype)
            edges = torch.as_tensor(numpy.asarray(edges)).long()
            if edges.shape[0] == 2 and edges.shape[1] != 2:
                edges = edges.T  # graph_bridge's (2, E) -> upstream's (E, 2)
            informer_inputs = {}
            offset = 0
            for name, width, components in field_specs:
                for column, component in enumerate(components):
                    informer_inputs[component] = prediction[
                        :, offset + column:offset + column + 1]  # (N, 1) per component
                offset += width
            informer_inputs["coordinates"] = coordinates
            # upstream's "nodes" are node IDS of shape (N, 1)
            informer_inputs["nodes"] = torch.arange(coordinates.shape[0]).reshape(-1, 1)
            informer_inputs["edges"] = edges
        elif grad_method in ("finite_difference", "spectral"):
            if not grid_shape:
                raise ValueError(
                    f"grad_method \"{grad_method}\" needs \"grid_shape\" to reshape the "
                    "flattened prediction rows into a grid.")
            informer_inputs = {}
            offset = 0
            for name, width, components in field_specs:
                for column, component in enumerate(components):
                    chunk = prediction[:, offset + column:offset + column + 1]
                    informer_inputs[component] = chunk.T.reshape((1, 1) + grid_shape)
                offset += width
        else:
            raise ValueError(
                f"grad_method \"{grad_method}\" is not supported as a training loss term.")

        residuals = informer.forward(informer_inputs)
        loss = None
        for name in residual_names:
            residual = residuals[name]
            if boundary_trim and grad_method in ("finite_difference", "spectral"):
                trim = (slice(None),) * (residual.ndim - len(grid_shape)) + tuple(
                    slice(boundary_trim, -boundary_trim) for _ in grid_shape)
                residual = residual[trim]
            contribution = residual.square().mean()
            loss = contribution if loss is None else loss + contribution
        return weight * loss

    return term
