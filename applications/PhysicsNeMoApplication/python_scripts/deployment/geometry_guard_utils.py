"""Out-of-distribution detection on the SHAPE, not the fields.

The shipped OOD guard (ood_guard_utils) looks at a model's INPUTS: the
field values it is about to consume. That catches a solve whose
temperatures drifted out of the training range, but not a surrogate trained
on one family of geometries being handed a different one - the fields can
look perfectly ordinary on a shape the model has never seen.

physicsnemo's experimental geometry guardrail fits a density model to
rotation- and translation-sensitive descriptors of a triangular SURFACE and
reports where a new geometry falls in that distribution. The surface it
wants is exactly what the mesh bridge's BoundarySurface already produces:
outward-oriented triangles covering the model part's boundary, with every
point kept.

Usage mirrors the field guard: fit on the training family once, save the
sidecar, and name it in a deployment process's "geometry_guard" block.

    guard = CreateGeometryGuard(Kratos.Parameters("{}"))
    FitGeometryGuard(guard, [part_a, part_b, ...])
    SaveGeometryGuard(guard, "surrogate.mdlus.geometry_guard.pt")

physicsnemo.experimental has no API stability guarantee. torch and
physicsnemo are imported lazily.
"""

import warnings

import numpy

import KratosMultiphysics as Kratos

_GUARD_POLICIES = ("advisory", "strict", "ignore")
_GUARD_STATUSES = ("OK", "WARN", "REJECT")


def _TryImportTorch():
    try:
        import torch
        return torch
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.geometry_guard_utils requires torch, which could "
            "not be imported. Install it with e.g. 'pip install torch'.") from e


def _TryImportGeometryGuardrail():
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # physicsnemo.experimental warns on import
            from physicsnemo.experimental.guardrails.geometry import (
                GeometryGuardrail, extract_features, validate_mesh)
        return GeometryGuardrail, extract_features, validate_mesh
    except ImportError as e:
        raise ImportError(
            "The geometry guardrail requires physicsnemo >= 2.2 (its experimental "
            "guardrails), which could not be imported. Install it with e.g. "
            "'pip install -U nvidia-physicsnemo'.") from e


def SurfaceOfModelPart(model_part: Kratos.ModelPart, source_container: str = "Elements",
                       tessellation_mode: str = "smallest_id_diagonal"):
    """The model part's boundary as a physicsnemo triangular surface Mesh.

    Volume parts are tessellated and their outward-oriented boundary
    extracted; a surface part passes straight through. This is the same
    BoundarySurface the signed-distance field is computed against, so a
    geometry the guardrail accepts is the geometry the SDF features describe.
    """
    from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import (
        domain_mesh_builder, spatial)

    mesh, _ = domain_mesh_builder.BuildMesh(
        model_part, source_container=source_container,
        tessellation_mode=tessellation_mode)
    return spatial.BoundarySurface(mesh)


def CreateGeometryGuard(settings: Kratos.Parameters):
    """Builds an unfitted GeometryGuardrail.

    Args:
        settings: Kratos Parameters; defaults:
            method ("gmm"), gmm_components (1), warn_pct (99.0),
            reject_pct (99.9), poly_degree (2), random_state (0),
            device ("cpu").

    Returns:
        A physicsnemo GeometryGuardrail.
    """
    default_settings = Kratos.Parameters("""{
        "method"         : "gmm",
        "gmm_components" : 1,
        "warn_pct"       : 99.0,
        "reject_pct"     : 99.9,
        "poly_degree"    : 2,
        "random_state"   : 0,
        "device"         : "cpu"
    }""")
    settings.ValidateAndAssignDefaults(default_settings)
    GeometryGuardrail, _, _ = _TryImportGeometryGuardrail()
    return GeometryGuardrail(
        method=settings["method"].GetString(),
        gmm_components=settings["gmm_components"].GetInt(),
        warn_pct=settings["warn_pct"].GetDouble(),
        reject_pct=settings["reject_pct"].GetDouble(),
        poly_degree=settings["poly_degree"].GetInt(),
        random_state=settings["random_state"].GetInt(),
        device=settings["device"].GetString())


def _AsMeshes(geometries, source_container: str = "Elements"):
    """Model parts become surfaces; meshes pass through."""
    meshes = []
    for geometry in geometries:
        if isinstance(geometry, Kratos.ModelPart):
            meshes.append(SurfaceOfModelPart(geometry, source_container))
        else:
            meshes.append(geometry)
    return meshes


def FitGeometryGuard(guard, geometries, source_container: str = "Elements"):
    """Fits the guardrail on the training family.

    Args:
        guard: A CreateGeometryGuard result.
        geometries: Kratos model parts and/or physicsnemo surface meshes.
        source_container: Which container a model part is tessellated from.
    """
    meshes = _AsMeshes(geometries, source_container)
    if len(meshes) < 2:
        raise ValueError(
            "A geometry guardrail needs a FAMILY to learn a distribution from; "
            f"got {len(meshes)} geometry. Fit it on the geometries the surrogate was "
            "trained across.")

    complaint = CheckFamilySize(meshes)
    if complaint:
        Kratos.Logger.PrintWarning("GeometryGuardrail", complaint)

    guard.fit(meshes)
    return guard


def FeatureWidth(mesh) -> int:
    """How many shape descriptors upstream extracts from one surface."""
    _, extract_features, _ = _TryImportGeometryGuardrail()
    return int(numpy.asarray(extract_features(mesh)).reshape(-1).size)


def CheckFamilySize(meshes):
    """Why a family is too small to fit a density model on, or None.

    A TRAP worth the extra feature extraction: upstream's descriptor vector
    is 22-wide, and a GMM fitted on fewer geometries than that is
    under-determined. Its verdicts then stop meaning anything - the first
    family tried here, six boxes, put EVERY query at percentile 100 and
    rejected them all, its own members included - and nothing upstream
    complains. A guard that rejects everything looks exactly like a guard
    that is working.
    """
    feature_width = FeatureWidth(meshes[0])
    if len(meshes) >= feature_width:
        return None
    return (
        f"Fitting on {len(meshes)} geometries with {feature_width} shape descriptors: "
        "the density model is under-determined and will flag geometries it should "
        "accept, its own training family included. Fit on at least "
        f"{feature_width} geometries, or expect false rejections.")


def QueryGeometry(guard, geometry, source_container: str = "Elements") -> dict:
    """The guardrail's verdict on one geometry.

    Returns:
        {"percentile": float, "status": "OK" | "WARN" | "REJECT"}.
    """
    return guard.query(_AsMeshes([geometry], source_container))[0]


def SaveGeometryGuard(guard, path) -> str:
    """Writes the fitted guardrail as a ``.npz`` sidecar.

    The extension is REQUIRED, and not as a style rule: upstream's save is
    a bare ``numpy.savez``, which appends ``.npz`` to any path lacking it,
    while its load is a bare ``numpy.load`` that does not. A guard saved as
    "guard.pt" therefore lands in "guard.pt.npz", reports success, and can
    never be loaded back under the name it was given - the same trap
    physicsnemo Module.save has with ``.mdlus``.
    """
    from pathlib import Path

    path = Path(str(path))
    if path.suffix != ".npz":
        raise ValueError(
            f"A geometry guard sidecar must be named \"*.npz\" [ path = {path} ]. "
            "Upstream saves with numpy.savez, which would silently append the "
            f"extension and write \"{path}.npz\" instead - a file its own loader "
            "would not find.")
    guard.save(path)
    return str(path)


def LoadGeometryGuard(path, device: str = "cpu"):
    """Reads a sidecar written by SaveGeometryGuard."""
    from pathlib import Path

    GeometryGuardrail, _, _ = _TryImportGeometryGuardrail()
    return GeometryGuardrail.load(Path(str(path)), device=device)


class GeometryGuardCheck:
    """Per-process geometry guardrail, with the field guard's policies.

    Settings:
        guard_file: The SaveGeometryGuard sidecar ("" disables it).
        policy: "advisory" (default; a Kratos warning per flagged
            geometry), "strict" (raise on REJECT) or "ignore".
        source_container: Which container the boundary is built from.
        reject_on_warn: Treat "WARN" as flagged too (default false, so
            only "REJECT" is flagged).

    The geometry is re-checked only when the node count changes, since
    rebuilding the surface every step would cost more than the inference.
    ``last_status`` holds the last verdict.
    """

    def __init__(self, settings: Kratos.Parameters) -> None:
        default_settings = Kratos.Parameters("""{
            "guard_file"       : "",
            "policy"           : "advisory",
            "source_container" : "Elements",
            "reject_on_warn"   : false
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        self.guard_file = settings["guard_file"].GetString()
        self.policy = settings["policy"].GetString()
        if self.policy not in _GUARD_POLICIES:
            raise ValueError(
                f"Unsupported geometry guard policy \"{self.policy}\". "
                f"Use one of {_GUARD_POLICIES}.")
        self.source_container = settings["source_container"].GetString()
        self.reject_on_warn = settings["reject_on_warn"].GetBool()
        self.last_status = None
        self.last_percentile = None
        self._guard = None
        self._checked_node_count = None

    @property
    def enabled(self) -> bool:
        return bool(self.guard_file) and self.policy != "ignore"

    def Check(self, model_part: Kratos.ModelPart, tag: str) -> bool:
        """Checks the model part's shape; True when flagged."""
        if not self.enabled:
            return False
        node_count = model_part.NumberOfNodes()
        if self._checked_node_count == node_count:
            return bool(self.last_status in self._FlaggedStatuses())
        if self._guard is None:
            self._guard = LoadGeometryGuard(self.guard_file)

        verdict = QueryGeometry(self._guard, model_part, self.source_container)
        self._checked_node_count = node_count
        self.last_status = verdict["status"]
        self.last_percentile = verdict["percentile"]
        flagged = self.last_status in self._FlaggedStatuses()
        if flagged:
            message = (
                f"geometry flagged as {self.last_status} at percentile "
                f"{self.last_percentile:.2f} of the training family")
            if self.policy == "strict":
                raise RuntimeError(
                    f"{tag}: {message}; the policy is \"strict\", so execution stops.")
            Kratos.Logger.PrintWarning(
                tag, f"The {message}; the guard is advisory, so execution continues.")
        return flagged

    def _FlaggedStatuses(self):
        return ("REJECT", "WARN") if self.reject_on_warn else ("REJECT",)
