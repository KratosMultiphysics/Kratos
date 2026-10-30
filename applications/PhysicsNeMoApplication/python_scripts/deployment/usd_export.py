"""OpenUSD (digital-twin) export of mesh series and fields.

Writes time-sampled USD stages - the interchange format Omniverse, usdview
and every USD-aware DCC tool read - so a running solve (or a deployed
surrogate writing into ordinary Kratos variables) becomes an interactive,
scrubbable 3D asset. This module is deliberately array-based and torch-free:
it takes numpy points/faces/fields and needs only ``usd-core`` (the
self-contained PyPI build of Pixar's USD - no Omniverse install), imported
lazily per the dependency policy.

The layout written is the plain UsdGeom one every viewer understands:

- one ``UsdGeomMesh`` (or ``UsdGeomPoints`` for meshless clouds) per prim,
- ``points`` time-sampled per step (deforming geometry plays back directly),
- topology (``faceVertexIndices``/``faceVertexCounts``) time-sampled ONLY on
  the steps where it actually changes, so adaptive-remeshing series stay
  valid and everything else stays compact,
- every exported field a time-sampled ``primvars:<NAME>`` with "vertex"
  interpolation (scalars as float[], 3-vectors as float3[], other widths as
  float[] with elementSize), which is what viewers color by.

UsdExportProcess (processes/export) drives this from the solution loop.
"""

import numpy


def _TryImportUsd():
    try:
        from pxr import Usd, UsdGeom, Sdf, Vt
        return Usd, UsdGeom, Sdf, Vt
    except ImportError as e:
        raise ImportError(
            "PhysicsNeMoApplication.usd_export requires OpenUSD, which could not be "
            "imported. Install it with e.g. 'pip install usd-core'.") from e


def CreateUsdStage(path, up_axis: str = "Z", meters_per_unit: float = 1.0,
                   time_codes_per_second: float = 1.0):
    """Creates (overwriting) a USD stage configured for a time-sampled export.

    Args:
        path: The stage file. ".usda" writes readable text, ".usd"/".usdc"
            the binary crate format - the extension picks the format.
        up_axis: "Z" (Kratos's own convention, the default) or "Y".
        meters_per_unit: Scene scale metadata (1.0 = meters).
        time_codes_per_second: Playback rate: how many time codes make one
            second. With the step number as the time code, this is the
            solver steps per second of playback; with solver TIME as the
            time code, 1.0 plays back in real time.

    Returns:
        The Usd.Stage. Callers write samples with WriteMeshTimeSample /
        WritePointsTimeSample and persist with SaveStage.
    """
    Usd, UsdGeom, _, _ = _TryImportUsd()

    if up_axis not in ("Y", "Z"):
        raise ValueError(f"\"up_axis\" must be \"Y\" or \"Z\", got \"{up_axis}\".")

    from pathlib import Path
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()  # CreateNew refuses to clobber an existing layer
    stage = Usd.Stage.CreateNew(str(path))
    UsdGeom.SetStageUpAxis(stage, UsdGeom.Tokens.z if up_axis == "Z" else UsdGeom.Tokens.y)
    UsdGeom.SetStageMetersPerUnit(stage, float(meters_per_unit))
    stage.SetTimeCodesPerSecond(float(time_codes_per_second))
    return stage


def _AsVec3fArray(Vt, points):
    points = numpy.ascontiguousarray(points, dtype=numpy.float32)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"Expected (N, 3) points, got shape {points.shape}.")
    return Vt.Vec3fArray.FromNumpy(points)


def _UpdateTimeRange(stage, time_code):
    start = stage.GetStartTimeCode() if stage.HasAuthoredTimeCodeRange() else time_code
    end = stage.GetEndTimeCode() if stage.HasAuthoredTimeCodeRange() else time_code
    stage.SetStartTimeCode(min(start, time_code))
    stage.SetEndTimeCode(max(end, time_code))


def _SetExtent(UsdGeom, Vt, geom, points, time_code):
    lo = points.min(axis=0).astype(numpy.float32)
    hi = points.max(axis=0).astype(numpy.float32)
    geom.GetExtentAttr().Set(Vt.Vec3fArray.FromNumpy(numpy.stack([lo, hi])), time_code)


def _WriteFieldPrimvars(UsdGeom, Sdf, Vt, geom, point_fields, n_points, time_code):
    primvars = UsdGeom.PrimvarsAPI(geom.GetPrim())
    for name, values in point_fields.items():
        values = numpy.ascontiguousarray(values, dtype=numpy.float32)
        if values.ndim == 1:
            values = values[:, None]
        if values.ndim != 2 or values.shape[0] != n_points:
            raise ValueError(
                f"Field \"{name}\" must be (n_points,) or (n_points, width) with "
                f"n_points = {n_points}, got shape {values.shape}.")
        width = values.shape[1]
        if width == 3:
            primvar = primvars.CreatePrimvar(
                name, Sdf.ValueTypeNames.Float3Array, UsdGeom.Tokens.vertex)
            primvar.Set(Vt.Vec3fArray.FromNumpy(values), time_code)
        else:
            primvar = primvars.CreatePrimvar(
                name, Sdf.ValueTypeNames.FloatArray, UsdGeom.Tokens.vertex)
            if width != 1:
                primvar.SetElementSize(width)
            primvar.Set(Vt.FloatArray.FromNumpy(values.reshape(-1)), time_code)


def WriteMeshTimeSample(stage, prim_path: str, points, triangles, point_fields,
                        time_code: float):
    """Writes one time sample of a triangle-surface mesh with nodal fields.

    The prim is defined on first use. Points, extent and field primvars are
    time-sampled every call; the topology attributes are time-sampled only
    when the triangle array differs from the previously written one, so a
    fixed mesh carries a single topology sample while an adaptively
    remeshed series stays correct frame by frame.

    Args:
        stage: A stage from CreateUsdStage.
        prim_path: Absolute prim path, e.g. "/Kratos/Structure".
        points: (N, 3) array; unreferenced points are allowed (a volume
            mesh's interior nodes keep field arrays aligned with N).
        triangles: (F, 3) int array of point indices.
        point_fields: {name: (N,) | (N, width) array} written as
            "vertex"-interpolated primvars (see the module docstring).
        time_code: The sample's time code (step number or solver time -
            the process decides, CreateUsdStage documents playback).
    """
    Usd, UsdGeom, Sdf, Vt = _TryImportUsd()

    points = numpy.ascontiguousarray(points, dtype=numpy.float64)
    triangles = numpy.ascontiguousarray(triangles, dtype=numpy.int64)
    if triangles.ndim != 2 or triangles.shape[1] != 3:
        raise ValueError(f"Expected (F, 3) triangles, got shape {triangles.shape}.")

    mesh = UsdGeom.Mesh.Define(stage, prim_path)
    mesh.GetPointsAttr().Set(_AsVec3fArray(Vt, points), time_code)
    _SetExtent(UsdGeom, Vt, mesh, points, time_code)

    topology_key = "physicsNemo:lastTopologyHash"
    topology_hash = str(hash(triangles.tobytes()))
    if mesh.GetPrim().GetCustomDataByKey(topology_key) != topology_hash:
        mesh.GetFaceVertexIndicesAttr().Set(
            Vt.IntArray.FromNumpy(triangles.reshape(-1).astype(numpy.int32)), time_code)
        mesh.GetFaceVertexCountsAttr().Set(
            Vt.IntArray.FromNumpy(numpy.full(len(triangles), 3, dtype=numpy.int32)),
            time_code)
        mesh.GetPrim().SetCustomDataByKey(topology_key, topology_hash)

    _WriteFieldPrimvars(UsdGeom, Sdf, Vt, mesh, point_fields, len(points), time_code)
    _UpdateTimeRange(stage, time_code)


def WritePointsTimeSample(stage, prim_path: str, points, point_fields,
                          time_code: float, widths=None):
    """Writes one time sample of a point cloud (UsdGeomPoints) with fields.

    The meshless counterpart of WriteMeshTimeSample, for particle clouds
    (MPM/SPH/DEM-style parts, ParticleInferenceProcess output).

    Args:
        stage: A stage from CreateUsdStage.
        prim_path: Absolute prim path.
        points: (N, 3) array.
        point_fields: {name: (N,) | (N, width) array} vertex primvars.
        time_code: The sample's time code.
        widths: Optional (N,) per-point display diameters.
    """
    Usd, UsdGeom, Sdf, Vt = _TryImportUsd()

    points = numpy.ascontiguousarray(points, dtype=numpy.float64)
    cloud = UsdGeom.Points.Define(stage, prim_path)
    cloud.GetPointsAttr().Set(_AsVec3fArray(Vt, points), time_code)
    _SetExtent(UsdGeom, Vt, cloud, points, time_code)
    if widths is not None:
        widths = numpy.ascontiguousarray(widths, dtype=numpy.float32).reshape(-1)
        if widths.shape[0] != points.shape[0]:
            raise ValueError(
                f"\"widths\" must have one entry per point, got {widths.shape[0]} "
                f"for {points.shape[0]} points.")
        cloud.GetWidthsAttr().Set(Vt.FloatArray.FromNumpy(widths), time_code)

    _WriteFieldPrimvars(UsdGeom, Sdf, Vt, cloud, point_fields, len(points), time_code)
    _UpdateTimeRange(stage, time_code)


def SaveStage(stage):
    """Persists the stage to its file (the root layer's own path)."""
    stage.GetRootLayer().Save()
