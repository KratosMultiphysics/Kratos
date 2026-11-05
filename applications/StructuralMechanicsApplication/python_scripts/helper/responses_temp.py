"""Stress responses for small-displacement, linear-elastic structural elements."""

from itertools import product
from math import comb, sqrt

import numpy as np

import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA


def get_von_mises(center_gp, center_average, center_interpolated, model_part):
    """Return one von Mises stress per element as ``{element_id: value}``.

    Call after the linear analysis has calculated the element displacements.
    Select exactly one of the three Boolean options:

    * center_gp: use the GP at the mapped parametric center; raise if absent.
    * center_average: take the arithmetic mean of the GP stress tensors first.
      This is a representative element value, not a volume-weighted average.
    * center_interpolated: use a center GP if present, otherwise recover the
      stress tensor at the parametric center using polynomial interpolation
      (least squares for oversampled simplex rules) in local coordinates.

    Supports standard SmallDisplacementElement and SmallDisplacementBbarElement
    continua with built-in linear-elastic laws, and homogeneous linear shells:
    ShellThinElement3D3N and LegacyThickShellElement3D4N.
    Shell stresses are evaluated at the middle surface. Plane-strain sigma_zz is
    reconstructed for isotropic linear elasticity with no initial material state.
    Current MITC/CS-DSG3 shells, 1D elements and nonlinear formulations are not
    supported. Results are returned without writing values to the model part.
    """
    modes = (center_gp, center_average, center_interpolated)
    if any(not isinstance(mode, bool) for mode in modes):
        raise TypeError("The three center options must be Booleans.")
    if sum(modes) != 1:
        raise ValueError("Select exactly one center option.")

    results = {}
    for element in model_part.Elements:
        name = KM.CompareElementsAndConditionsUtility.GetRegisteredName(element)
        stresses = _get_stress_tensors(element, name, model_part.ProcessInfo)
        if center_average:
            stress = np.mean(stresses, axis=0)
        else:
            geometry = element.GetGeometry()
            center = _parametric_center(geometry)
            gp_coordinates = np.asarray(element.GetIntegrationPoints())
            if len(gp_coordinates) != len(stresses):
                raise RuntimeError(
                    f"Element {element.Id} ({name}): GP coordinates and stress "
                    "output have different lengths."
                )
            center_global = _global_coordinates(geometry, center)
            tolerance = 1e-10 * _geometry_scale(geometry)
            indices = np.flatnonzero(
                np.linalg.norm(gp_coordinates - center_global, axis=1) <= tolerance
            )
            if len(indices):
                stress = stresses[indices[0]]
            elif center_gp:
                raise ValueError(
                    f"Element {element.Id} ({name}) has no GP at its center. "
                    "Select center_interpolated or center_average."
                )
            else:
                local_points = _local_coordinates(geometry, gp_coordinates, center)
                stress = _interpolate_at_center(
                    geometry, local_points, center, stresses
                )
        xx, yy, zz = np.diag(stress)
        xy, yz, xz = stress[0, 1], stress[1, 2], stress[0, 2]
        results[element.Id] = sqrt(
            0.5 * ((xx - yy)**2 + (yy - zz)**2 + (zz - xx)**2)
            + 3.0 * (xy**2 + yz**2 + xz**2)
        )
    return results


def _get_stress_tensors(element, name, process_info):
    geometry = element.GetGeometry()
    dimensions = (geometry.LocalSpaceDimension(), geometry.WorkingSpaceDimension())
    plane_strain = False
    if name.startswith(("SmallDisplacementElement", "SmallDisplacementBbarElement")):
        law = element.Properties[KM.CONSTITUTIVE_LAW]
        if dimensions == (2, 2):
            if not isinstance(law, (SMA.LinearElasticPlaneStress2DLaw,
                                    SMA.LinearElasticPlaneStrain2DLaw)):
                raise NotImplementedError(
                    f"Element {element.Id}: only the built-in isotropic linear "
                    "plane-stress and plane-strain laws are supported in 2D."
                )
            plane_strain = isinstance(law, SMA.LinearElasticPlaneStrain2DLaw)
        elif dimensions == (3, 3):
            if not isinstance(law, (SMA.LinearElastic3DLaw,
                                    SMA.FlexibleLinearElastic3DLaw,
                                    SMA.UserProvidedLinearElastic3DLaw)):
                raise NotImplementedError(
                    f"Element {element.Id}: unsupported 3D material law. "
                    "Use a built-in linear-elastic law."
                )
        else:
            raise NotImplementedError(f"Unsupported continuum geometry: {name}.")
        variable = KM.CAUCHY_STRESS_TENSOR
    elif name in ("ShellThinElement3D3N", "LegacyThickShellElement3D4N"):
        if (element.Properties.Has(SMA.SHELL_ORTHOTROPIC_LAYERS)
                or not isinstance(element.Properties[KM.CONSTITUTIVE_LAW],
                                  (SMA.LinearElasticPlaneStress2DLaw,
                                   SMA.LinearElastic3DLaw))):
            raise NotImplementedError(
                f"Element {element.Id}: only homogeneous isotropic linear shells "
                "are supported."
            )
        variable = SMA.SHELL_STRESS_MIDDLE_SURFACE_GLOBAL
    elif name in ("MITCThickShellElement3D4N", "ShellThickElement3D4N",
                  "CSDSG3ThickShellLinearElement3D3N"):
        raise NotImplementedError(
            f"Element {element.Id} ({name}) does not expose the physical stress "
            "tensor required for this Python-only response. Its scalar von Mises "
            "output cannot be used for tensor averaging or interpolation."
        )
    else:
        raise NotImplementedError(
            f"Element {element.Id} ({name}) is unsupported; 1D elements and "
            "nonlinear formulations are not implemented."
        )

    raw_stresses = element.CalculateOnIntegrationPoints(variable, process_info)
    if not raw_stresses:
        raise RuntimeError(f"Element {element.Id} returned no GP stress tensors.")
    stresses = np.zeros((len(raw_stresses), 3, 3))
    for index, tensor in enumerate(raw_stresses):
        size = tensor.Size1()
        if size not in (2, 3) or tensor.Size2() != size:
            raise RuntimeError(
                f"Element {element.Id}: expected a 2x2 or 3x3 stress tensor."
            )
        stresses[index, :size, :size] = [
            [tensor[i, j] for j in range(size)] for i in range(size)
        ]
        if plane_strain:
            stresses[index, 2, 2] = element.Properties[KM.POISSON_RATIO] * (
                stresses[index, 0, 0] + stresses[index, 1, 1]
            )
    if (not np.all(np.isfinite(stresses))
            or not np.allclose(stresses, stresses.transpose(0, 2, 1))):
        raise RuntimeError(f"Element {element.Id} returned invalid stress tensors.")
    return stresses


def _parametric_center(geometry):
    family = KM.GeometryData.KratosGeometryFamily
    centers = {
        family.Kratos_Triangle: (1.0 / 3.0, 1.0 / 3.0, 0.0),
        family.Kratos_Quadrilateral: (0.0, 0.0, 0.0),
        family.Kratos_Tetrahedra: (0.25, 0.25, 0.25),
        family.Kratos_Hexahedra: (0.0, 0.0, 0.0),
        family.Kratos_Prism: (1.0 / 3.0, 1.0 / 3.0, 0.0),
        family.Kratos_Pyramid: (0.0, 0.0, -0.5),
    }
    try:
        return np.array(centers[geometry.GetGeometryFamily()])
    except KeyError:
        raise NotImplementedError(
            f"Center recovery is not implemented for {geometry.Info()}."
        ) from None


def _global_coordinates(geometry, local_point):
    return np.array(geometry.GlobalCoordinates(KM.Array3(local_point.tolist())))


def _geometry_scale(geometry):
    coordinates = np.array([[node.X, node.Y, node.Z] for node in geometry])
    return np.linalg.norm(np.ptp(coordinates, axis=0))


def _local_coordinates(geometry, global_points, center):
    """Invert GlobalCoordinates; the general C++ inverse is not bound to Python."""
    dimension = geometry.LocalSpaceDimension()
    tolerance = 1e-10 * _geometry_scale(geometry)
    local_points = []
    step = 1e-5
    for global_point in global_points:
        local = center.copy()
        for _ in range(25):
            residual = global_point - _global_coordinates(geometry, local)
            if np.linalg.norm(residual) <= tolerance:
                local_points.append(local[:dimension].copy())
                break
            columns = []
            for axis in range(dimension):
                offset = np.zeros(3)
                offset[axis] = step
                columns.append((
                    _global_coordinates(geometry, local + offset)
                    - _global_coordinates(geometry, local - offset)
                ) / (2.0 * step))
            increment, _, rank, _ = np.linalg.lstsq(
                np.column_stack(columns), residual, rcond=None
            )
            if rank != dimension:
                raise ValueError("Cannot invert a degenerate element geometry.")
            local[:dimension] += increment
        else:
            raise ValueError("Mapping a GP to local coordinates did not converge.")
    return np.array(local_points)


def _simplex_powers(dimension, number_of_points):
    degree = 0
    while comb(degree + 1 + dimension, dimension) <= number_of_points:
        degree += 1
    return [powers for powers in product(range(degree + 1), repeat=dimension)
            if sum(powers) <= degree]


def _interpolate_at_center(geometry, local_points, center, stresses):
    number_of_points, dimension = local_points.shape
    if number_of_points == 1:
        raise ValueError("A single off-center GP cannot determine a center value.")
    family = KM.GeometryData.KratosGeometryFamily
    geometry_family = geometry.GetGeometryFamily()
    if geometry_family in (family.Kratos_Quadrilateral, family.Kratos_Hexahedra):
        counts = [len(np.unique(np.round(local_points[:, axis], 8)))
                  for axis in range(dimension)]
        if np.prod(counts) != number_of_points:
            raise NotImplementedError("Expected a tensor-product GP arrangement.")
        powers = list(product(*(range(count) for count in counts)))
    elif geometry_family == family.Kratos_Prism:
        z_count = len(np.unique(np.round(local_points[:, 2], 8)))
        if number_of_points % z_count:
            raise NotImplementedError("Expected a triangle-times-line GP arrangement.")
        powers = [(*xy, z) for xy in _simplex_powers(2, number_of_points // z_count)
                  for z in range(z_count)]
    else:
        powers = _simplex_powers(dimension, number_of_points)

    relative_points = local_points - center[:dimension]
    matrix = np.array([
        np.prod(relative_points**np.array(exponent), axis=1)
        for exponent in powers
    ]).T
    coefficients, _, rank, _ = np.linalg.lstsq(
        matrix, stresses.reshape(number_of_points, 9), rcond=None
    )
    if rank != len(powers) or len(powers) < dimension + 1:
        raise ValueError("The GP arrangement cannot determine the center polynomial.")
    # At the shifted origin only the constant monomial (the first column) survives.
    return coefficients[0].reshape(3, 3)
