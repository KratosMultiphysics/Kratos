#!/usr/bin/env python3
"""Static Stokes convergence driver for the FEM Gap-SBM workflow."""

from __future__ import annotations

import math
import os
from pathlib import Path
import sys
import time


EXAMPLE_DIRECTORY = Path(__file__).resolve().parent
BASE_EXAMPLE_DIRECTORY = (
    EXAMPLE_DIRECTORY.parent / "fluid_fem_sbm_manufactured_solution"
)
GEOMETRY = os.environ.get("GAP_SBM_GEOMETRY", "circle").lower()
if GEOMETRY not in ("circle", "diamond"):
    raise ValueError("GAP_SBM_GEOMETRY must be either 'circle' or 'diamond'.")
RESULTS_DIRECTORY = EXAMPLE_DIRECTORY / (
    "results" if GEOMETRY == "circle" else "results_diamond"
)
RESULTS_DIRECTORY.mkdir(parents=True, exist_ok=True)

os.environ["OMP_NUM_THREADS"] = "4"
os.environ["OMP_THREAD_LIMIT"] = "4"
sys.path.insert(0, str(BASE_EXAMPLE_DIRECTORY))

import KratosMultiphysics as KM

import run_convergence as study


# Reuse the verified divergence-free sin/cos manufactured solution at t=0.
# Its analytical acceleration is zero at that instant, hence the same forcing
# exactly represents the static Stokes equations when force_steady_state=true.
_transient_exact_velocity = study.exact_velocity
_transient_exact_velocity_gradient = study.exact_velocity_gradient
_transient_exact_acceleration = study.exact_acceleration
_transient_exact_pressure = study.exact_pressure
_transient_body_force = study.manufactured_body_force
_background_error_norms = study._calculate_error_norms
_base_calculate_rates = study._calculate_rates
_base_plot_convergence = study._plot_convergence
_circle_true_skin_factory = study._create_true_skin

DIAMOND_CENTER = (0.5, 0.5)
DIAMOND_HALF_DIAGONAL = 0.25
DIAMOND_FLUID_AREA = 1.0 - 2.0 * DIAMOND_HALF_DIAGONAL**2

# Filled once per mesh level by _calculate_full_domain_error_norms. Keeping
# these values separate avoids changing the CSV format of the reusable
# classical-SBM convergence driver.
_BACKGROUND_ONLY_ERRORS: list[dict[str, float]] = []
_ERROR_COMPARISON_PLOT: Path | None = None
_INTERFACE_DIAGNOSTIC_PLOT: Path | None = None
_FULL_DOMAIN_AREAS: list[float] = []
_DIAGNOSTIC_ROWS: list[dict[str, float]] = []
_ENABLE_INTERFACE_DIAGNOSTICS = False


def _consume_interface_diagnostics_flag() -> None:
    """Enable the optional E_Gamma/J_u/J_t study and hide its local flag."""
    global _ENABLE_INTERFACE_DIAGNOSTICS

    flag = "--interface-diagnostics"
    _ENABLE_INTERFACE_DIAGNOSTICS = flag in sys.argv
    if _ENABLE_INTERFACE_DIAGNOSTICS:
        sys.argv[:] = [argument for argument in sys.argv if argument != flag]


def _static_exact_velocity(x: float, y: float, _: float):
    return _transient_exact_velocity(x, y, 0.0)


def _static_exact_velocity_gradient(x: float, y: float, _: float):
    return _transient_exact_velocity_gradient(x, y, 0.0)


def _static_exact_acceleration(x: float, y: float, _: float):
    return _transient_exact_acceleration(x, y, 0.0)


def _static_exact_pressure(x: float, y: float, _: float):
    return _transient_exact_pressure(x, y, 0.0)


def _static_body_force(x: float, y: float, _: float):
    return _transient_body_force(x, y, 0.0)


def _set_gap_manufactured_solution(name: str) -> None:
    """Select the analytical fields and keep each study's output separate."""
    global RESULTS_DIRECTORY

    if name == "sin-cos":
        study.exact_velocity = _static_exact_velocity
        study.exact_velocity_gradient = _static_exact_velocity_gradient
        study.exact_acceleration = _static_exact_acceleration
        study.exact_pressure = _static_exact_pressure
        study.manufactured_body_force = _static_body_force
    elif name == "linear":
        study.exact_velocity = study.linear_exact_velocity
        study.exact_velocity_gradient = study.linear_exact_velocity_gradient
        study.exact_acceleration = study.linear_exact_acceleration
        study.exact_pressure = study.linear_exact_pressure
        study.manufactured_body_force = study.linear_manufactured_body_force
    else:
        raise ValueError(f"Unknown manufactured solution: {name}")

    directory_suffix = "" if name == "sin-cos" else "_linear"
    geometry_suffix = "" if GEOMETRY == "circle" else "_diamond"
    RESULTS_DIRECTORY = EXAMPLE_DIRECTORY / (
        f"results{geometry_suffix}{directory_suffix}"
    )
    RESULTS_DIRECTORY.mkdir(parents=True, exist_ok=True)
    study.RESULTS_DIRECTORY = RESULTS_DIRECTORY
    study.PROBLEM_LABEL = (
        f"Static Stokes Gap-SBM ({GEOMETRY}, {name})"
    )


def _create_diamond_true_skin(model: KM.Model) -> KM.ModelPart:
    """Create an exact four-segment diamond hole."""
    skin = model.CreateModelPart("TrueSkin")
    properties = skin.CreateNewProperties(0)
    center_x, center_y = DIAMOND_CENTER
    half_diagonal = DIAMOND_HALF_DIAGONAL
    coordinates = (
        (center_x, center_y - half_diagonal),
        (center_x + half_diagonal, center_y),
        (center_x, center_y + half_diagonal),
        (center_x - half_diagonal, center_y),
    )
    for node_id, (x, y) in enumerate(coordinates, start=1):
        skin.CreateNewNode(node_id, x, y, 0.0)
    for condition_id in range(1, 5):
        skin.CreateNewCondition(
            "LineCondition2D2N",
            condition_id,
            [condition_id, condition_id % 4 + 1],
            properties,
        )
    return skin


def _create_selected_true_skin(model: KM.Model) -> KM.ModelPart:
    if GEOMETRY == "diamond":
        return _create_diamond_true_skin(model)
    return _circle_true_skin_factory(model)


def _reconstruct_boundary_endpoints(boundary_quadrature):
    """Recover the Line2 endpoints from its two Gauss-point coordinates."""
    inverse_gauss_coordinate = math.sqrt(3.0)
    gauss_points = (
        (boundary_quadrature[0, 0], boundary_quadrature[0, 1]),
        (boundary_quadrature[1, 0], boundary_quadrature[1, 1]),
    )
    midpoint = tuple(
        0.5 * (gauss_points[0][d] + gauss_points[1][d]) for d in range(2)
    )
    delta = tuple(
        inverse_gauss_coordinate * (gauss_points[1][d] - gauss_points[0][d])
        for d in range(2)
    )
    return (
        tuple(midpoint[d] - 0.5 * delta[d] for d in range(2)),
        tuple(midpoint[d] + 0.5 * delta[d] for d in range(2)),
    )


def _evaluate_extended_parent_field(geometry, x: float, y: float):
    """Evaluate the parent Q1 velocity, pressure and total Cauchy stress."""
    x_min = min(node.X for node in geometry)
    x_max = max(node.X for node in geometry)
    y_min = min(node.Y for node in geometry)
    y_max = max(node.Y for node in geometry)
    cell_size_x = x_max - x_min
    cell_size_y = y_max - y_min
    if cell_size_x <= 0.0 or cell_size_y <= 0.0:
        raise RuntimeError("Gap diagnostics require a valid Cartesian Quad4 parent.")

    xi = 2.0 * (x - x_min) / cell_size_x - 1.0
    eta = 2.0 * (y - y_min) / cell_size_y - 1.0
    shape_functions = (
        0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta),
    )
    derivatives_x = (
        -0.5 * (1.0 - eta) / cell_size_x,
        0.5 * (1.0 - eta) / cell_size_x,
        0.5 * (1.0 + eta) / cell_size_x,
        -0.5 * (1.0 + eta) / cell_size_x,
    )
    derivatives_y = (
        -0.5 * (1.0 - xi) / cell_size_y,
        -0.5 * (1.0 + xi) / cell_size_y,
        0.5 * (1.0 + xi) / cell_size_y,
        0.5 * (1.0 - xi) / cell_size_y,
    )
    velocity_x = [
        geometry[i].GetSolutionStepValue(KM.VELOCITY_X) for i in range(4)
    ]
    velocity_y = [
        geometry[i].GetSolutionStepValue(KM.VELOCITY_Y) for i in range(4)
    ]
    pressure_values = [
        geometry[i].GetSolutionStepValue(KM.PRESSURE) for i in range(4)
    ]
    velocity = (
        sum(shape_functions[i] * velocity_x[i] for i in range(4)),
        sum(shape_functions[i] * velocity_y[i] for i in range(4)),
    )
    pressure = sum(
        shape_functions[i] * pressure_values[i] for i in range(4)
    )
    du_x_dx = sum(derivatives_x[i] * velocity_x[i] for i in range(4))
    du_x_dy = sum(derivatives_y[i] * velocity_x[i] for i in range(4))
    du_y_dx = sum(derivatives_x[i] * velocity_y[i] for i in range(4))
    du_y_dy = sum(derivatives_y[i] * velocity_y[i] for i in range(4))

    # Match Newtonian2DLaw exactly, including its three-dimensional deviatoric
    # correction applied to the 2D strain-rate trace.
    viscosity = study.DYNAMIC_VISCOSITY
    volumetric_part = (du_x_dx + du_y_dy) / 3.0
    viscous_xx = 2.0 * viscosity * (du_x_dx - volumetric_part)
    viscous_yy = 2.0 * viscosity * (du_y_dy - volumetric_part)
    viscous_xy = viscosity * (du_x_dy + du_y_dx)
    total_stress = (
        (viscous_xx - pressure, viscous_xy),
        (viscous_xy, viscous_yy - pressure),
    )
    return velocity, total_stress


def _distance(first, second) -> float:
    return math.hypot(first[0] - second[0], first[1] - second[1])


def _calculate_boundary_and_interface_diagnostics(
    fluid: KM.ModelPart, final_time: float
) -> dict[str, float]:
    """Measure the true-boundary error and jumps between gap extensions."""
    gap_elements = list(fluid.GetSubModelPart("GapElements").Elements)
    surrogate_conditions = list(
        fluid.GetSubModelPart(study.SURROGATE_SUB_MODEL_PART_NAME).Conditions
    )
    if len(gap_elements) != len(surrogate_conditions):
        raise RuntimeError(
            "Gap diagnostics require one gap element per surrogate condition."
        )

    boundary_error_squared = 0.0
    boundary_length = 0.0
    endpoint_records: dict[tuple[float, float], list[dict[str, object]]] = {}
    for element, condition in zip(gap_elements, surrogate_conditions):
        geometry = element.GetGeometry()
        face = condition.GetValue(KM.SURROGATE_BOUNDARY_FACE_COORDINATES)
        boundary_quadrature = condition.GetValue(
            KM.SURROGATE_BOUNDARY_PROJECTION
        )
        if boundary_quadrature.Size1() != 2 or boundary_quadrature.Size2() != 7:
            raise RuntimeError("Expected Gap-SBM 2x7 boundary quadrature data.")

        true_endpoints = _reconstruct_boundary_endpoints(boundary_quadrature)
        surrogate_endpoints = (
            (face[0, 0], face[0, 1]),
            (face[1, 0], face[1, 1]),
        )
        for endpoint in range(2):
            key = tuple(round(value, 12) for value in surrogate_endpoints[endpoint])
            endpoint_records.setdefault(key, []).append(
                {
                    "true_point": true_endpoints[endpoint],
                    "geometry": geometry,
                }
            )

        for g in range(2):
            weight = boundary_quadrature[g, 6]
            if weight <= 0.0:
                continue
            x = boundary_quadrature[g, 0]
            y = boundary_quadrature[g, 1]
            velocity, _ = _evaluate_extended_parent_field(geometry, x, y)
            prescribed_velocity = study.exact_velocity(x, y, final_time)
            boundary_error_squared += weight * (
                (velocity[0] - prescribed_velocity[0]) ** 2
                + (velocity[1] - prescribed_velocity[1]) ** 2
            )
            boundary_length += weight

    velocity_jump_squared = 0.0
    traction_jump_squared = 0.0
    interface_length = 0.0
    gauss_points, gauss_weights = study.np.polynomial.legendre.leggauss(2)
    tolerance = 1.0e-9
    interfaces = fluid.GetSubModelPart("GapInterfaces").Conditions
    for condition in interfaces:
        connector = condition.GetValue(KM.SURROGATE_BOUNDARY_FACE_COORDINATES)
        connector_points = (
            (connector[0, 0], connector[0, 1]),
            (connector[1, 0], connector[1, 1]),
        )

        matched_records = None
        for surrogate_endpoint in range(2):
            key = tuple(
                round(value, 12)
                for value in connector_points[surrogate_endpoint]
            )
            candidates = endpoint_records.get(key, [])
            true_point = connector_points[1 - surrogate_endpoint]
            compatible = [
                record
                for record in candidates
                if _distance(record["true_point"], true_point) <= tolerance
            ]
            if len(compatible) == 2:
                matched_records = compatible
                break
        if matched_records is None:
            raise RuntimeError(
                f"Could not associate gap interface condition #{condition.Id} "
                "with its two parent extensions."
            )

        tangent = (
            connector_points[1][0] - connector_points[0][0],
            connector_points[1][1] - connector_points[0][1],
        )
        length = math.hypot(tangent[0], tangent[1])
        if length <= 0.0:
            raise RuntimeError(
                f"Gap interface condition #{condition.Id} has zero length."
            )
        normal = (tangent[1] / length, -tangent[0] / length)
        for local_coordinate, gauss_weight in zip(gauss_points, gauss_weights):
            first_weight = 0.5 * (1.0 - local_coordinate)
            second_weight = 0.5 * (1.0 + local_coordinate)
            x = (
                first_weight * connector_points[0][0]
                + second_weight * connector_points[1][0]
            )
            y = (
                first_weight * connector_points[0][1]
                + second_weight * connector_points[1][1]
            )
            integration_weight = 0.5 * length * gauss_weight
            side_fields = [
                _evaluate_extended_parent_field(record["geometry"], x, y)
                for record in matched_records
            ]
            velocity_jump = (
                side_fields[0][0][0] - side_fields[1][0][0],
                side_fields[0][0][1] - side_fields[1][0][1],
            )
            tractions = []
            for _, stress in side_fields:
                tractions.append(
                    (
                        stress[0][0] * normal[0] + stress[0][1] * normal[1],
                        stress[1][0] * normal[0] + stress[1][1] * normal[1],
                    )
                )
            traction_jump = (
                tractions[0][0] - tractions[1][0],
                tractions[0][1] - tractions[1][1],
            )
            velocity_jump_squared += integration_weight * sum(
                component**2 for component in velocity_jump
            )
            traction_jump_squared += integration_weight * sum(
                component**2 for component in traction_jump
            )
            interface_length += integration_weight

    return {
        "boundary_trace_l2": math.sqrt(max(0.0, boundary_error_squared)),
        "interface_velocity_jump_l2": math.sqrt(
            max(0.0, velocity_jump_squared)
        ),
        "interface_traction_jump_l2": math.sqrt(
            max(0.0, traction_jump_squared)
        ),
        "reconstructed_boundary_length": boundary_length,
        "gap_interface_length": interface_length,
    }


def _calculate_gap_error_integrals(
    fluid: KM.ModelPart, final_time: float
) -> dict[str, float]:
    """Integrate error contributions over all ruled gap patches."""
    integral_names = (
        "area",
        "velocity_error_squared",
        "velocity_norm_squared",
        "pressure_error_squared",
        "pressure_error_integral",
        "pressure_norm_squared",
        "h1_error_squared",
        "h1_norm_squared",
    )
    integrals = {name: 0.0 for name in integral_names}
    gap_elements = list(fluid.GetSubModelPart("GapElements").Elements)
    surrogate_conditions = list(
        fluid.GetSubModelPart(study.SURROGATE_SUB_MODEL_PART_NAME).Conditions
    )
    if len(gap_elements) != len(surrogate_conditions):
        raise RuntimeError(
            "Gap-SBM requires one gap element per surrogate condition."
        )

    gauss_points, gauss_weights = study.np.polynomial.legendre.leggauss(4)
    inverse_boundary_gauss_coordinate = math.sqrt(3.0)

    # Both containers are ID ordered. The modeler creates gap elements in the
    # same face order and additionally stores the explicit C++ association in
    # NEIGHBOUR_CONDITIONS.
    for element, condition in zip(gap_elements, surrogate_conditions):
        geometry = element.GetGeometry()
        face = condition.GetValue(KM.SURROGATE_BOUNDARY_FACE_COORDINATES)
        boundary_quadrature = condition.GetValue(
            KM.SURROGATE_BOUNDARY_PROJECTION
        )
        if boundary_quadrature.Size1() != 2 or boundary_quadrature.Size2() != 7:
            raise RuntimeError("Expected Gap-SBM 2x7 boundary quadrature data.")

        surrogate_points = (
            (face[0, 0], face[0, 1]),
            (face[1, 0], face[1, 1]),
        )
        boundary_gauss_points = (
            (boundary_quadrature[0, 0], boundary_quadrature[0, 1]),
            (boundary_quadrature[1, 0], boundary_quadrature[1, 1]),
        )
        boundary_midpoint = tuple(
            0.5 * (boundary_gauss_points[0][d] + boundary_gauss_points[1][d])
            for d in range(2)
        )
        boundary_delta = tuple(
            inverse_boundary_gauss_coordinate
            * (boundary_gauss_points[1][d] - boundary_gauss_points[0][d])
            for d in range(2)
        )
        true_points = (
            tuple(boundary_midpoint[d] - 0.5 * boundary_delta[d] for d in range(2)),
            tuple(boundary_midpoint[d] + 0.5 * boundary_delta[d] for d in range(2)),
        )

        x_min = min(node.X for node in geometry)
        x_max = max(node.X for node in geometry)
        y_min = min(node.Y for node in geometry)
        y_max = max(node.Y for node in geometry)
        cell_size_x = x_max - x_min
        cell_size_y = y_max - y_min
        if cell_size_x <= 0.0 or cell_size_y <= 0.0:
            raise RuntimeError("Gap-SBM error integration requires a valid Quad4 parent.")

        nodal_velocity_x = [
            geometry[i].GetSolutionStepValue(KM.VELOCITY_X) for i in range(4)
        ]
        nodal_velocity_y = [
            geometry[i].GetSolutionStepValue(KM.VELOCITY_Y) for i in range(4)
        ]
        nodal_pressure = [
            geometry[i].GetSolutionStepValue(KM.PRESSURE) for i in range(4)
        ]

        for xi, weight_x in zip(gauss_points, gauss_weights):
            first_xi = 0.5 * (1.0 - xi)
            second_xi = 0.5 * (1.0 + xi)
            surrogate_point = tuple(
                first_xi * surrogate_points[0][d]
                + second_xi * surrogate_points[1][d]
                for d in range(2)
            )
            true_point = tuple(
                first_xi * true_points[0][d] + second_xi * true_points[1][d]
                for d in range(2)
            )
            surrogate_tangent = tuple(
                0.5 * (surrogate_points[1][d] - surrogate_points[0][d])
                for d in range(2)
            )
            true_tangent = tuple(
                0.5 * (true_points[1][d] - true_points[0][d])
                for d in range(2)
            )

            for eta, weight_y in zip(gauss_points, gauss_weights):
                first_eta = 0.5 * (1.0 - eta)
                second_eta = 0.5 * (1.0 + eta)
                x = first_eta * surrogate_point[0] + second_eta * true_point[0]
                y = first_eta * surrogate_point[1] + second_eta * true_point[1]
                derivative_xi = tuple(
                    first_eta * surrogate_tangent[d]
                    + second_eta * true_tangent[d]
                    for d in range(2)
                )
                derivative_eta = (
                    0.5 * (true_point[0] - surrogate_point[0]),
                    0.5 * (true_point[1] - surrogate_point[1]),
                )
                patch_jacobian = abs(
                    derivative_xi[0] * derivative_eta[1]
                    - derivative_xi[1] * derivative_eta[0]
                )
                integration_weight = weight_x * weight_y * patch_jacobian

                parent_xi = 2.0 * (x - x_min) / cell_size_x - 1.0
                parent_eta = 2.0 * (y - y_min) / cell_size_y - 1.0
                shape_functions = (
                    0.25 * (1.0 - parent_xi) * (1.0 - parent_eta),
                    0.25 * (1.0 + parent_xi) * (1.0 - parent_eta),
                    0.25 * (1.0 + parent_xi) * (1.0 + parent_eta),
                    0.25 * (1.0 - parent_xi) * (1.0 + parent_eta),
                )
                derivatives_x = (
                    -0.5 * (1.0 - parent_eta) / cell_size_x,
                    0.5 * (1.0 - parent_eta) / cell_size_x,
                    0.5 * (1.0 + parent_eta) / cell_size_x,
                    -0.5 * (1.0 + parent_eta) / cell_size_x,
                )
                derivatives_y = (
                    -0.5 * (1.0 - parent_xi) / cell_size_y,
                    -0.5 * (1.0 + parent_xi) / cell_size_y,
                    0.5 * (1.0 + parent_xi) / cell_size_y,
                    0.5 * (1.0 - parent_xi) / cell_size_y,
                )

                velocity_x = sum(
                    shape_functions[i] * nodal_velocity_x[i] for i in range(4)
                )
                velocity_y = sum(
                    shape_functions[i] * nodal_velocity_y[i] for i in range(4)
                )
                pressure = sum(
                    shape_functions[i] * nodal_pressure[i] for i in range(4)
                )
                numerical_gradient = (
                    sum(derivatives_x[i] * nodal_velocity_x[i] for i in range(4)),
                    sum(derivatives_y[i] * nodal_velocity_x[i] for i in range(4)),
                    sum(derivatives_x[i] * nodal_velocity_y[i] for i in range(4)),
                    sum(derivatives_y[i] * nodal_velocity_y[i] for i in range(4)),
                )
                exact_u = study.exact_velocity(x, y, final_time)
                exact_p = study.exact_pressure(x, y, final_time)
                exact_gradient = study.exact_velocity_gradient(x, y, final_time)

                integrals["area"] += integration_weight
                integrals["velocity_error_squared"] += integration_weight * (
                    (velocity_x - exact_u[0]) ** 2
                    + (velocity_y - exact_u[1]) ** 2
                )
                integrals["velocity_norm_squared"] += integration_weight * (
                    exact_u[0] ** 2 + exact_u[1] ** 2
                )
                integrals["pressure_error_squared"] += integration_weight * (
                    pressure - exact_p
                ) ** 2
                integrals["pressure_error_integral"] += integration_weight * (
                    pressure - exact_p
                )
                integrals["pressure_norm_squared"] += integration_weight * exact_p**2
                integrals["h1_error_squared"] += integration_weight * sum(
                    (numerical_gradient[i] - exact_gradient[i]) ** 2
                    for i in range(4)
                )
                integrals["h1_norm_squared"] += integration_weight * sum(
                    value**2 for value in exact_gradient
                )

    return integrals


def _recover_exact_norm_squared(error: float, relative_error: float) -> float:
    if relative_error <= 1.0e-30:
        return 0.0
    return (error / relative_error) ** 2


def _calculate_full_domain_error_norms(
    fluid: KM.ModelPart, final_time: float
) -> dict[str, float]:
    """Combine retained-cell and gap-patch contributions."""
    active = _background_error_norms(fluid, final_time)
    gap = _calculate_gap_error_integrals(fluid, final_time)
    diagnostics = {}
    if _ENABLE_INTERFACE_DIAGNOSTICS:
        diagnostics = _calculate_boundary_and_interface_diagnostics(
            fluid, final_time
        )

    _BACKGROUND_ONLY_ERRORS.append(
        {
            "velocity_l2": active["velocity_l2"],
            "pressure_l2": active["pressure_l2"],
            "velocity_h1": active["velocity_h1"],
        }
    )

    area = active["area"] + gap["area"]
    velocity_error_squared = active["velocity_l2"] ** 2 + gap[
        "velocity_error_squared"
    ]
    pressure_error_squared = active["pressure_l2_with_gauge"] ** 2 + gap[
        "pressure_error_squared"
    ]
    pressure_error_integral = (
        active["pressure_mean_error"] * active["area"]
        + gap["pressure_error_integral"]
    )
    h1_error_squared = active["velocity_h1"] ** 2 + gap["h1_error_squared"]
    velocity_norm_squared = _recover_exact_norm_squared(
        active["velocity_l2"], active["relative_velocity_l2"]
    ) + gap["velocity_norm_squared"]
    pressure_norm_squared = _recover_exact_norm_squared(
        active["pressure_l2"], active["relative_pressure_l2"]
    ) + gap["pressure_norm_squared"]
    h1_norm_squared = _recover_exact_norm_squared(
        active["velocity_h1"], active["relative_velocity_h1"]
    ) + gap["h1_norm_squared"]

    pressure_mean_error = pressure_error_integral / area
    pressure_error_squared_mean_free = max(
        0.0,
        pressure_error_squared - pressure_error_integral**2 / area,
    )
    velocity_l2 = math.sqrt(velocity_error_squared)
    pressure_l2 = math.sqrt(pressure_error_squared_mean_free)
    velocity_h1 = math.sqrt(h1_error_squared)
    _FULL_DOMAIN_AREAS.append(area)
    errors = {
        "area": area,
        "velocity_l2": velocity_l2,
        "pressure_l2": pressure_l2,
        "pressure_l2_with_gauge": math.sqrt(pressure_error_squared),
        "pressure_mean_error": pressure_mean_error,
        "velocity_h1": velocity_h1,
        "relative_velocity_l2": velocity_l2 / math.sqrt(velocity_norm_squared),
        "relative_pressure_l2": pressure_l2 / math.sqrt(pressure_norm_squared),
        "relative_velocity_h1": velocity_h1 / math.sqrt(h1_norm_squared),
    }
    errors.update(diagnostics)
    return errors


def _calculate_gap_rates(rows: list[dict[str, float]]) -> None:
    """Add convergence rates for the boundary and interface diagnostics."""
    _base_calculate_rates(rows)
    if not _ENABLE_INTERFACE_DIAGNOSTICS:
        return
    diagnostic_names = (
        "boundary_trace_l2",
        "interface_velocity_jump_l2",
        "interface_traction_jump_l2",
    )
    for row in rows:
        for name in diagnostic_names:
            row[f"{name}_rate"] = math.nan
    for index in range(1, len(rows)):
        mesh_ratio = rows[index - 1]["h"] / rows[index]["h"]
        for name in diagnostic_names:
            previous_value = rows[index - 1][name]
            current_value = rows[index][name]
            if previous_value > 0.0 and current_value > 0.0:
                rows[index][f"{name}_rate"] = math.log(
                    previous_value / current_value
                ) / math.log(mesh_ratio)


def _plot_boundary_and_interface_diagnostics(rows) -> Path:
    """Plot E_Gamma, J_u and J_t against the background mesh size."""
    h = study.np.asarray([row["h"] for row in rows], dtype=float)
    quantities = (
        (
            "boundary_trace_l2",
            r"$E_\Gamma = \|u_h-g\|_{L^2(\Gamma)}$",
        ),
        (
            "interface_velocity_jump_l2",
            r"$J_u = \|[u_h]\|_{L^2(\Gamma_I)}$",
        ),
        (
            "interface_traction_jump_l2",
            r"$J_t = \|[\sigma_h n]\|_{L^2(\Gamma_I)}$",
        ),
    )
    figure, axes = study.plt.subplots(
        1, 3, figsize=(14.5, 4.5), constrained_layout=True
    )
    for axis, (name, title) in zip(axes, quantities):
        values = study.np.asarray([row[name] for row in rows], dtype=float)
        positive = values > 0.0
        if study.np.any(positive):
            axis.loglog(h[positive], values[positive], "o-")
        else:
            axis.plot(h, values, "o-")
        axis.invert_xaxis()
        axis.grid(True, which="both", alpha=0.3)
        axis.set_xlabel("background mesh size h")
        axis.set_title(title)
    axes[0].set_ylabel("integrated diagnostic")
    figure.suptitle(
        "Static Stokes Gap-SBM: boundary enforcement and gap continuity"
    )
    output = RESULTS_DIRECTORY / "boundary_and_interface_diagnostics.png"
    figure.savefig(output, dpi=220)
    if "--no-show" in sys.argv:
        study.plt.close(figure)
    return output


def _plot_error_with_and_without_gap(rows, problem_label="Static Stokes Gap-SBM"):
    """Plot full-domain errors against errors on retained cells only."""
    global _ERROR_COMPARISON_PLOT
    global _INTERFACE_DIAGNOSTIC_PLOT
    global _DIAGNOSTIC_ROWS

    convergence_output, convergence_figure = _base_plot_convergence(
        rows, problem_label
    )
    if len(_BACKGROUND_ONLY_ERRORS) != len(rows):
        raise RuntimeError(
            "Expected one background-only error entry for every mesh level."
        )

    h = study.np.asarray([row["h"] for row in rows], dtype=float)
    quantities = (
        ("velocity_l2", r"$\|u-u_h\|_{L^2}$"),
        ("pressure_l2", r"$\|p-p_h\|_{L^2}$"),
        ("velocity_h1", r"$|u-u_h|_{H^1}$"),
    )
    comparison_figure, axes = study.plt.subplots(
        1, 3, figsize=(14.5, 4.5), constrained_layout=True
    )
    for axis, (name, title) in zip(axes, quantities):
        full_domain_error = study.np.asarray(
            [row[name] for row in rows], dtype=float
        )
        retained_cells_error = study.np.asarray(
            [entry[name] for entry in _BACKGROUND_ONLY_ERRORS], dtype=float
        )
        axis.loglog(
            h,
            full_domain_error,
            "o-",
            label="with gap: full physical domain",
        )
        axis.loglog(
            h,
            retained_cells_error,
            "s--",
            label="without gap: retained cells only",
        )
        axis.invert_xaxis()
        axis.grid(True, which="both", alpha=0.3)
        axis.set_xlabel("background mesh size h")
        axis.set_title(title)
    axes[0].set_ylabel("absolute error")
    axes[0].legend()
    comparison_figure.suptitle(
        "Static Stokes Gap-SBM: influence of the gap-region error"
    )
    _ERROR_COMPARISON_PLOT = RESULTS_DIRECTORY / "error_with_and_without_gap.png"
    comparison_figure.savefig(_ERROR_COMPARISON_PLOT, dpi=220)
    if _ENABLE_INTERFACE_DIAGNOSTICS:
        _INTERFACE_DIAGNOSTIC_PLOT = (
            _plot_boundary_and_interface_diagnostics(rows)
        )
        _DIAGNOSTIC_ROWS = [dict(row) for row in rows]
    if "--no-show" in sys.argv:
        study.plt.close(comparison_figure)

    return convergence_output, convergence_figure


def _create_gap_modeler_settings(
    true_skin_name: str, number_of_cells: int
) -> KM.Parameters:
    return KM.Parameters(
        f"""
        {{
            "echo_level" : 0,
            "domain_size" : 2,
            "input_model_part_name" : "{true_skin_name}",
            "output_model_part_name" : "FluidModelPart",
            "volume_sub_model_part_name" : "FluidDomain",
            "surrogate_sub_model_part_name" : "{study.SURROGATE_SUB_MODEL_PART_NAME}",
            "gap_element_sub_model_part_name" : "GapElements",
            "gap_interface_sub_model_part_name" : "GapInterfaces",
            "lower_point" : [{study.DOMAIN_MIN[0]}, {study.DOMAIN_MIN[1]}, 0.0],
            "upper_point" : [{study.DOMAIN_MAX[0]}, {study.DOMAIN_MAX[1]}, 0.0],
            "number_of_elements" : [{number_of_cells}, {number_of_cells}, 1],
            "physical_domain" : "outside",
            "lambda" : 1.0,
            "classification_samples_per_direction" : 3,
            "element_name" : "SymbolicStokes2D4N",
            "condition_name" : "SbmFluidDirichletCondition2D4N",
            "gap_element_name" : "SymbolicStokes2D4N",
            "gap_interface_condition_name" : "SbmFluidGapInterfaceCondition2D",
            "properties_id" : 1,
            "number_of_boundary_gauss_points" : 2,
            "gap_integration_order" : 2
        }}
        """
    )


def _run_gap_modeler(
    model: KM.Model, number_of_cells: int
) -> tuple[KM.ModelPart, list[study.Segment2D], float]:
    true_skin = study._create_true_skin(model)
    true_segments = study._extract_segments(true_skin)
    fluid = (
        model.GetModelPart("FluidModelPart")
        if model.HasModelPart("FluidModelPart")
        else model.CreateModelPart("FluidModelPart")
    )
    if not fluid.HasProperties(0):
        fluid.CreateNewProperties(0)
    if not fluid.HasProperties(1):
        fluid.CreateNewProperties(1)

    modeler = KM.CreateModeler(
        "GapSurrogateBoundaryModeler",
        model,
        _create_gap_modeler_settings(true_skin.Name, number_of_cells),
    )
    start_time = time.perf_counter()
    modeler.SetupModelPart()
    elapsed_time = time.perf_counter() - start_time
    return fluid, true_segments, elapsed_time


class StaticGapSbmAnalysis(study.ManufacturedSbmAnalysis):
    """Adds gap volume and internal interfaces to the computing ModelPart."""

    def _InitializeInternals(self) -> None:
        computing_model_part = self._GetSolver().GetComputingModelPart()
        root = self.model.GetModelPart("FluidModelPart")
        gap_elements = root.GetSubModelPart("GapElements")
        gap_interfaces = root.GetSubModelPart("GapInterfaces")
        computing_model_part.AddElements(
            [element.Id for element in gap_elements.Elements]
        )
        computing_model_part.AddConditions(
            [condition.Id for condition in gap_interfaces.Conditions]
        )
        super()._InitializeInternals()


def _configure_study() -> None:
    """Configure the reusable standard-SBM driver for this Gap-SBM case."""
    study.EXAMPLE_DIRECTORY = EXAMPLE_DIRECTORY
    study.RESULTS_DIRECTORY = RESULTS_DIRECTORY
    study.SURROGATE_SUB_MODEL_PART_NAME = (
        "DiamondSbm" if GEOMETRY == "diamond" else "CylinderSbm"
    )
    study._set_manufactured_solution = _set_gap_manufactured_solution
    _set_gap_manufactured_solution("sin-cos")
    study._create_true_skin = _create_selected_true_skin
    study._run_surrogate_modeler = _run_gap_modeler
    study._calculate_error_norms = _calculate_full_domain_error_norms
    study._calculate_rates = _calculate_gap_rates
    study._plot_convergence = _plot_error_with_and_without_gap
    study.ManufacturedSbmAnalysis = StaticGapSbmAnalysis


def main() -> None:
    # Patch only the reusable driver's extension points. The original standard
    # SBM example is left untouched and remains available as a reference.
    global _ERROR_COMPARISON_PLOT
    global _INTERFACE_DIAGNOSTIC_PLOT
    _consume_interface_diagnostics_flag()
    _BACKGROUND_ONLY_ERRORS.clear()
    _FULL_DOMAIN_AREAS.clear()
    _DIAGNOSTIC_ROWS.clear()
    _ERROR_COMPARISON_PLOT = None
    _INTERFACE_DIAGNOSTIC_PLOT = None
    _configure_study()
    study.main()
    if _ERROR_COMPARISON_PLOT is not None:
        print(f"  Error with/without gap : {_ERROR_COMPARISON_PLOT}")
    if _INTERFACE_DIAGNOSTIC_PLOT is not None:
        print(f"  Boundary/interface plot: {_INTERFACE_DIAGNOSTIC_PLOT}")
    if _DIAGNOSTIC_ROWS:
        print("\nBoundary and gap-interface diagnostics")
        print(
            "  n       E_Gamma       rate    J_u           rate    "
            "J_t           rate"
        )
        for row in _DIAGNOSTIC_ROWS:
            rates = []
            for name in (
                "boundary_trace_l2_rate",
                "interface_velocity_jump_l2_rate",
                "interface_traction_jump_l2_rate",
            ):
                rates.append(
                    "  -  " if math.isnan(row[name]) else f"{row[name]:5.2f}"
                )
            print(
                f"  {int(row['n']):<5d} "
                f"{row['boundary_trace_l2']:<13.4e} {rates[0]}  "
                f"{row['interface_velocity_jump_l2']:<13.4e} {rates[1]}  "
                f"{row['interface_traction_jump_l2']:<13.4e} {rates[2]}"
            )
    if GEOMETRY == "diamond" and _FULL_DOMAIN_AREAS:
        finest_area = _FULL_DOMAIN_AREAS[-1]
        print(f"  Exact diamond-fluid area: {DIAMOND_FLUID_AREA:.12f}")
        print(f"  Finest integrated area  : {finest_area:.12f}")
        print(
            "  Finest absolute area error: "
            f"{abs(finest_area - DIAMOND_FLUID_AREA):.4e}"
        )


if __name__ == "__main__":
    main()
