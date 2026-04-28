from pathlib import Path
from typing import Dict, List, Optional, Union

import test_helper


NodalValue = Union[float, List[float]]
NodalSeries = Dict[float, Dict[int, NodalValue]]


def _parse_gid_ascii_nodal_variable(filepath: Path, variable_name: str) -> NodalSeries:
    result = test_helper.get_nodal_variable_from_ascii(str(filepath), variable_name)

    if not result:
        raise RuntimeError(f"Variable '{variable_name}' not found in {filepath}")

    return result


def _get_tolerance(lhs_value: float, rhs_value: float, abs_tol: float, rel_tol: float) -> float:
    return max(abs_tol, rel_tol * max(abs(lhs_value), abs(rhs_value), 1.0))


def _node_key_mismatch_error(lhs_nodes: Dict[int, NodalValue], rhs_nodes: Dict[int, NodalValue], time: float) -> List[str]:
    lhs_keys = set(lhs_nodes.keys())
    rhs_keys = set(rhs_nodes.keys())
    if lhs_keys == rhs_keys:
        return []

    only_lhs = sorted(lhs_keys - rhs_keys)
    only_rhs = sorted(rhs_keys - lhs_keys)
    return [
        f"Node id mismatch at time {time}: only in lhs={only_lhs[:10]}, "
        f"only in rhs={only_rhs[:10]}"
    ]


def _compare_scalar_values(
    node_id: int,
    lhs_value: float,
    rhs_value: float,
    abs_tol: float,
    rel_tol: float,
) -> List[str]:
    diff = abs(lhs_value - rhs_value)
    tol = _get_tolerance(lhs_value, rhs_value, abs_tol, rel_tol)
    if diff > tol:
        return [
            f"Node {node_id}: lhs={lhs_value}, rhs={rhs_value}, diff={diff}, tol={tol}"
        ]
    return []


def _compare_vector_values(
    node_id: int,
    lhs_value: List[float],
    rhs_value: List[float],
    abs_tol: float,
    rel_tol: float,
) -> List[str]:
    if len(lhs_value) != len(rhs_value):
        return [
            f"Component count mismatch at node {node_id}: lhs={lhs_value}, rhs={rhs_value}"
        ]

    errors: List[str] = []
    for comp, (lv, rv) in enumerate(zip(lhs_value, rhs_value)):
        diff = abs(lv - rv)
        tol = _get_tolerance(lv, rv, abs_tol, rel_tol)
        if diff > tol:
            errors.append(
                f"Node {node_id}, comp {comp}: lhs={lv}, rhs={rv}, diff={diff}, tol={tol}"
            )
    return errors


def _compare_node_values(
    node_id: int,
    lhs_value: NodalValue,
    rhs_value: NodalValue,
    abs_tol: float,
    rel_tol: float,
) -> List[str]:
    if isinstance(lhs_value, list) != isinstance(rhs_value, list):
        return [
            f"Type mismatch at node {node_id}: lhs={lhs_value}, rhs={rhs_value}"
        ]

    if isinstance(lhs_value, list):
        return _compare_vector_values(
            node_id, lhs_value, rhs_value, abs_tol, rel_tol
        )

    return _compare_scalar_values(
        node_id, lhs_value, rhs_value, abs_tol, rel_tol
    )


def _compare_nodal_series_at_time(
    lhs: NodalSeries,
    rhs: NodalSeries,
    time: float,
    abs_tol: float,
    rel_tol: float,
) -> List[str]:
    lhs_nodes = lhs[time]
    rhs_nodes = rhs[time]

    key_mismatch_error = _node_key_mismatch_error(lhs_nodes, rhs_nodes, time)
    if key_mismatch_error:
        return key_mismatch_error

    errors: List[str] = []
    for node_id in sorted(lhs_nodes.keys()):
        lhs_value = lhs_nodes[node_id]
        rhs_value = rhs_nodes[node_id]
        errors.extend(
            _compare_node_values(
                node_id, lhs_value, rhs_value, abs_tol, rel_tol
            )
        )

    return errors


def _compare_case_outputs(
    full_run_stage3: Path,
    checkpoint_stage3: Path,
    abs_tol: float,
    rel_tol: float,
    variables: Optional[List[str]] = None,
) -> None:
    if variables is None:
        variables = ["TOTAL_DISPLACEMENT", "WATER_PRESSURE"]

    failures: List[str] = []

    for variable in variables:
        full_run_series = _parse_gid_ascii_nodal_variable(
            full_run_stage3, variable
        )
        checkpoint_series = _parse_gid_ascii_nodal_variable(
            checkpoint_stage3, variable
        )
        # compare at the final time present in the full-run series
        compare_time = max(full_run_series.keys())
        errors = _compare_nodal_series_at_time(
            full_run_series, checkpoint_series, compare_time, abs_tol, rel_tol
        )

        if errors:
            failures.append(
                f"Variable {variable} mismatch at time {compare_time}. First 20 differences:\n"
                + "\n".join(errors[:20])
            )

    if failures:
        raise AssertionError("\n\n".join(failures))
