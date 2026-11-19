import json
import math
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo


def _Magnitude(components):
    return math.sqrt(sum(component * component for component in components))


class IncrementalCPhiReductionDiagnostics:
    def __init__(self, output_path="incremental_c_phi_reduction_results.json",
                 monitored_node_ids=()):
        self.output_path = Path(output_path)
        self.monitored_node_ids = list(monitored_node_ids)
        self.trial_results = []
        self.factor_displacement_points = []
        self.monitored_displacement_points = []
        # DISPLACEMENT of each monitored node at the last converged trial
        self.last_converged_node_displacements = {
            node_id: [0.0, 0.0, 0.0] for node_id in self.monitored_node_ids}

    def RecordTrial(self, model_part, factor, converged, elapsed_time):
        process_info = model_part.ProcessInfo
        iterations = process_info[Kratos.NL_ITERATION_NUMBER]
        residual = (process_info[Kratos.RESIDUAL_NORM]
                    if process_info.Has(Kratos.RESIDUAL_NORM) else None)
        max_displacement, node_id, components = self._GetMaximumDisplacement(model_part)
        monitored_nodes = self._GetMonitoredNodeDisplacements(model_part)

        has_invalid_values = any(
            value is None or not math.isfinite(value)
            for value in (factor, max_displacement, elapsed_time))
        has_invalid_values = has_invalid_values or (
            residual is not None and not math.isfinite(residual))
        has_invalid_values = has_invalid_values or any(
            not math.isfinite(component) for component in components)
        has_invalid_values = has_invalid_values or any(
            not math.isfinite(component)
            for monitored_node in monitored_nodes
            for key in ("displacement_components", "total_displacement_components")
            for component in monitored_node[key])

        slope, slope_growth_indicator = None, None
        if converged and not has_invalid_values:
            self.factor_displacement_points.append((factor, max_displacement))
            point = {}
            for monitored_node in monitored_nodes:
                prefix = f"node_{monitored_node['node_id']}"
                point[f"{prefix}_displacement"] = monitored_node["displacement"]
                point[f"{prefix}_total_displacement"] = monitored_node["total_displacement"]
                point[f"{prefix}_incremental_displacement"] = monitored_node["incremental_displacement"]
                self.last_converged_node_displacements[monitored_node["node_id"]] = (
                    monitored_node["displacement_components"])
            self.monitored_displacement_points.append(point)
            slope, slope_growth_indicator = self._CalculateSlopeIndicators()

        result = {
            "trial": len(self.trial_results) + 1,
            "factor": factor,
            "converged": converged,
            "iterations": iterations,
            "max_displacement": max_displacement,
            "max_displacement_node_id": node_id,
            "displacement_components": list(components),
            "residual": residual,
            "elapsed_time": elapsed_time,
            "slope": slope,
            "slope_growth_indicator": slope_growth_indicator
        }
        if self.monitored_node_ids:
            result["monitored_nodes"] = monitored_nodes
        self.trial_results.append(result)
        self._WriteResults(model_part)

        Kratos.Logger.PrintInfo(
            "Incremental c-phi trial",
            f"F={factor}, converged={converged}, iterations={iterations}, "
            f"u_max={max_displacement} at node {node_id}, components={components}, "
            + "".join(f"u_node_{n['node_id']}={n['displacement']}, " for n in monitored_nodes) +
            f"residual={residual}, elapsed_time={elapsed_time:.6f} s, slope={slope}, "
            f"slope_growth_indicator={slope_growth_indicator}")
        return has_invalid_values

    def _CalculateSlopeIndicators(self):
        slope = None
        slope_growth_indicator = None
        if len(self.factor_displacement_points) >= 2:
            factor_0, displacement_0 = self.factor_displacement_points[-2]
            factor_1, displacement_1 = self.factor_displacement_points[-1]
            factor_difference = factor_1 - factor_0
            if factor_difference != 0.0:
                slope = (displacement_1 - displacement_0) / factor_difference

        if len(self.factor_displacement_points) >= 3 and slope is not None:
            factor_0, displacement_0 = self.factor_displacement_points[-3]
            factor_1, displacement_1 = self.factor_displacement_points[-2]
            factor_difference = factor_1 - factor_0
            if factor_difference != 0.0:
                previous_slope = (displacement_1 - displacement_0) / factor_difference
                if previous_slope != 0.0:
                    slope_growth_indicator = abs(slope) / abs(previous_slope)

        return slope, slope_growth_indicator

    @staticmethod
    def _GetMaximumDisplacement(model_part):
        max_displacement = -1.0
        max_node_id = None
        max_components = (0.0, 0.0, 0.0)

        for node in model_part.Nodes:
            displacement = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            components = tuple(float(displacement[i]) for i in range(3))
            magnitude = math.sqrt(sum(component * component for component in components))
            if magnitude > max_displacement:
                max_displacement = magnitude
                max_node_id = node.Id
                max_components = components

        if max_node_id is None:
            return None, None, max_components
        return max_displacement, max_node_id, max_components

    def _GetMonitoredNodeDisplacements(self, model_part):
        monitored_nodes = []
        for node_id in self.monitored_node_ids:
            if not model_part.HasNode(node_id):
                raise RuntimeError(f"Monitored node {node_id} is not part of {model_part.Name}")
            node = model_part.GetNode(node_id)
            displacement = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            total_displacement = node.GetSolutionStepValue(KratosGeo.TOTAL_DISPLACEMENT)
            components = [float(displacement[i]) for i in range(3)]
            total_components = [float(total_displacement[i]) for i in range(3)]
            previous_components = self.last_converged_node_displacements[node_id]
            incremental_components = [
                current - previous for current, previous in zip(components, previous_components)]
            monitored_nodes.append({
                "node_id": node_id,
                "displacement": _Magnitude(components),
                "displacement_components": components,
                "total_displacement": _Magnitude(total_components),
                "total_displacement_components": total_components,
                "incremental_displacement": _Magnitude(incremental_components),
                "incremental_displacement_components": incremental_components
            })
        return monitored_nodes

    def _WriteResults(self, model_part):
        data_communicator = model_part.GetCommunicator().GetDataCommunicator()
        if data_communicator.Rank() != 0:
            return

        output = {
            "trials": self.trial_results,
            "factor_displacement_points": [
                {"factor": factor, "max_displacement": displacement, **monitored}
                for (factor, displacement), monitored in zip(
                    self.factor_displacement_points, self.monitored_displacement_points)
            ]
        }
        if self.monitored_node_ids:
            output["monitored_node_ids"] = self.monitored_node_ids
        with open(self.output_path, "w") as output_file:
            json.dump(output, output_file, indent=2, allow_nan=False)
