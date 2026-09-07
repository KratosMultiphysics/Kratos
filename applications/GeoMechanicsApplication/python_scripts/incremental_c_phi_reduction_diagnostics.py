import json
import math
from pathlib import Path

import KratosMultiphysics as Kratos


class IncrementalCPhiReductionDiagnostics:
    def __init__(self, output_path="incremental_c_phi_reduction_results.json"):
        self.output_path = Path(output_path)
        self.trial_results = []
        self.factor_displacement_points = []

    def RecordTrial(self, model_part, factor, converged, elapsed_time):
        process_info = model_part.ProcessInfo
        iterations = process_info[Kratos.NL_ITERATION_NUMBER]
        residual = (process_info[Kratos.RESIDUAL_NORM]
                    if process_info.Has(Kratos.RESIDUAL_NORM) else None)
        max_displacement, node_id, components = self._GetMaximumDisplacement(model_part)

        has_invalid_values = any(
            value is None or not math.isfinite(value)
            for value in (factor, max_displacement, elapsed_time))
        has_invalid_values = has_invalid_values or (
            residual is not None and not math.isfinite(residual))
        has_invalid_values = has_invalid_values or any(
            not math.isfinite(component) for component in components)

        slope, slope_growth_indicator = None, None
        if converged and not has_invalid_values:
            self.factor_displacement_points.append((factor, max_displacement))
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
        self.trial_results.append(result)
        self._WriteResults(model_part)

        Kratos.Logger.PrintInfo(
            "Incremental c-phi trial",
            f"F={factor}, converged={converged}, iterations={iterations}, "
            f"u_max={max_displacement} at node {node_id}, components={components}, "
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

    def _WriteResults(self, model_part):
        data_communicator = model_part.GetCommunicator().GetDataCommunicator()
        if data_communicator.Rank() != 0:
            return

        output = {
            "trials": self.trial_results,
            "factor_displacement_points": [
                {"factor": factor, "max_displacement": displacement}
                for factor, displacement in self.factor_displacement_points
            ]
        }
        with open(self.output_path, "w") as output_file:
            json.dump(output, output_file, indent=2, allow_nan=False)
