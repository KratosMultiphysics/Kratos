"""Test-only probe process for the lightweight initial-state workflow.

Records the model-part state at lifecycle boundaries without modifying it. It is
used to prove what the destination contains before the first solve and after
each completed solution step during a real ``DamAnalysis.Run()``.
"""
import json

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as Dam


def Factory(settings, model):
    return StateProbeProcess(model, settings["Parameters"])


class StateProbeProcess(KM.Process):
    def __init__(self, model, settings):
        super().__init__()
        self.model_part = model[settings["model_part_name"].GetString()]
        self.output_file_name = settings["output_file_name"].GetString()
        self.records = {}

    def ExecuteBeforeSolutionLoop(self):
        self.records["initial"] = self._Snapshot()

    def ExecuteFinalizeSolutionStep(self):
        self.records.setdefault("steps", []).append(self._Snapshot())

    def ExecuteFinalize(self):
        with open(self.output_file_name, "w", encoding="utf-8") as stream:
            json.dump(self.records, stream)

    def _Snapshot(self):
        part = self.model_part
        record = {
            "time": float(part.ProcessInfo[KM.TIME]),
            "step": int(part.ProcessInfo[KM.STEP]),
            "nodes": [node.Id for node in part.Nodes],
            "displacement": [list(node.GetSolutionStepValue(KM.DISPLACEMENT))
                             for node in part.Nodes],
        }
        if part.HasNodalSolutionStepVariable(KM.REACTION):
            record["reaction"] = [list(node.GetSolutionStepValue(KM.REACTION))
                                  for node in part.Nodes]
        if part.HasNodalSolutionStepVariable(KM.TEMPERATURE):
            record["temperature"] = [node.GetSolutionStepValue(KM.TEMPERATURE)
                                     for node in part.Nodes]
            record["temperature_buffer"] = [
                [node.GetSolutionStepValue(KM.TEMPERATURE, step)
                 for step in range(node.GetBufferSize())]
                for node in part.Nodes]
        if part.HasNodalSolutionStepVariable(Dam.NODAL_REFERENCE_TEMPERATURE):
            record["reference_temperature"] = [
                node.GetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE)
                for node in part.Nodes]
            record["reference_temperature_buffer"] = [
                [node.GetSolutionStepValue(Dam.NODAL_REFERENCE_TEMPERATURE, step)
                 for step in range(node.GetBufferSize())]
                for node in part.Nodes]
        if part.HasNodalSolutionStepVariable(KM.VOLUME_ACCELERATION):
            record["volume_acceleration"] = [
                list(node.GetSolutionStepValue(KM.VOLUME_ACCELERATION))
                for node in part.Nodes]
        record["cauchy_stress"] = {}
        record["thermal_strain"] = {}
        residual = {}
        thermal_strain_variable = KM.KratosGlobals.GetVariable("THERMAL_STRAIN_VECTOR")
        for element in part.Elements:
            values = element.CalculateOnIntegrationPoints(
                KM.CAUCHY_STRESS_VECTOR, part.ProcessInfo)
            record["cauchy_stress"][str(element.Id)] = [list(value) for value in values]
            if part.HasNodalSolutionStepVariable(Dam.NODAL_REFERENCE_TEMPERATURE):
                record["thermal_strain"][str(element.Id)] = [
                    list(value) for value in element.CalculateOnIntegrationPoints(
                        thermal_strain_variable, part.ProcessInfo)]
            lhs, rhs = KM.Matrix(), KM.Vector()
            element.CalculateLocalSystem(lhs, rhs, part.ProcessInfo)
            for local_index, node in enumerate(element.GetGeometry()):
                for component in range(2):
                    key = (node.Id, component)
                    residual[key] = residual.get(key, 0.0) + rhs[2 * local_index + component]
        node_map = {node.Id: node for node in part.Nodes}
        free_residual = [
            abs(value) for (node_id, component), value in residual.items()
            if not node_map[node_id].IsFixed(
                KM.DISPLACEMENT_X if component == 0 else KM.DISPLACEMENT_Y)]
        record["free_residual_max"] = (max(free_residual) if free_residual else 0.0)
        return record
