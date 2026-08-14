import csv
import os

import KratosMultiphysics
import KratosMultiphysics.PoromechanicsApplication as KratosPoro


class BofangLifecycleRecorder:
    """Records nodal temperatures at lifecycle stages and per-step results.

    All files are written as CSV so that variants can be compared robustly.
    """

    def __init__(self, output_dir, variant, model, model_part_name, monitored_node_ids):
        self.output_dir = output_dir
        self.variant = variant
        self.model = model
        self.model_part_name = model_part_name
        self.monitored_node_ids = sorted(monitored_node_ids)
        os.makedirs(output_dir, exist_ok=True)

        self.lifecycle_file = os.path.join(output_dir, "lifecycle_temperatures.csv")
        self.results_file = os.path.join(output_dir, "results.csv")
        self.stress_file = os.path.join(output_dir, "stress_results.csv")

        with open(self.lifecycle_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["variant", "lifecycle_stage", "node_id", "x", "y", "temperature"])

        with open(self.results_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "variant", "time", "node_id", "temperature",
                "displacement_x", "displacement_y", "displacement_norm"
            ])

        with open(self.stress_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "variant", "time", "node_id", "stress_xx", "stress_yy", "stress_xy"
            ])

    def _temperature(self, node, step=0):
        return node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE, step)

    def _displacement(self, node, step=0):
        disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, step)
        return disp[0], disp[1]

    def _nodal_stress(self, node, step=0):
        stress = node.GetSolutionStepValue(KratosPoro.NODAL_CAUCHY_STRESS_TENSOR, step)
        if stress.Size1() < 2 or stress.Size2() < 2:
            return 0.0, 0.0, 0.0
        return stress[0, 0], stress[1, 1], stress[0, 1]

    def record_lifecycle_stage(self, stage):
        model_part = self.model[self.model_part_name]
        with open(self.lifecycle_file, "a", newline="") as f:
            writer = csv.writer(f)
            for nid in self.monitored_node_ids:
                node = model_part.GetNode(nid)
                writer.writerow([
                    self.variant, stage, nid, node.X, node.Y,
                    self._temperature(node)
                ])

    def record_step(self):
        model_part = self.model[self.model_part_name]
        time = model_part.ProcessInfo[KratosMultiphysics.TIME]
        with open(self.results_file, "a", newline="") as f:
            writer = csv.writer(f)
            for nid in self.monitored_node_ids:
                node = model_part.GetNode(nid)
                ux, uy = self._displacement(node)
                disp_norm = (ux * ux + uy * uy) ** 0.5
                writer.writerow([
                    self.variant, time, nid, self._temperature(node), ux, uy, disp_norm
                ])

        with open(self.stress_file, "a", newline="") as f:
            writer = csv.writer(f)
            for nid in self.monitored_node_ids:
                node = model_part.GetNode(nid)
                sxx, syy, sxy = self._nodal_stress(node)
                writer.writerow([self.variant, time, nid, sxx, syy, sxy])
