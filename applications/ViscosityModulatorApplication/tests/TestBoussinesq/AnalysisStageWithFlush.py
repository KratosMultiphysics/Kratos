import os
import sys
import time

import numpy as np

import scripts.processes as prcs
import scripts.utils as utils


def CreateAnalysisStageWithFlushInstance(cls, global_model, parameters, config):
    class AnalysisStageWithFlush(cls):

        def __init__(self, model, project_parameters, flush_frequency=10.0):
            super().__init__(model, project_parameters)
            self.flush_frequency = flush_frequency
            self.last_flush = time.time()
            sys.stdout.flush()

            self.config = config

        def get_inlet_model_part(self):
            return self.model["FluidModelPart"].GetSubModelPart("AutomaticInlet3D_Automatic_inlet_velocity_Auto1")

        def Initialize(self):
            super().Initialize()
            sys.stdout.flush()

            prcs.impose_inlet_solution(self.get_inlet_model_part(), self.config["zi_layers"], self.config["ui_layers"])

            # Compute lumped nodal areas (volumes in 3D) for L2 norm calculation
            import KratosMultiphysics as KM
            model_part = self.model["FluidModelPart"]
            domain_size = model_part.ProcessInfo[KM.DOMAIN_SIZE]
            nodal_area_process = KM.CalculateNodalAreaProcess(model_part, domain_size)
            nodal_area_process.Execute()

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()

            if self.parallel_type == "OpenMP":
                now = time.time()
                if now - self.last_flush > self.flush_frequency:
                    sys.stdout.flush()
                    self.last_flush = now

        def Finalize(self):
            super().Finalize()
            if self.config["export_results"]:
                self._export_solution_data()

        def _export_solution_data(self):
            """Export mesh + nodal solution + nodal volumes to HDF5 and MDPA files."""
            import h5py
            import KratosMultiphysics as KM

            h = self.config["h"]
            output_dir = os.path.join(utils.get_case_dir(), "convergence_results")
            os.makedirs(output_dir, exist_ok=True)
            h_str = str(h).replace(".", "_")

            model_part = self.model["FluidModelPart"]

            # Collect node data
            node_ids = []
            coords = []
            temperatures = []
            velocities = []
            pressures = []
            nodal_areas = []

            for node in model_part.Nodes:
                node_ids.append(node.Id)
                coords.append([node.X, node.Y, node.Z])
                temperatures.append(node.GetSolutionStepValue(KM.TEMPERATURE))
                vel = node.GetSolutionStepValue(KM.VELOCITY)
                velocities.append([vel[0], vel[1], vel[2]])
                pressures.append(node.GetSolutionStepValue(KM.PRESSURE))
                nodal_areas.append(node.GetSolutionStepValue(KM.NODAL_AREA))

            # Collect element connectivity
            elem_ids = []
            connectivity = []
            for elem in model_part.Elements:
                elem_ids.append(elem.Id)
                connectivity.append([n.Id for n in elem.GetNodes()])

            # Write HDF5 file
            h5_path = os.path.join(output_dir, f"solution_h{h_str}.h5")
            with h5py.File(h5_path, 'w') as f:
                f.create_dataset("nodes/ids", data=np.array(node_ids, dtype=np.int64))
                f.create_dataset("nodes/coordinates", data=np.array(coords, dtype=np.float64))
                f.create_dataset("elements/ids", data=np.array(elem_ids, dtype=np.int64))
                f.create_dataset("elements/connectivity", data=np.array(connectivity, dtype=np.int64))
                f.create_dataset("solution/TEMPERATURE", data=np.array(temperatures, dtype=np.float64))
                f.create_dataset("solution/VELOCITY", data=np.array(velocities, dtype=np.float64))
                f.create_dataset("solution/PRESSURE", data=np.array(pressures, dtype=np.float64))
                f.create_dataset("nodal_areas", data=np.array(nodal_areas, dtype=np.float64))
                f.attrs["h"] = h
                f.attrs["num_nodes"] = len(node_ids)
                f.attrs["num_elements"] = len(elem_ids)

            print(f"[ConvergenceExport] Saved solution for h={h} to {h5_path}")
            print(f"  Nodes: {len(node_ids)}, Elements: {len(elem_ids)}")

            # Also write a results MDPA file with nodal solution data
            mdpa_path = os.path.join(output_dir, f"solution_h{h_str}.mdpa")
            self._write_solution_mdpa(mdpa_path, node_ids, coords,
                                      temperatures, velocities, pressures,
                                      nodal_areas, elem_ids, connectivity)

        def _write_solution_mdpa(self, filepath, node_ids, coords,
                                 temperatures, velocities, pressures,
                                 nodal_areas, elem_ids, connectivity):
            """Write an MDPA file containing mesh and nodal solution data."""
            with open(filepath, 'w') as f:
                f.write("Begin ModelPartData\nEnd ModelPartData\n\n")
                f.write("Begin Properties 0\nEnd Properties\n\n")

                # Nodes
                f.write("Begin Nodes\n")
                for i, nid in enumerate(node_ids):
                    f.write(f"    {nid}    {coords[i][0]:.10e}    {coords[i][1]:.10e}    {coords[i][2]:.10e}\n")
                f.write("End Nodes\n\n")

                # Elements (tetrahedra)
                if connectivity:
                    f.write("Begin Elements Element3D4N\n")
                    for i, eid in enumerate(elem_ids):
                        nodes_str = "  ".join(str(n) for n in connectivity[i])
                        f.write(f"    {eid}     0  {nodes_str}\n")
                    f.write("End Elements\n\n")

                # Nodal solution data: TEMPERATURE
                f.write("Begin NodalData TEMPERATURE\n")
                for i, nid in enumerate(node_ids):
                    f.write(f"    {nid}  0  {temperatures[i]:.10e}\n")
                f.write("End NodalData\n\n")

                # Nodal solution data: VELOCITY
                f.write("Begin NodalData VELOCITY_X\n")
                for i, nid in enumerate(node_ids):
                    f.write(f"    {nid}  0  {velocities[i][0]:.10e}\n")
                f.write("End NodalData\n\n")

                f.write("Begin NodalData VELOCITY_Y\n")
                for i, nid in enumerate(node_ids):
                    f.write(f"    {nid}  0  {velocities[i][1]:.10e}\n")
                f.write("End NodalData\n\n")

                f.write("Begin NodalData VELOCITY_Z\n")
                for i, nid in enumerate(node_ids):
                    f.write(f"    {nid}  0  {velocities[i][2]:.10e}\n")
                f.write("End NodalData\n\n")

                # Nodal solution data: PRESSURE
                f.write("Begin NodalData PRESSURE\n")
                for i, nid in enumerate(node_ids):
                    f.write(f"    {nid}  0  {pressures[i]:.10e}\n")
                f.write("End NodalData\n\n")

                # Nodal data: NODAL_AREA (lumped volume)
                f.write("Begin NodalData NODAL_AREA\n")
                for i, nid in enumerate(node_ids):
                    f.write(f"    {nid}  0  {nodal_areas[i]:.10e}\n")
                f.write("End NodalData\n\n")

            print(f"[ConvergenceExport] Saved results MDPA to {filepath}")

    return AnalysisStageWithFlush(global_model, parameters)