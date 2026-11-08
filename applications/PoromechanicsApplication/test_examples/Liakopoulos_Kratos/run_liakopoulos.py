import os
import sys
import numpy as np
import pandas as pd

# Path configuration for Kratos
kratos_path = "/Users/janchauncey/Documents/Repos/Kratos/bin/Release"
libs_path = "/Users/janchauncey/Documents/Repos/Kratos/bin/Release/libs"
if kratos_path not in sys.path:
    sys.path.insert(0, kratos_path)
if "DYLD_LIBRARY_PATH" not in os.environ or libs_path not in os.environ["DYLD_LIBRARY_PATH"]:
    os.environ["DYLD_LIBRARY_PATH"] = libs_path + ":" + os.environ.get("DYLD_LIBRARY_PATH", "")

import KratosMultiphysics as KM
import KratosMultiphysics.PoromechanicsApplication as KP
from KratosMultiphysics.PoromechanicsApplication.poromechanics_analysis import PoromechanicsAnalysis

case_dir = "/Users/janchauncey/Documents/Repos/Kratos/applications/PoromechanicsApplication/test_examples/Liakopoulos_Kratos"
scratch_dir = "/Users/janchauncey/.gemini/antigravity-ide/brain/39aa9aee-cf80-43bf-8ecb-3d0d2a05e231/scratch"
os.chdir(case_dir)

class LiakopoulosAnalysis(PoromechanicsAnalysis):
    def __init__(self, model, parameters):
        super().__init__(model, parameters)
        self.target_times = [120.0, 300.0, 4800.0, 7200.0]
        self.snapshots = {}
        self.dt = 0.001
        self.dt_max = 30.0

    def Initialize(self):
        super().Initialize()
        mp = self._GetSolver().GetComputingModelPart()
        stress_init = KM.Matrix(2, 2)
        stress_init[0, 0] = 101325.0
        stress_init[1, 1] = 101325.0
        for nd in mp.Nodes:
            nd.SetSolutionStepValue(KP.INITIAL_STRESS_TENSOR, stress_init)

    def RunSolutionLoop(self):
        """Custom solution loop with adaptive time stepping hitting target snapshots exactly."""
        target_idx = 0
        step = 0
        time = self.time

        while time < self.end_time - 1e-8:
            step += 1
            # Adjust dt to hit target snapshot exactly
            if target_idx < len(self.target_times):
                next_target = self.target_times[target_idx]
                if time < next_target and (time + self.dt > next_target):
                    dt_step = next_target - time
                else:
                    dt_step = self.dt
            else:
                dt_step = self.dt

            # Set delta time
            comp_mp = self._GetSolver().GetComputingModelPart()
            comp_mp.ProcessInfo[KM.DELTA_TIME] = dt_step

            # Advance in time
            time = self._GetSolver().AdvanceInTime(time)
            self.time = time

            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            converged = self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()

            # Check if target snapshot was reached
            if target_idx < len(self.target_times) and abs(time - self.target_times[target_idx]) < 1e-6:
                tt = self.target_times[target_idx]
                print(f"[KRATOS] Captured snapshot at t = {tt:.1f} s (step {step})")
                self._extract_snapshot(round(tt))
                target_idx += 1

            # Grow time step adaptively if converged
            if converged:
                self.dt = min(self.dt * 1.25, self.dt_max)
            else:
                self.dt = max(self.dt * 0.5, 0.0001)

    def _extract_snapshot(self, t_val):
        mp = self._GetSolver().GetComputingModelPart()
        nodes = list(mp.Nodes)
        # Group nodes by y coordinate (2 nodes per height, x=0 and x=0.1)
        # Sort by y coordinate
        nodes.sort(key=lambda nd: nd.Y)
        
        y_dict = {}
        for nd in nodes:
            y = round(nd.Y, 6)
            if y not in y_dict:
                y_dict[y] = []
            y_dict[y].append(nd)

        y_list = sorted(y_dict.keys())
        y_vals = []
        uy_vals = []
        pl_vals = []
        pg_vals = []
        pc_vals = []
        sl_vals = []

        for y in y_list:
            nd_pair = y_dict[y]
            uy = np.mean([nd.GetSolutionStepValue(KM.DISPLACEMENT_Y) for nd in nd_pair])
            pl = np.mean([nd.GetSolutionStepValue(KP.LIQUID_PRESSURE) for nd in nd_pair])
            pg = np.mean([nd.GetSolutionStepValue(KP.GAS_PRESSURE) for nd in nd_pair])
            pc = pg - pl

            # Liakopoulos saturation:
            if pc <= 0.0:
                sl = 1.0
            else:
                sl = max(0.2, 1.0 - 1.9722e-11 * (pc ** 2.4279))

            y_vals.append(y)
            uy_vals.append(uy)
            pl_vals.append(pl)
            pg_vals.append(pg)
            pc_vals.append(pc)
            sl_vals.append(sl)

        df = pd.DataFrame({
            'y': y_vals,
            'uy': uy_vals,
            'pl': pl_vals,
            'pg': pg_vals,
            'pc': pc_vals,
            'Sl': sl_vals
        })
        
        csv_path = os.path.join(scratch_dir, f"kratos_snapshot_{t_val}s.csv")
        df.to_csv(csv_path, index=False)
        print(f"[KRATOS] Saved {csv_path}")

def main():
    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    analysis = LiakopoulosAnalysis(model, parameters)
    analysis.Run()

if __name__ == "__main__":
    main()
