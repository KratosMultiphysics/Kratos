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

case_dir = "/Users/janchauncey/Documents/Repos/Kratos/applications/PoromechanicsApplication/test_examples/Buckley_Leverett_Kratos"
os.chdir(case_dir)

class BuckleyLeverettAnalysis(PoromechanicsAnalysis):
    def __init__(self, model, parameters):
        super().__init__(model, parameters)
        day = 86400.0
        self.dt = 0.001 * day
        self.dt_max = 1.0 * day
        self.dt_min = 0.001 * day
        self.target_time = 500.0 * day

    def Initialize(self):
        super().Initialize()
        mp = self._GetSolver().GetComputingModelPart()
        for nd in mp.Nodes:
            if nd.X < 1e-5:
                nd.SetSolutionStepValue(KP.GAS_PRESSURE, 0, 0.0)
                nd.SetSolutionStepValue(KP.GAS_PRESSURE, 1, 0.0)
                nd.Fix(KP.GAS_PRESSURE)
                nd.SetSolutionStepValue(KP.LIQUID_PRESSURE, 0, 0.0)
                nd.SetSolutionStepValue(KP.LIQUID_PRESSURE, 1, 0.0)
                nd.Fix(KP.LIQUID_PRESSURE)
            else:
                nd.SetSolutionStepValue(KP.GAS_PRESSURE, 0, 1.0)
                nd.SetSolutionStepValue(KP.GAS_PRESSURE, 1, 1.0)
                nd.SetSolutionStepValue(KP.LIQUID_PRESSURE, 0, 0.0)
                nd.SetSolutionStepValue(KP.LIQUID_PRESSURE, 1, 0.0)
        print("[KRATOS BL] Correctly initialized initial conditions and Dirichlet BCs.")

    def RunSolutionLoop(self):
        step = 0
        time = self.time

        while time < self.end_time - 1e-4:
            step += 1
            if time + self.dt > self.end_time:
                dt_step = self.end_time - time
            else:
                dt_step = self.dt

            comp_mp = self._GetSolver().GetComputingModelPart()
            comp_mp.ProcessInfo[KM.DELTA_TIME] = dt_step

            time = self._GetSolver().AdvanceInTime(time)
            self.time = time

            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            converged = self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()

            if step % 25 == 0 or time >= self.end_time - 1e-4:
                days = time / 86400.0
                print(f"[KRATOS BL] Step {step}: Time = {days:.2f} days / 500 days (dt = {dt_step/86400.0:.2f} d, converged={converged})")

            if converged:
                self.dt = min(self.dt * 1.2, self.dt_max)
            else:
                self.dt = max(self.dt * 0.5, self.dt_min)

        print("[KRATOS BL] Simulation finished! Extracting results...")
        self.extract_results()

    def extract_results(self):
        mp = self._GetSolver().GetComputingModelPart()
        data = []
        for nd in mp.Nodes:
            nid = nd.Id
            x = nd.X
            y = nd.Y
            pl = nd.GetSolutionStepValue(KP.LIQUID_PRESSURE)
            pg = nd.GetSolutionStepValue(KP.GAS_PRESSURE)
            pc = pg - pl
            
            # Use PiecewiseLinearLaw formula for Sl
            # Sl = 0.8 - 0.6 * pc, clamped between 0.2 and 0.8
            if pc <= 0.0:
                sl = 0.8
            elif pc >= 1.0:
                sl = 0.2
            else:
                sl = 0.8 - 0.6 * pc
            data.append([nid, x, y, pl, pg, pc, sl])

        df = pd.DataFrame(data, columns=['NodeId', 'x', 'y', 'Pl', 'Pg', 'Pc', 'Sl'])
        df = df.sort_values(by=['x', 'y']).reset_index(drop=True)
        csv_path = os.path.join(case_dir, "kratos_bl_results.csv")
        df.to_csv(csv_path, index=False)
        print(f"[KRATOS BL] Results successfully saved to {csv_path}")

if __name__ == "__main__":
    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    simulation = BuckleyLeverettAnalysis(model, parameters)
    simulation.Run()

