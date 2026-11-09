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

case_dir = "/Users/janchauncey/Documents/Repos/Kratos/applications/PoromechanicsApplication/test_examples/DNAPL_Infiltration_2D_Kratos"
os.chdir(case_dir)

class DNAPLAnalysis(PoromechanicsAnalysis):
    def __init__(self, model, parameters):
        super().__init__(model, parameters)
        self.dt = 1.0
        self.dt_max = 20.0
        self.dt_min = 0.01

    def Initialize(self):
        super().Initialize()
        mp = self._GetSolver().GetComputingModelPart()
        Lx = 0.90
        tol = 1.0e-5

        # Initialize hydrostatic profile and fix lateral boundaries exactly as PorousLab
        for nd in mp.Nodes:
            y = nd.Y
            pg = 7886.5 - 9810.0 * y
            pl = pg - 755.0

            nd.SetSolutionStepValue(KP.GAS_PRESSURE, 0, pg)
            nd.SetSolutionStepValue(KP.GAS_PRESSURE, 1, pg)
            nd.SetSolutionStepValue(KP.LIQUID_PRESSURE, 0, pl)
            nd.SetSolutionStepValue(KP.LIQUID_PRESSURE, 1, pl)

            if abs(nd.X) < tol or abs(nd.X - Lx) < tol:
                nd.Fix(KP.LIQUID_PRESSURE)
                nd.Fix(KP.GAS_PRESSURE)

        print("[KRATOS DNAPL] Initialized hydrostatic profile and fixed lateral boundaries.")

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

            if step % 20 == 0 or time >= self.end_time - 1e-4:
                print(f"[KRATOS DNAPL] Step {step}: Time = {time:.2f} s / 2000.0 s (dt = {dt_step:.2f} s, converged={converged})")

            if converged:
                self.dt = min(self.dt * 1.25, self.dt_max)
            else:
                self.dt = max(self.dt * 0.5, self.dt_min)

        print("[KRATOS DNAPL] Simulation finished! Extracting results...")
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
            
            # Brooks and Corey Law:
            # Pb = 755, lambda = 2.7, Slr = 0.1, Sgr = 0.0
            # If Pc <= Pb: Sl = 1.0, Sg = 0.0
            # Else: Se = (Pb / Pc)^lambda, Sl = Slr + Se * (1 - Slr), Sg = 1 - Sl
            if pc <= 755.0:
                sl = 1.0
                sg = 0.0
            else:
                se = (755.0 / pc) ** 2.7
                sl = 0.1 + se * 0.9
                sg = 1.0 - sl

            data.append([nid, x, y, pl, pg, pc, sl, sg])

        df = pd.DataFrame(data, columns=['NodeId', 'x', 'y', 'Pl', 'Pg', 'Pc', 'Sl', 'Sg'])
        df = df.sort_values(by=['x', 'y']).reset_index(drop=True)
        csv_path = os.path.join(case_dir, "kratos_dnapl_results.csv")
        df.to_csv(csv_path, index=False)
        print(f"[KRATOS DNAPL] Results successfully saved to {csv_path}")

if __name__ == "__main__":
    with open("ProjectParameters.json", 'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    simulation = DNAPLAnalysis(model, parameters)
    simulation.Run()

