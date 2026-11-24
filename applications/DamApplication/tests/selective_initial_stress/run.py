"""Child harness for the Dam initial-state workflow tests.

The child process isolates serializer failures and uses the DamApplication
scripts from the source tree under test.
"""
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as Dam
Dam.__path__.insert(0, sys.argv[1])
from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis
from test_dam_initial_state_workflow import snapshot

with open("parameters.json") as stream:
    parameters = KM.Parameters(stream.read())
simulation = DamAnalysis(KM.Model(), parameters)
simulation.Initialize()
part = simulation.main_model_part.GetSubModelPart("Parts_Concrete")
initial = snapshot(part)
steps = []
solve = simulation.solver.Solve


def recording_solve():
    solve()
    steps.append(snapshot(part))


simulation.solver.Solve = recording_solve
simulation.RunMainTemporalLoop()
simulation.Finalize()
with open("snapshots.json", "w") as stream:
    json.dump({"initial": initial, "steps": steps}, stream)
