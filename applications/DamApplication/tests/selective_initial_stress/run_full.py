"""Child harness that executes a complete, unmodified DamAnalysis.Run().

Used by the lightweight initial-state end-to-end test so that no solver or
process callback is monkey-patched.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as Dam

Dam.__path__.insert(0, sys.argv[1])
from KratosMultiphysics.DamApplication.dam_analysis import DamAnalysis

with open(sys.argv[2], "r", encoding="utf-8") as stream:
    parameters = KM.Parameters(stream.read())

simulation = DamAnalysis(KM.Model(), parameters)
simulation.Run()
