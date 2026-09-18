#  ____             _       |
# |  _ \ _   _  ___| | __   |
# | |_) | | | |/ __| |/ /   |
# |  __/| |_| | (__|   <    |
# |_|    \__,_|\___|_|\_\   ANALYSIS
#
#  Main authors:    Lucas Rimpl
#  Co-authors:      Tobias Siemer
#
import numpy as np
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.analysis_result import AnalysisResult
from KratosMultiphysics.StructuralMechanicsApplication.handbook_methods.method_base import HandbookMethod
import KratosMultiphysics.StructuralMechanicsApplication as SMA

class PuckAnalysis:

    def IsApplicable(self):
        return True

    def Evaluate(self, structural_component):
        self.PuckAnalysis()
        self.PuckFF()
