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
