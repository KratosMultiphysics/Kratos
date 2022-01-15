import os
import sys
import math
import time as timer

import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as DEM
from   KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEM_analysis_stage

BaseAnalysis = DEM_analysis_stage.DEMAnalysisStage

class ThermalDEMAnalysis(BaseAnalysis):

    def __init__(self, model, parameters):
        pass