# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Additional imports
from analyzer_internal import KratosInternalAnalyzer

# ==============================================================================
def CreateAnalyzer( project_parameters, model_part_controller, external_anlyzer ):
    return CombinedAnalyzer(project_parameters, model_part_controller, external_anlyzer)

# ==============================================================================
class CombinedAnalyzer:

    # --------------------------------------------------------------------------
    def __init__( self, project_parameters, model_part_controller, external_anlyzer ):
        self.external_anlyzer = external_anlyzer
        self.internal_analyzer = KratosInternalAnalyzer(project_parameters, model_part_controller.GetOptimizationModelPart())

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop( self ):
        self.internal_analyzer.InitializeBeforeOptimizationLoop()
        if self.external_anlyzer is not None:
            self.external_anlyzer.InitializeBeforeOptimizationLoop()

    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator( self, current_design, unique_iterator, communicator ):
        self.internal_analyzer.AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)
        if self.external_anlyzer is not None:
            self.external_anlyzer.AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop( self ):
        self.internal_analyzer.FinalizeAfterOptimizationLoop()
        if self.external_anlyzer is not None:
            self.external_anlyzer.FinalizeAfterOptimizationLoop()

# ==============================================================================
