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

# ==============================================================================
def CreateAnalyzer(optimization_settings, model_part_controller, external_analyzer):
    return Analyzer(optimization_settings, model_part_controller, external_analyzer)

# ==============================================================================
class Analyzer:
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, model_part_controller, external_analyzer):
        self.external_analyzer = external_analyzer

        if self.__IsInternalAnalyzerRequired(optimization_settings):
            from analyzer_internal import KratosInternalAnalyzer
            self.internal_analyzer = KratosInternalAnalyzer(optimization_settings, model_part_controller.GetOptimizationModelPart())
        else:
            from analyzer_empty import EmptyAnalyzer
            self.internal_analyzer = EmptyAnalyzer()
            if isinstance(external_analyzer, EmptyAnalyzer):
                raise RuntimeError("Neither an internal nor an external analyzer is defined!")

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop(self):
        self.external_analyzer.InitializeBeforeOptimizationLoop()
        self.internal_analyzer.InitializeBeforeOptimizationLoop()

    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, unique_iterator, communicator):
        self.external_analyzer.AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)
        self.internal_analyzer.AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop(self):
        self.external_analyzer.FinalizeAfterOptimizationLoop()
        self.internal_analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __IsInternalAnalyzerRequired(self, optimization_settings):
        for objective_number in range(optimization_settings["objectives"].size()):
            if optimization_settings["objectives"][objective_number]["use_kratos"].GetBool():
                return True

        for constraint_number in range(optimization_settings["constraints"].size()):
            if optimization_settings["constraints"][constraint_number]["use_kratos"].GetBool():
                return True
        return False

# ==============================================================================
