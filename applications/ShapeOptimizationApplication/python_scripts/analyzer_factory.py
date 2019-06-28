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

# additional imports
from KratosMultiphysics.ShapeOptimizationApplication.analyzer_internal import KratosInternalAnalyzer
from KratosMultiphysics.ShapeOptimizationApplication.analyzer_empty import EmptyAnalyzer

# ==============================================================================
def CreateAnalyzer(optimization_settings, model_part_controller, external_analyzer):
    if _IsInternalAnalyzerRequired(optimization_settings):
        internal_analyzer = _CreateInternalAnalyzer(optimization_settings, model_part_controller)
    else:
        internal_analyzer = EmptyAnalyzer()

    return Analyzer(internal_analyzer, model_part_controller, external_analyzer)

# --------------------------------------------------------------------------
def _IsInternalAnalyzerRequired(optimization_settings):
    for itr in range(optimization_settings["objectives"].size()):
        if optimization_settings["objectives"][itr]["use_kratos"].GetBool():
            return True

    for itr in range(optimization_settings["constraints"].size()):
        if optimization_settings["constraints"][itr]["use_kratos"].GetBool():
            return True
    return False

# ------------------------------------------------------------------------------
def _CreateInternalAnalyzer(optimization_settings, model_part_controller):
    specified_response_functions = {}

    for itr in range(optimization_settings["objectives"].size()):
        objective_settings = optimization_settings["objectives"][itr]

        if objective_settings["use_kratos"].GetBool():
            identifier = objective_settings["identifier"].GetString()
            kratos_settings = objective_settings["kratos_response_settings"]
            specified_response_functions[identifier] = kratos_settings

    for itr in range(optimization_settings["constraints"].size()):
        constraint_settings = optimization_settings["constraints"][itr]

        if constraint_settings["use_kratos"].GetBool():
            identifier = constraint_settings["identifier"].GetString()
            kratos_settings = constraint_settings["kratos_response_settings"]
            specified_response_functions[identifier] = kratos_settings

    return KratosInternalAnalyzer(specified_response_functions, model_part_controller)

# ==============================================================================
class Analyzer:
    # --------------------------------------------------------------------------
    def __init__(self, internal_analyzer, model_part_controller, external_analyzer):
        self.model_part_controller = model_part_controller
        self.internal_analyzer = internal_analyzer
        self.external_analyzer = external_analyzer

        if internal_analyzer.IsEmpty() and external_analyzer.IsEmpty():
            raise RuntimeError("Neither an internal nor an external analyzer is defined!")

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop(self):
        self.internal_analyzer.InitializeBeforeOptimizationLoop()
        self.external_analyzer.InitializeBeforeOptimizationLoop()

    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, unique_iterator, communicator):
        self.internal_analyzer.AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)
        self.external_analyzer.AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)

        self.__ResetPossibleShapeModificationsFromAnalysis()

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop(self):
        self.internal_analyzer.FinalizeAfterOptimizationLoop()
        self.external_analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __ResetPossibleShapeModificationsFromAnalysis( self ):
        self.model_part_controller.SetMeshToReferenceMesh()
        self.model_part_controller.SetDeformationVariablesToZero()

# ==============================================================================
