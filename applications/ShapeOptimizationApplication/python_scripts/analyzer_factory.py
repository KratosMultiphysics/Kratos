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
    objectives = optimization_settings["objectives"]
    constraints = optimization_settings["constraints"]

    internal_objectives = _IdentifyInternalResponsesRecursively(objectives)
    internal_constraints = _IdentifyInternalResponsesRecursively(constraints)
    internal_responses = internal_objectives + internal_constraints

    if len(internal_responses) > 0:
        internal_analyzer = KratosInternalAnalyzer(internal_responses, model_part_controller)
    else:
        internal_analyzer = EmptyAnalyzer()

    return Analyzer(internal_analyzer, model_part_controller, external_analyzer)

# ------------------------------------------------------------------------------
def _IdentifyInternalResponsesRecursively(responses):
    internal_responses = []

    for itr in range(responses.size()):
        response_i = responses[itr]

        if response_i.Has("is_combined"):
            if response_i["is_combined"].GetBool():
                internal_responses += _IdentifyInternalResponsesRecursively(response_i["combined_responses"])

        if response_i["use_kratos"].GetBool():
            identifier = response_i["identifier"].GetString()
            kratos_settings = response_i["kratos_response_settings"]
            internal_responses.append([identifier, kratos_settings])

    return internal_responses

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
