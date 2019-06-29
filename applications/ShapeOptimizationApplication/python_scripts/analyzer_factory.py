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
    # internal_responses = [ (response_id_1, kratos_settings_1),
    #                        (response_id_2, kratos_settings_2),
    #                        ... ]

    if len(internal_responses) > 0:
        internal_analyzer = KratosInternalAnalyzer(internal_responses, model_part_controller)
    else:
        internal_analyzer = EmptyAnalyzer()

    dependency_graph = _CreateDependencyGraphRecursively(objectives)
    # dependency_graph = [ (response_id_1, [], weight_1),
    #                      (response_id_2, [], weight_2),
    #                      (response_id_3, [ (response_id_3a, [], weight_3a),
    #                                        (response_id_3b, [], weight_3b),
    #                                        (response_id_3c, [...], weight_3c) ], weight_3),
    #                      ... ]

    return Analyzer(internal_analyzer, model_part_controller, external_analyzer, dependency_graph)

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

# ------------------------------------------------------------------------------
def _CreateDependencyGraphRecursively(responses):
    dependency_graph = []

    for itr in range(responses.size()):
        response_i = responses[itr]
        identifier = response_i["identifier"].GetString()
        weight = response_i["weight"].GetDouble()

        if response_i.Has("is_combined"):
            if response_i["is_combined"].GetBool():
                sub_dependency_graph = _CreateDependencyGraphRecursively(response_i["combined_responses"])
                dependency_graph.append((identifier, sub_dependency_graph, weight))
            else:
                dependency_graph.append((identifier, [], weight))

    return dependency_graph

                combinations[identifier] = combi_list

                _IdentifyAndAddResponseCombinationsRecursively(combined_responses, combinations)

# ==============================================================================
class Analyzer:
    # --------------------------------------------------------------------------
    def __init__(self, internal_analyzer, model_part_controller, external_analyzer, dependency_graph = []):
        self.model_part_controller = model_part_controller
        self.internal_analyzer = internal_analyzer
        self.external_analyzer = external_analyzer
        self.dependency_graph = dependency_graph

        if internal_analyzer.IsEmpty() and external_analyzer.IsEmpty():
            raise RuntimeError("Neither an internal nor an external analyzer is defined!")

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop(self):
        self.internal_analyzer.InitializeBeforeOptimizationLoop()
        self.external_analyzer.InitializeBeforeOptimizationLoop()

    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, unique_iterator, communicator):

        for combination_id, functions in self.response_combinations.items():
            if communicator.isRequestingValueOf(combination_id):
                self.__RequestResponsesValuesRecursively(functions, communicator)

            if communicator.isRequestingGradientOf(combination_id):
                self.__RequestResponsesGradientsRecursively(functions, communicator)

        # Analyze design w.r.t to all given response functions
        self.internal_analyzer.AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)
        self.external_analyzer.AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)

        # Combine response functions if necessary
        for combination_id, functions in self.response_combinations.items():
            if communicator.isRequestingValueOf(combination_id):
                combined_value = self.__ComputeCombinationValueRecursively(functions, communicator)
                communicator.reportValue(combination_id, combined_value)

            # if communicator.isRequestingGradientOf(combination_id):
            #     combined_gradient = None
            #     for (function_id, weight) in functions:
            #         if combined_gradient is None:
            #             combined_gradient = communicator.getStandardizedGradient(function_id)
            #         else:
            #             combined_to_add = communicator.getStandardizedGradient(function_id)
            #             combined_gradient.update(...)

            #     communicator.reportGradient(function_id, combined_gradient)

        self.__ResetPossibleShapeModificationsFromAnalysis()

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop(self):
        self.internal_analyzer.FinalizeAfterOptimizationLoop()
        self.external_analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __RequestResponsesValuesRecursively(self, functions, communicator):
        for function_id, _ in functions:
            if function_id in self.response_combinations:
                communicator.requestValueOf(function_id)
                self.__RequestResponsesValuesRecursively(self.response_combinations[function_id], communicator)
            else:
                communicator.requestValueOf(function_id)

    # --------------------------------------------------------------------------
    def __RequestResponsesGradientsRecursively(self, functions, communicator):
        for function_id, _ in functions:
            if function_id in self.response_combinations:
                communicator.requestGradientOf(function_id)
                self.__RequestResponsesGradientsRecursively(self.response_combinations[function_id], communicator)
            else:
                communicator.requestGradientOf(function_id)

    # # --------------------------------------------------------------------------
    # def __RequestResponsesFromCombinationsRecursively(self, communicator):
    #         if combination_value_required:
    #             for function_id, _ in functions:
    #                 if function_id in self.response_combinations:
    #                     communicator.requestValueOf(function_id)
    #                     self.__RequestResponsesFromCombinationsRecursively(communicator)
    #                 else:
    #                     communicator.requestValueOf(function_id)

    #         if combination_gradient_required:
    #             for function_id, _ in functions:
    #                 if function_id in self.response_combinations:
    #                     communicator.requestGradientOf(function_id)
    #                     self.__RequestResponsesFromCombinationsRecursively(communicator)
    #                 else:
    #                     communicator.requestGradientOf(function_id)

    # --------------------------------------------------------------------------
    def __ComputeCombinationValueRecursively(self, functions, communicator):
        combined_value = 0.0

        for (function_id, weight) in functions:
            if function_id in self.response_combinations:
                value = self.__ComputeCombinationValueRecursively(self.response_combinations[function_id], communicator)
                communicator.reportValue(function_id, value)
            else:
                value = communicator.getStandardizedValue(function_id)
            combined_value += weight*value

        return combined_value

    # --------------------------------------------------------------------------
    def __ResetPossibleShapeModificationsFromAnalysis( self ):
        self.model_part_controller.SetMeshToReferenceMesh()
        self.model_part_controller.SetDeformationVariablesToZero()

# ==============================================================================
