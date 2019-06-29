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

        if len(self.dependency_graph) != 0:
            self.__RequestResponsesFromDependencies(communicator)

        self.internal_analyzer.AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)
        self.external_analyzer.AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)

        if len(self.dependency_graph) != 0:
            self.__CombineResponsesAccordingDependencies(communicator)

        self.__ResetPossibleShapeModificationsFromAnalysis()

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop(self):
        self.internal_analyzer.FinalizeAfterOptimizationLoop()
        self.external_analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __RequestResponsesFromDependencies(self, communicator):
        for response_id, dependencies, _ in self.dependency_graph:
            if len(dependencies) > 0:
                if communicator.isRequestingValueOf(response_id):
                    self.__RequestResponsesValuesRecursively(dependencies, communicator)

                if communicator.isRequestingGradientOf(response_id):
                    self.__RequestResponsesGradientsRecursively(dependencies, communicator)

    # --------------------------------------------------------------------------
    def __RequestResponsesValuesRecursively(self, dependencies, communicator):
        for response_id, dependencies, _ in dependencies:
            if len(dependencies) > 0:
                self.__RequestResponsesValuesRecursively(dependencies, communicator)
            else:
                communicator.requestValueOf(response_id)

    # --------------------------------------------------------------------------
    def __RequestResponsesGradientsRecursively(self, dependencies, communicator):
        for response_id, dependencies, _ in dependencies:
            if len(dependencies) > 0:
                self.__RequestResponsesGradientsRecursively(dependencies, communicator)
            else:
                communicator.requestGradientOf(response_id)

    # --------------------------------------------------------------------------
    def __CombineResponsesAccordingDependencies(self, communicator):
        for response_id, dependencies, _ in self.dependency_graph:
            if len(dependencies) > 0:
                if communicator.isRequestingValueOf(response_id):
                    combined_value = self.__ComputeCombinationValueRecursively(dependencies, communicator)
                    communicator.reportValue(response_id, combined_value)

        for response_id, dependencies, _ in self.dependency_graph:
            if len(dependencies) > 0:
                if communicator.isRequestingGradientOf(response_id):
                    combined_gradient = self.__ComputeCombinationGradientRecursively(dependencies, communicator)
                    communicator.reportGradient(response_id, combined_gradient)

    # --------------------------------------------------------------------------
    def __ComputeCombinationValueRecursively(self, dependencies, communicator):
        combined_value = 0.0

        for response_id, dependencies, weight in dependencies:
            if len(dependencies) > 0:
                value = self.__ComputeCombinationValueRecursively(dependencies, communicator)
                communicator.reportValue(response_id, value)
            else:
                value = communicator.getStandardizedValue(response_id)
            combined_value += weight*value

        return combined_value

    # --------------------------------------------------------------------------
    def __ComputeCombinationGradientRecursively(self, dependencies, communicator):
        combined_gradient = None

        for response_id, dependencies, weight in dependencies:
            if len(dependencies) > 0:
                gradient = self.__ComputeCombinationGradientRecursively(dependencies, communicator)
                communicator.reportGradient(response_id, gradient)
            else:
                gradient = communicator.getStandardizedGradient(response_id)

            if combined_gradient is None:
                combined_gradient = gradient
                combined_gradient.update({key: [weight*value[0],weight*value[1],weight*value[2]] for key, value in gradient.items()})
            else:
                # Perform nodal sum
                update = {key_a: [a+b for a, b in zip(list_a, list_b)] for ((key_a, list_a),(key_b, list_b)) in zip(combined_gradient.items(), gradient.items())}
                combined_gradient.update( update )

        return combined_gradient

    # --------------------------------------------------------------------------
    def __ResetPossibleShapeModificationsFromAnalysis( self ):
        self.model_part_controller.SetMeshToReferenceMesh()
        self.model_part_controller.SetDeformationVariablesToZero()

# ==============================================================================
