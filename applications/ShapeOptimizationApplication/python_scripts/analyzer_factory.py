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
import KratosMultiphysics.kratos_utilities as kratos_utilities
import csv, math

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

        self.exist_dependencies = False
        for _, dependencies, _ in dependency_graph:
            if len(dependencies) > 0:
                self.exist_dependencies = True

        self.response_values_filename = "response_values_from_analyzer.csv"
        self.gradients_norms = {}

        if internal_analyzer.IsEmpty() and external_analyzer.IsEmpty():
            raise RuntimeError("Neither an internal nor an external analyzer is defined!")

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop(self):
        if len(self.dependency_graph) != 0:
            self.__InitializeOutputOfCombinedResponses()

        self.internal_analyzer.InitializeBeforeOptimizationLoop()
        self.external_analyzer.InitializeBeforeOptimizationLoop()

    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, unique_iterator, communicator):
        if self.exist_dependencies:
            self.__RequestResponsesAccordingDependencies(communicator)

        self.internal_analyzer.AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)
        self.external_analyzer.AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)

        if self.exist_dependencies:
            self.__CombineResponsesAccordingDependencies(communicator)
            self.__ComputeGradientNormsRecursively(self.dependency_graph, communicator)
            self.__WriteResultsOfCombinedResponses(unique_iterator,communicator)

        self.__ResetPossibleShapeModificationsFromAnalysis()

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop(self):
        self.internal_analyzer.FinalizeAfterOptimizationLoop()
        self.external_analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __InitializeOutputOfCombinedResponses(self):
        kratos_utilities.DeleteFileIfExisting(self.response_values_filename)

        with open(self.response_values_filename, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=',',quotechar='|',escapechar=' ',quoting=csv.QUOTE_MINIMAL)
            identifiers = self.__ReadIdentifiersFromDependencyGraphRecursively(self.dependency_graph)

            row = ["iteration"] + identifiers
            writer.writerow(row)

    # --------------------------------------------------------------------------
    def __RequestResponsesAccordingDependencies(self, communicator):
        for response_id, dependencies, _ in self.dependency_graph:
            if len(dependencies) > 0:
                if communicator.isRequestingValueOf(response_id):
                    self.__RequestValuesRecursively(dependencies, communicator)

                if communicator.isRequestingGradientOf(response_id):
                    self.__RequestGradientsRecursively(dependencies, communicator)

    # --------------------------------------------------------------------------
    def __RequestValuesRecursively(self, dependencies, communicator):
        for response_id, dependencies, _ in dependencies:
            if len(dependencies) > 0:
                self.__RequestValuesRecursively(dependencies, communicator)
            else:
                communicator.requestValueOf(response_id)

    # --------------------------------------------------------------------------
    def __RequestGradientsRecursively(self, dependencies, communicator):
        for response_id, dependencies, _ in dependencies:
            if len(dependencies) > 0:
                self.__RequestGradientsRecursively(dependencies, communicator)
            else:
                communicator.requestGradientOf(response_id)

    # --------------------------------------------------------------------------
    def __CombineResponsesAccordingDependencies(self, communicator):
        for response_id, dependencies, _ in self.dependency_graph:
            if len(dependencies) > 0:
                if communicator.isRequestingValueOf(response_id):
                    combined_value = self.__ComputeCombinedValuesRecursively(dependencies, communicator)
                    communicator.reportValue(response_id, combined_value)

        for response_id, dependencies, _ in self.dependency_graph:
            if len(dependencies) > 0:
                if communicator.isRequestingGradientOf(response_id):
                    combined_gradient = self.__ComputeCombinedGradientsRecursively(dependencies, communicator)
                    communicator.reportGradient(response_id, combined_gradient)

    # --------------------------------------------------------------------------
    def __ComputeCombinedValuesRecursively(self, dependencies, communicator):
        combined_value = 0.0

        for response_id, dependencies, weight in dependencies:
            if len(dependencies) > 0:
                value = self.__ComputeCombinedValuesRecursively(dependencies, communicator)
                communicator.reportValue(response_id, value)
            else:
                value = communicator.getStandardizedValue(response_id)
            combined_value += weight*value

        return combined_value

    # --------------------------------------------------------------------------
    def __ComputeCombinedGradientsRecursively(self, dependencies, communicator):
        combined_gradient = None

        for response_id, dependencies, weight in dependencies:
            if len(dependencies) > 0:
                gradient = self.__ComputeCombinedGradientsRecursively(dependencies, communicator)
                communicator.reportGradient(response_id, gradient)
            else:
                gradient = communicator.getStandardizedGradient(response_id)

            gradient = {key: [weight*value[0],weight*value[1],weight*value[2]] for key, value in gradient.items()}
            if combined_gradient is None:
                combined_gradient = gradient
            else:
                # Perform nodal sum
                combined_gradient = {key_a: [a+b for a, b in zip(list_a, list_b)] for ((key_a, list_a),(key_b, list_b)) in zip(combined_gradient.items(), gradient.items())}

        return combined_gradient

    # --------------------------------------------------------------------------
    def __ComputeGradientNormsRecursively(self, dependencies, communicator):
        for response_id, dependencies, _ in dependencies:
            gradient = communicator.getStandardizedGradient(response_id)

            nodal_norms = [ entry[0]**2 + entry[1]**2 + entry[2]**2 for entry in gradient.values() ]
            max_norm = math.sqrt(max(nodal_norms))
            self.gradients_norms[response_id] = max_norm

            if len(dependencies) > 0:
                self.__ComputeGradientNormsRecursively(dependencies, communicator)

    # --------------------------------------------------------------------------
    def __WriteResultsOfCombinedResponses(self, iteration, communicator):
        with open(self.response_values_filename, 'a') as csvfile:
            writer = csv.writer(csvfile, delimiter=',',quotechar='|',escapechar=' ',quoting=csv.QUOTE_MINIMAL)

            values = self.__ReadValuesFromDependencyGraphRecursively(self.dependency_graph, communicator)
            values = ["{:>.5E}".format(entry) for entry in values]

            row = [iteration] + values
            writer.writerow(row)

    # --------------------------------------------------------------------------
    def __ReadIdentifiersFromDependencyGraphRecursively(self, dependencies):
        identifiers = []

        for response_id, dependencies, weight in dependencies:
            identifiers += [response_id]
            if len(dependencies) > 0:
                sub_identifiers = self.__ReadIdentifiersFromDependencyGraphRecursively(dependencies)
                identifiers += sub_identifiers

        return identifiers

    # --------------------------------------------------------------------------
    def __ReadValuesFromDependencyGraphRecursively(self, dependencies, communicator):
        values = []

        for response_id, dependencies, weight in dependencies:
            values += [communicator.getValue(response_id)]
            if len(dependencies) > 0:
                sub_values = self.__ReadValuesFromDependencyGraphRecursively(dependencies, communicator)
                values += sub_values

        return values

    # --------------------------------------------------------------------------
    def __ResetPossibleShapeModificationsFromAnalysis( self ):
        self.model_part_controller.SetMeshToReferenceMesh()
        self.model_part_controller.SetDeformationVariablesToZero()

# ==============================================================================
