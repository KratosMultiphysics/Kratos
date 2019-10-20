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
from KratosMultiphysics.ShapeOptimizationApplication.analyzer_base import AnalyzerBaseClass
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

    dependency_graph, exist_dependencies = _CreateDependencyGraphRecursively(objectives)
    # dependency_graph = [ (response_id_1, [], weight_1),
    #                      (response_id_2, [], weight_2),
    #                      (response_id_3, [ (response_id_3a, [], weight_3a),
    #                                        (response_id_3b, [], weight_3b),
    #                                        (response_id_3c, [...], weight_3c) ], weight_3),
    #                      ... ]

    if exist_dependencies:
        return AnalyzerWithDependencies(internal_analyzer, model_part_controller, external_analyzer, dependency_graph)
    else:
        return Analyzer(internal_analyzer, model_part_controller, external_analyzer)

# ------------------------------------------------------------------------------
def _IdentifyInternalResponsesRecursively(responses):
    internal_responses = []

    for itr in range(responses.size()):
        response_i = responses[itr]

        if response_i.Has("is_combined"):
            if response_i["is_combined"].GetBool():
                internal_responses += _IdentifyInternalResponsesRecursively(response_i["combined_responses"])

        if response_i["analyzer"].GetString() == "kratos":
            identifier = response_i["identifier"].GetString()
            kratos_settings = response_i["response_settings"]
            internal_responses.append([identifier, kratos_settings])

    return internal_responses

# ------------------------------------------------------------------------------
def _CreateDependencyGraphRecursively(responses):
    dependency_graph = []
    exist_dependencies = False

    for itr in range(responses.size()):
        response_i = responses[itr]
        identifier = response_i["identifier"].GetString()
        weight = response_i["weight"].GetDouble()

        if response_i.Has("is_combined"):
            if response_i["is_combined"].GetBool():
                exist_dependencies = True
                sub_dependency_graph, _ = _CreateDependencyGraphRecursively(response_i["combined_responses"])
                dependency_graph.append((identifier, sub_dependency_graph, weight))
            else:
                dependency_graph.append((identifier, [], weight))

    return dependency_graph, exist_dependencies

# ==============================================================================
class Analyzer(AnalyzerBaseClass):
    # --------------------------------------------------------------------------
    def __init__(self, internal_analyzer, model_part_controller, external_analyzer):
        self.model_part_controller = model_part_controller
        self.internal_analyzer = internal_analyzer
        self.external_analyzer = external_analyzer

        if isinstance(internal_analyzer, EmptyAnalyzer) and  isinstance(external_analyzer, EmptyAnalyzer):
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
class AnalyzerWithDependencies(Analyzer):
    # --------------------------------------------------------------------------
    def __init__(self, internal_analyzer, model_part_controller, external_analyzer, dependency_graph):
        super().__init__(internal_analyzer, model_part_controller, external_analyzer)

        self.dependency_graph = dependency_graph
        self.response_combination_filename = "response_combination.csv"
        self.gradient_max_norms = {}

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop(self):
        super().InitializeBeforeOptimizationLoop()

        self.__InitializeOutputOfResponses()

    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, unique_iterator, communicator):
        self.__RequestResponsesAccordingDependencies(communicator)

        super().AnalyzeDesignAndReportToCommunicator(current_design, unique_iterator, communicator)

        self.__CombineResponsesAccordingDependencies(communicator)
        self.__ComputeGradientNorms(communicator)
        self.__WriteResultsOfCombinedResponses(unique_iterator,communicator)

    # --------------------------------------------------------------------------
    def __InitializeOutputOfResponses(self):
        kratos_utilities.DeleteFileIfExisting(self.response_combination_filename)

        with open(self.response_combination_filename, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            response_ids = self.__GetIdentifiersRecursively(self.dependency_graph)

            writer.writerow(["---------------------------------"])
            for itr, identifier in enumerate(response_ids):
                writer.writerow(["f"+str(itr)+": "+identifier])
            writer.writerow(["---------------------------------"])

            row = []
            row.append("{:>4s}".format("itr"))
            for itr, identifier in enumerate(response_ids):
                row.append("{:>13s}".format("f"+str(itr)+"_value"))
            for itr, identifier in enumerate(response_ids):
                row.append("{:>13s}".format("||df"+str(itr)+"dx_st||"))
            writer.writerow(row)

    # --------------------------------------------------------------------------
    def __GetIdentifiersRecursively(self, dependencies):
        response_ids = []

        for response_id, dependencies, _ in dependencies:
            response_ids += [response_id]
            if len(dependencies) > 0:
                sub_identifiers = self.__GetIdentifiersRecursively(dependencies)
                response_ids += sub_identifiers

        return response_ids

    # --------------------------------------------------------------------------
    def __RequestResponsesAccordingDependencies(self, communicator):
        for response_id, dependencies, _ in self.dependency_graph:
            sub_response_ids = self.__GetIdentifiersRecursively(dependencies)

            if communicator.isRequestingValueOf(response_id):
                for sub_response_id in sub_response_ids:
                    communicator.requestValueOf(sub_response_id)

            if communicator.isRequestingGradientOf(response_id):
                for sub_response_id in sub_response_ids:
                    communicator.requestGradientOf(sub_response_id)

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

            for vector in gradient.values():
                vector[0] *= weight
                vector[1] *= weight
                vector[2] *= weight

            if combined_gradient is None:
                combined_gradient = gradient
            else:
                # Perform nodal sum
                for key_a, vector_a in combined_gradient.items():
                    vector_b = gradient[key_a]
                    vector_a[0] += vector_b[0]
                    vector_a[1] += vector_b[1]
                    vector_a[2] += vector_b[2]

        return combined_gradient

    # --------------------------------------------------------------------------
    def __ComputeGradientNorms(self, communicator):
        response_ids = self.__GetIdentifiersRecursively(self.dependency_graph)

        for response_id in response_ids:
            gradient = communicator.getStandardizedGradient(response_id)

            nodal_norms = [ entry[0]**2 + entry[1]**2 + entry[2]**2 for entry in gradient.values() ]
            max_norm = math.sqrt(max(nodal_norms))
            self.gradient_max_norms[response_id] = max_norm

    # --------------------------------------------------------------------------
    def __WriteResultsOfCombinedResponses(self, iteration, communicator):
        with open(self.response_combination_filename, 'a') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')

            identifiers = self.__GetIdentifiersRecursively(self.dependency_graph)
            values = []
            for response_id in identifiers:
                values.append(communicator.getValue(response_id))

            row = []
            row.append("{:>4d}".format(iteration))
            for identifier, value in zip(identifiers, values):
                row.append(" {:> .5E}".format(value))
            for identifier, value in zip(identifiers, values):
                row.append(" {:> .5E}".format(self.gradient_max_norms[identifier]))

            writer.writerow(row)

# ==============================================================================