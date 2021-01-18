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
def CreateCommunicator(optimization_settings):
    return Communicator(optimization_settings)

# ==============================================================================
class Communicator:

    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings):
        self.optimization_settings = optimization_settings
        self.supported_objective_types = ["minimization", "maximization"]
        self.supported_constraint_types = ["=", "<", ">", "<=", ">="]
        self.supported_constraint_references = ["initial_value", "specified_value"]

        objective_settings = self.__ExtractResponseSettingsRecursively(optimization_settings["objectives"])
        constraint_settings = self.__ExtractResponseSettingsRecursively(optimization_settings["constraints"])
        all_response_settings = objective_settings + constraint_settings

        self.list_of_requests = self.__initializeListOfRequests(all_response_settings)
        self.list_of_responses = self.__initializeListOfObjectives(objective_settings)
        self.list_of_responses.update(self.__initializeListOfConstraints(constraint_settings))

    # --------------------------------------------------------------------------
    def initializeCommunication(self):
        self.__deleteAllRequests()
        self.__deleteAllReportedValues()

    # --------------------------------------------------------------------------
    def updateRequest(self, response_id, **kwargs):
        if "request_value" in kwargs.keys():
            self.list_of_requests[response_id]["calculateValue"] = kwargs["request_value"]
        if "request_gradient" in kwargs.keys():
            self.list_of_requests[response_id]["calculateGradient"] = kwargs["request_gradient"]

    # --------------------------------------------------------------------------
    def requestValueOf(self, response_id):
        self.list_of_requests[response_id]["calculateValue"] = True

    # --------------------------------------------------------------------------
    def requestGradientOf(self, response_id):
        self.list_of_requests[response_id]["calculateGradient"] = True

    # --------------------------------------------------------------------------
    def isRequestingValueOf(self, response_id):
        if response_id not in self.list_of_requests.keys():
            return False
        else:
            return self.list_of_requests[response_id]["calculateValue"]

    # --------------------------------------------------------------------------
    def isRequestingGradientOf(self, response_id):
        if response_id not in self.list_of_requests.keys():
            return False
        else:
            return self.list_of_requests[response_id]["calculateGradient"]

    # --------------------------------------------------------------------------
    def reportValue(self, response_id, value):
        self.__storeValue(response_id, value)
        if self.__isResponseWaitingForInitialValueAsReference(response_id):
            self.__setValueAsReference(response_id, value)

        standardized_value = self.__translateValueToStandardForm(response_id, value)
        self.__storeStandardizedValue(response_id, standardized_value)

    # --------------------------------------------------------------------------
    def reportGradient(self, response_id, gradient):
        dimension = len(next(iter(gradient.values())))
        if dimension == 1:
            gradient = {key: [value[0], 0.0, 0.0] for key, value in gradient.items()}
        elif dimension == 2:
            gradient = {key: [value[0], value[1], 0.0] for key, value in gradient.items()}

        standardized_gradient = self.__translateGradientToStandardForm(response_id, gradient)
        self.__storeStandardizedGradient(response_id, standardized_gradient)

    # --------------------------------------------------------------------------
    def getValue(self, response_id):
        if self.list_of_responses[response_id]["value"] is None:
            raise RuntimeError("The requested value for ", response_id, " is None!")
        return self.list_of_responses[response_id]["value"]

    # --------------------------------------------------------------------------
    def getReferenceValue(self, response_id):
        if self.list_of_responses[response_id]["reference_value"] is None:
            raise RuntimeError("The requested reference_value for ", response_id, " is None!")
        return self.list_of_responses[response_id]["reference_value"]

    # --------------------------------------------------------------------------
    def getStandardizedValue(self, response_id):
        if self.list_of_responses[response_id]["standardized_value"] is None:
            raise RuntimeError("The requested standardized_value for ", response_id, " is None!")
        return self.list_of_responses[response_id]["standardized_value"]

    # --------------------------------------------------------------------------
    def getStandardizedGradient(self, response_id):
        if self.list_of_responses[response_id]["standardized_gradient"] is None:
            raise RuntimeError("The requested standardized_gradient for ", response_id, " is None!")
        return self.list_of_responses[response_id]["standardized_gradient"]

    # --------------------------------------------------------------------------
    @classmethod
    def __ExtractResponseSettingsRecursively(cls, response_settings):
        list_of_settings = []

        for itr in range(response_settings.size()):
            response_i = response_settings[itr]
            if response_i.Has("is_combined"):
                if response_i["is_combined"].GetBool():
                    list_of_settings += cls.__ExtractResponseSettingsRecursively(response_i["combined_responses"])

            list_of_settings += [response_i]

        return list_of_settings

    # --------------------------------------------------------------------------
    @classmethod
    def __initializeListOfRequests(cls, response_settings):
        list_of_requests = {}

        for response in response_settings:
            response_id = response["identifier"].GetString()
            list_of_requests[response_id] = {"calculateValue": False, "calculateGradient": False}

        return list_of_requests

    # --------------------------------------------------------------------------
    def __initializeListOfObjectives(self, objective_settings):
        list_of_objectives = {}

        for objective in objective_settings:
            objective_id =  objective["identifier"].GetString()

            if objective["type"].GetString() not in self.supported_objective_types:
                raise RuntimeError("Unsupported type defined for the following objective: " + objective_id)

            list_of_objectives[objective_id] = { "type"                 : objective["type"].GetString(),
                                                 "value"                : None,
                                                 "scaling_factor"       : objective["scaling_factor"].GetDouble(),
                                                 "standardized_value"   : None,
                                                 "standardized_gradient": None }

        return list_of_objectives

    # --------------------------------------------------------------------------
    def __initializeListOfConstraints(self, constraint_settings):
        list_of_constraints = {}

        for constraint in constraint_settings:
            constraint_id = constraint["identifier"].GetString()

            if constraint["type"].GetString() not in self.supported_constraint_types:
                raise RuntimeError("Unsupported type defined for the following constraint: " + constraint_id)

            if constraint["reference"].GetString() not in self.supported_constraint_references:
                raise RuntimeError("Unsupported reference defined for the following constraint: " + constraint_id)

            if  constraint["reference"].GetString() == "specified_value":
                list_of_constraints[constraint_id] = { "type"                 : constraint["type"].GetString(),
                                                        "value"                : None,
                                                        "scaling_factor"       : constraint["scaling_factor"].GetDouble(),
                                                        "standardized_value"   : None,
                                                        "standardized_gradient": None,
                                                        "reference_value"      : constraint["reference_value"].GetDouble() }
            elif constraint["reference"].GetString() == "initial_value":
                list_of_constraints[constraint_id] = { "type"                 : constraint["type"].GetString(),
                                                        "value"                : None,
                                                        "scaling_factor"       : constraint["scaling_factor"].GetDouble(),
                                                        "standardized_value"   : None,
                                                        "standardized_gradient": None,
                                                        "reference_value"      : None }
            else:
                raise RuntimeError("Unsupported reference defined for the following constraint: " + constraint_id)

        return list_of_constraints

    # --------------------------------------------------------------------------
    def __deleteAllRequests(self):
        for response_id in self.list_of_requests:
            self.list_of_requests[response_id]["calculateValue"] = False
            self.list_of_requests[response_id]["calculateGradient"] = False

    # --------------------------------------------------------------------------
    def __deleteAllReportedValues(self):
        for response_id in self.list_of_responses:
            self.list_of_responses[response_id]["value"] = None
            self.list_of_responses[response_id]["standardized_value"] = None
            self.list_of_responses[response_id]["standardized_gradient"] = None

    # --------------------------------------------------------------------------
    def __isResponseWaitingForInitialValueAsReference(self, response_id):
        response = self.list_of_responses[response_id]
        if "reference_value" in response:
            is_reference_value_missing = (response["reference_value"] is None)
            if is_reference_value_missing:
                return True
        return False

    # --------------------------------------------------------------------------
    def __setValueAsReference(self, response_id, value):
        self.list_of_responses[response_id]["reference_value"] = value

    # --------------------------------------------------------------------------
    def __translateValueToStandardForm(self, response_id, value):
        response_type = self.list_of_responses[response_id]["type"]
        scaling_factor = self.list_of_responses[response_id]["scaling_factor"]

        if response_type in self.supported_objective_types:
            if response_type == "maximization":
                return -scaling_factor*value
            else:
                return scaling_factor*value
        else:
            reference_value = self.list_of_responses[response_id]["reference_value"]

            if response_type == ">" or response_type == ">=":
                return scaling_factor*(reference_value-value)
            else:
                return scaling_factor*(value-reference_value)

    # --------------------------------------------------------------------------
    def __translateGradientToStandardForm(self, response_id, gradient):
        response_type = self.list_of_responses[response_id]["type"]
        scaling_factor = self.list_of_responses[response_id]["scaling_factor"]

        if response_type == "maximization" or response_type == ">" or response_type == ">=":
            scaling_factor *= -1

        for vector in gradient.values():
            vector[0] *= scaling_factor
            vector[1] *= scaling_factor
            vector[2] *= scaling_factor

        return gradient

    # --------------------------------------------------------------------------
    def __storeValue(self, response_id, value):
        self.list_of_responses[response_id]["value"] = value

    # --------------------------------------------------------------------------
    def __storeStandardizedValue(self, response_id, value):
        self.list_of_responses[response_id]["standardized_value"] = value

    # --------------------------------------------------------------------------
    def __storeStandardizedGradient(self, response_id, gradient):
        self.list_of_responses[response_id]["standardized_gradient"] = gradient

# ==============================================================================
