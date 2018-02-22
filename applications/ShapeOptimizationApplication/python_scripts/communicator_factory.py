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

# ==============================================================================
def CreateCommunicator( optimization_settings ):
    return Communicator( optimization_settings )

# ==============================================================================
class Communicator:

    # --------------------------------------------------------------------------
    def __init__( self, optimization_settings ):
        self.optimization_settings = optimization_settings
        self.__initializeListOfRequests()
        self.__initializeListOfResponses()

    # --------------------------------------------------------------------------
    def initializeCommunication( self ):
        self.__deleteAllRequests()
        self.__deleteAllReportedValues()

    # --------------------------------------------------------------------------
    def requestValueOf( self, response_id ):
        self.list_of_requests[response_id]["calculateValue"] = True

    # --------------------------------------------------------------------------
    def requestGradientOf( self, response_id ):
        self.list_of_requests[response_id]["calculateGradient"] = True

    # --------------------------------------------------------------------------
    def isRequestingValueOf( self, response_id ):
        if response_id not in self.list_of_requests.keys():
            return False
        else:
            return self.list_of_requests[response_id]["calculateValue"]

    # --------------------------------------------------------------------------
    def isRequestingGradientOf( self, response_id ):
        if response_id not in self.list_of_requests.keys():
            return False
        else:
            return self.list_of_requests[response_id]["calculateGradient"]

    # --------------------------------------------------------------------------
    def reportValue( self, response_id, value ):
        self.__setFirstReportedValueAsReferenceIfSpecified( response_id, value )
        translated_value = self.__translateValueToStandardForm( response_id, value )
        self.__storeStandardizedValue( response_id, translated_value )

    # --------------------------------------------------------------------------
    def reportGradient( self, response_id, gradient ):
        translated_gradient = self.__translateGradientToStandardForm( response_id, gradient )
        self.__storeStandardizedGradient( response_id, translated_gradient )

    # --------------------------------------------------------------------------
    def getStandardizedValue( self, response_id ):
        return self.list_of_responses[response_id]["standardized_value"]

    # --------------------------------------------------------------------------
    def getReferenceValue( self, response_id ):
        return self.list_of_responses[response_id]["reference_value"]

    # --------------------------------------------------------------------------
    def getStandardizedGradient( self, response_id ):
        return self.list_of_responses[response_id]["standardized_gradient"]

    # --------------------------------------------------------------------------
    def __initializeListOfRequests( self ):
        self.list_of_requests = {}

        for objective_number in range( self.optimization_settings["objectives"].size()):
            objective_id = self.optimization_settings["objectives"][objective_number]["identifier"].GetString()
            self.list_of_requests[objective_id] = {"calculateValue": False, "calculateGradient": False}

        for constraint_number in range(self.optimization_settings["constraints"].size()):
            constraint_id = self.optimization_settings["constraints"][constraint_number]["identifier"].GetString()
            self.list_of_requests[constraint_id] = {"calculateValue": False, "calculateGradient": False}

    # --------------------------------------------------------------------------
    def __initializeListOfResponses( self ):
        self.list_of_responses = {}

        self.__addObjectivesToListOfResponses()
        self.__addConstraintsToListOfResponses()

    # --------------------------------------------------------------------------
    def __addObjectivesToListOfResponses( self ):
        for objective_number in range(self.optimization_settings["objectives"].size()):
            objective = self.optimization_settings["objectives"][objective_number]
            objective_id =  objective["identifier"].GetString()
            objective_task = objective["task"].GetString()

            if objective_task == "minimize":
                self.list_of_responses[objective_id] = { "type"                 : "minimize_objective",
                                                         "standardized_value"   : None,
                                                         "standardized_gradient": None }
            elif objective_task == "maximize":
                self.list_of_responses[objective_id] = { "type"                 : "maximize_objective",
                                                         "standardized_value"   : None,
                                                         "standardized_gradient": None}
            else:
                raise RuntimeError(">\nImproper definition of task for the following objective function: " + objective_id )

    # --------------------------------------------------------------------------
    def __addConstraintsToListOfResponses( self ):
        for constraint_number in range( self.optimization_settings["constraints"].size()):
            constraint = self.optimization_settings["constraints"][constraint_number]
            constraint_id = self.optimization_settings["constraints"][constraint_number]["identifier"].GetString()
            constraint_reference = constraint["reference"].GetString()

            if  constraint_reference == "specified_value":
                self.list_of_responses[constraint_id] = { "type"                 : "constraint",
                                                          "standardized_value"   : None,
                                                          "standardized_gradient": None,
                                                          "reference_value"      : constraint["reference_value"].GetDouble() }
            elif constraint_reference == "initial_value":
                self.list_of_responses[constraint_id] = {"type"                 : "constraint",
                                                         "standardized_value"   : None,
                                                         "standardized_gradient": None,
                                                         "reference_value"      : "waiting_for_initial_value" }
            else:
                raise RuntimeError(">\nImproper definition of reference for the following constraint function: " + constraint_id )

    # --------------------------------------------------------------------------
    def __deleteAllRequests( self ):
        for response_id in self.list_of_requests:
            self.list_of_requests[response_id]["calculateValue"] = False
            self.list_of_requests[response_id]["calculateGradient"] = False

    # --------------------------------------------------------------------------
    def __deleteAllReportedValues( self ):
        for response_id in self.list_of_responses:
            self.list_of_responses[response_id]["standardized_value"] = None
            self.list_of_responses[response_id]["standardized_gradient"] = None

    # --------------------------------------------------------------------------
    def __setFirstReportedValueAsReferenceIfSpecified( self, response_id, value ):
        response = self.list_of_responses[response_id]
        if response["type"] == "constraint" and response["reference_value"] == "waiting_for_initial_value":
            response["reference_value"] = value

    # --------------------------------------------------------------------------
    def __translateValueToStandardForm( self, response_id, value ):
        response = self.list_of_responses[response_id]
        if response["type"] == "constraint":
            value = value - response["reference_value"]
            return value
        else:
            return value

    # --------------------------------------------------------------------------
    def __translateGradientToStandardForm( self, response_id, gradient ):
        response = self.list_of_responses[response_id]
        if response["type"] == "maximize_objective":
            for local_gradient in gradient.values():
                local_gradient[0] = -local_gradient[0]
                local_gradient[1] = -local_gradient[1]
                local_gradient[2] = -local_gradient[2]
            return gradient
        else:
            return gradient

    # --------------------------------------------------------------------------
    def __storeStandardizedValue( self, response_id, value ):
        self.list_of_responses[response_id]["standardized_value"] = value

    # --------------------------------------------------------------------------
    def __storeStandardizedGradient( self, response_id, gradient ):
        self.list_of_responses[response_id]["standardized_gradient"] = gradient

# ==============================================================================
