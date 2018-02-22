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
        self.__initializeListOfRequests( optimization_settings )
        self.__initializeListOfResponses( optimization_settings )

    # --------------------------------------------------------------------------
    def __initializeListOfRequests( self, optimization_settings ):
        self.list_of_requests = {}

        for objective_number in range( optimization_settings["objectives"].size()):
            objective_id = optimization_settings["objectives"][objective_number]["identifier"].GetString()
            self.list_of_requests[objective_id] = {"calculateValue": False, "calculateGradient": False}

        for constraint_number in range(optimization_settings["constraints"].size()):
            constraint_id = optimization_settings["constraints"][constraint_number]["identifier"].GetString()
            self.list_of_requests[constraint_id] = {"calculateValue": False, "calculateGradient": False}

    # --------------------------------------------------------------------------
    def __initializeListOfResponses( self, optimization_settings ):
        self.list_of_responses = {}

        for objective_number in range(optimization_settings["objectives"].size()):
            objective = optimization_settings["objectives"][objective_number]
            objective_id =  objective["identifier"].GetString()
            objective_task = objective["task"].GetString()

            if objective_task == "minimize":
                self.list_of_responses[objective_id] = {"type": "min_objective", "value": None, "gradient": None}
            elif objective_task == "maximize":
                self.list_of_responses[objective_id] = {"type": "max_objective", "value": None, "gradient": None}
            else:
                raise RuntimeError(">\nImproper definition of task for the following objective function: " + objective_id )

        for constraint_number in range( optimization_settings["constraints"].size()):
            constraint = optimization_settings["constraints"][constraint_number]
            constraint_id = optimization_settings["constraints"][constraint_number]["identifier"].GetString()
            constraint_reference = constraint["reference"].GetString()

            if  constraint_reference == "specified_value":
                reference_value =  constraint["reference_value"].GetDouble()
                self.list_of_responses[constraint_id] = {"type": "constraint", "value": None, "gradient": None, "referenceValue": reference_value}
            elif constraint_reference == "initial_value":
                self.list_of_responses[constraint_id] = {"type": "constraint", "value": None, "gradient": None, "referenceValue": "waitingForInitialValue"}
            else:
                raise RuntimeError(">\nImproper definition of reference for the following constraint function: " + constraint_id )

    # --------------------------------------------------------------------------
    def initializeCommunication( self ):
        self.__deleteAllRequests()
        self.__deleteAllReportedValues()

    # --------------------------------------------------------------------------
    def requestValue( self, response_id ):
        self.list_of_requests[response_id]["calculateValue"] = True

    # --------------------------------------------------------------------------
    def requestGradient( self, response_id ):
        self.list_of_requests[response_id]["calculateGradient"] = True

    # --------------------------------------------------------------------------
    def isRequestingValueOf( self, response_id ):
        if response_id not in self.list_of_requests.keys():
            return False
        return self.list_of_requests[response_id]["calculateValue"]

    # --------------------------------------------------------------------------
    def isRequestingGradientOf( self, response_id ):
        if response_id not in self.list_of_requests.keys():
            return False
        return self.list_of_requests[response_id]["calculateGradient"]

    # --------------------------------------------------------------------------
    def reportValue( self, response_id, value ):
        self.__setReferenceValueIfNecessary( response_id, value )
        self.__storeResponseValue( response_id, self.__translateValueToStandardForm( response_id, value ) )

    # --------------------------------------------------------------------------
    def reportGradient( self, response_id, gradient ):
        print("gradient[1001] = ", gradient[1001])
        translated_gradient = self.__translateGradientToStandardForm( response_id, gradient )
        print("translated_gradient[1001] = ", translated_gradient[1001])
        self.__storeResponseGradient( response_id, translated_gradient )

    # --------------------------------------------------------------------------
    def getValueInStandardForm( self, response_id ):
        return self.list_of_responses[response_id]["value"]

    # --------------------------------------------------------------------------
    def getReferenceValue( self, response_id ):
        return self.list_of_responses[response_id]["referenceValue"]

    # --------------------------------------------------------------------------
    def getGradientInStandardForm( self, response_id ):
        return self.list_of_responses[response_id]["gradient"]

    # --------------------------------------------------------------------------
    def __deleteAllRequests( self ):
        for response_id in self.list_of_requests:
            self.list_of_requests[response_id]["calculateValue"] = False
            self.list_of_requests[response_id]["calculateGradient"] = False

    # --------------------------------------------------------------------------
    def __deleteAllReportedValues( self ):
        for response_id in self.list_of_responses:
            self.list_of_responses[response_id]["value"] = None
            self.list_of_responses[response_id]["gradient"] = None

    # --------------------------------------------------------------------------
    def __setReferenceValueIfNecessary( self, response_id, value ):
        response = self.list_of_responses[response_id]
        if response["type"] == "constraint" and response["referenceValue"] == "waitingForInitialValue":
            response["referenceValue"] = value

    # --------------------------------------------------------------------------
    def __translateValueToStandardForm( self, response_id, value ):
        response = self.list_of_responses[response_id]
        if response["type"] == "constraint":
            value = value - response["referenceValue"]
            return value
        else:
            return value

    # --------------------------------------------------------------------------
    def __translateGradientToStandardForm( self, response_id, gradient ):
        response = self.list_of_responses[response_id]
        if response["type"] == "max_objective":
            gradient = -1*gradient
            return gradient
        else:
            return gradient

    # --------------------------------------------------------------------------
    def __storeResponseValue( self, response_id, value ):
        self.list_of_responses[response_id]["value"] = value

    # --------------------------------------------------------------------------
    def __storeResponseGradient( self, response_id, gradient ):
        self.list_of_responses[response_id]["gradient"] = gradient

# ==============================================================================
