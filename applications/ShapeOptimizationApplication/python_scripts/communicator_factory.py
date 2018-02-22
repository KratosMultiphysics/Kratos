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
            objective_id =  optimization_settings["objectives"][objective_number]["identifier"].GetString()
            self.list_of_responses[objective_id] = {"value": None, "referenceValue": None, "gradient": None}

        for constraint_number in range( optimization_settings["constraints"].size()):
            constraint = optimization_settings["constraints"][constraint_number]
            constraint_id = optimization_settings["constraints"][constraint_number]["identifier"].GetString()
            reference_value = None
            if  constraint["reference"].GetString() == "specified_value":
                reference_value =  constraint["reference_value"].GetDouble()
            elif constraint["reference"].GetString() == "initial_value":
                pass
            else:
                raise RuntimeError(">\nImproper definition of reference for the following response function: " + constraint_id )
            self.list_of_responses[constraint_id] = {"value": None, "referenceValue": reference_value, "gradient": None}

    # --------------------------------------------------------------------------
    def initializeCommunication( self ):
        self.__deleteAllRequests()
        self.__deleteAllReportedValues()

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
    def requestFunctionValueOf( self, response_id ):
        self.list_of_requests[response_id]["calculateValue"] = True

    # --------------------------------------------------------------------------
    def requestGradientOf( self, response_id ):
        self.list_of_requests[response_id]["calculateGradient"] = True

    # --------------------------------------------------------------------------
    def isRequestingFunctionValueOf( self, response_id ):
        if response_id not in self.list_of_requests.keys(): return False
        return self.list_of_requests[response_id]["calculateValue"]

    # --------------------------------------------------------------------------
    def isRequestingGradientOf( self, response_id ):
        if response_id not in self.list_of_requests.keys(): return False
        return self.list_of_requests[response_id]["calculateGradient"]

    # --------------------------------------------------------------------------
    def reportFunctionValue( self, response_id, functionValue ):
        if response_id in self.list_of_responses.keys():
            self.list_of_responses[response_id]["value"] = functionValue
        else:
            raise NameError("Reported function is not specified: " + response_id)

    # --------------------------------------------------------------------------
    def reportGradient( self, response_id, gradient ):
        if response_id in self.list_of_responses.keys():
            self.list_of_responses[response_id]["gradient"] = gradient
        else:
            raise NameError("Reported function is not specified: " + response_id)

    # --------------------------------------------------------------------------
    def setReferenceValue( self, response_id, reference_value ):
        if response_id in self.list_of_responses.keys():
            self.list_of_responses[response_id]["referenceValue"] = reference_value
        else:
            raise NameError("Function for which reference values is reported is not specified: " + response_id)

    # --------------------------------------------------------------------------
    def getFunctionValue( self, response_id ):
        return self.list_of_responses[response_id]["value"]

    # --------------------------------------------------------------------------
    def getReferenceValue( self, response_id ):
        if self.list_of_responses[response_id]["referenceValue"] is not None:
            return self.list_of_responses[response_id]["referenceValue"]
        else:
            raise NameError("Reference value required but not yet reported for the following response function: " + response_id)

    # --------------------------------------------------------------------------
    def getGradient( self, response_id ):
        return self.list_of_responses[response_id]["gradient"]

# ==============================================================================
