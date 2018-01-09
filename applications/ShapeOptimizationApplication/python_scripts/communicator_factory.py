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
def CreateCommunicator( optimizationSettings ):
    return Communicator( optimizationSettings )

# ==============================================================================
class Communicator:
    
    # --------------------------------------------------------------------------
    def __init__( self, optimizationSettings ):
        self.listOfRequests = self.__initializeListOfRequests( optimizationSettings )
        self.listOfResponses = self.__initializeListOfResponses( optimizationSettings )

    # --------------------------------------------------------------------------
    def __initializeListOfRequests( self, optimizationSettings ):

        listOfRequests = {}
        numberOfObjectives = optimizationSettings["objectives"].size()
        numberOfConstraints = optimizationSettings["constraints"].size()

        for objectiveNumber in range(numberOfObjectives):
            objectiveId = optimizationSettings["objectives"][objectiveNumber]["identifier"].GetString()
            listOfRequests[objectiveId] = {"calculateValue": False, "calculateGradient": False}

        for constraintNumber in range(numberOfConstraints):
            constraintId = optimizationSettings["constraints"][constraintNumber]["identifier"].GetString()
            listOfRequests[constraintId] = {"calculateValue": False, "calculateGradient": False}

        return listOfRequests

    # --------------------------------------------------------------------------
    def __initializeListOfResponses( self, optimizationSettings ):
        
        listOfResponses = {}  
        numberOfObjectives = optimizationSettings["objectives"].size()
        numberOfConstraints = optimizationSettings["constraints"].size()

        for objectiveNumber in range(numberOfObjectives):
            objectiveId = optimizationSettings["objectives"][objectiveNumber]["identifier"].GetString()
            listOfResponses[objectiveId] = {}

        for constraintNumber in range(numberOfConstraints):
            constraintId = optimizationSettings["constraints"][constraintNumber]["identifier"].GetString()
            listOfResponses[constraintId] = {}  

        for responseId in listOfResponses:
            listOfResponses[responseId] = {"value": None, "referenceValue": None, "gradient": None} 

        return listOfResponses
    
    # --------------------------------------------------------------------------
    def initializeCommunication( self ):
        self.__deleteAllRequests()
        self.__deleteAllReportedValues()

    # --------------------------------------------------------------------------
    def __deleteAllRequests( self ):
        for responseId in self.listOfRequests:
            self.listOfRequests[responseId]["CalculateValue"] = False
            self.listOfRequests[responseId]["calculateGradient"] = False

    # --------------------------------------------------------------------------
    def __deleteAllReportedValues( self ):
        for responseId in self.listOfResponses:
            self.listOfResponses[responseId]["value"] = None 
            self.listOfResponses[responseId]["gradient"] = None 
            
    # --------------------------------------------------------------------------
    def requestFunctionValueOf( self, responseId ):
        self.listOfRequests[responseId]["calculateValue"] = True

    # --------------------------------------------------------------------------
    def requestGradientOf( self, responseId ):
        self.listOfRequests[responseId]["calculateGradient"] = True

    # --------------------------------------------------------------------------
    def isRequestingFunctionValueOf( self, responseId ):
        if responseId not in self.listOfRequests.keys(): return False
        return self.listOfRequests[responseId]["calculateValue"]

    # --------------------------------------------------------------------------
    def isRequestingGradientOf( self, responseId ):
        if responseId not in self.listOfRequests.keys(): return False
        return self.listOfRequests[responseId]["calculateGradient"]

    # --------------------------------------------------------------------------
    def reportFunctionValue( self, responseId, functionValue ):
        if responseId in self.listOfResponses.keys():
            self.listOfResponses[responseId]["value"] = functionValue
        else:
            raise NameError("Reported function is not specified: " + responseId)

    # --------------------------------------------------------------------------
    def reportGradient( self, responseId, gradient ):
        if responseId in self.listOfResponses.keys():        
            self.listOfResponses[responseId]["gradient"] = gradient
        else:
            raise NameError("Reported function is not specified: " + responseId)

    # --------------------------------------------------------------------------
    def setFunctionReferenceValue( self, responseId, functionReferenceValue ):
        if responseId in self.listOfResponses.keys():
            self.listOfResponses[responseId]["referenceValue"] = functionReferenceValue
        else:
            raise NameError("Reported function is not specified: " + responseId)

    # --------------------------------------------------------------------------
    def getReportedFunctionValueOf( self, responseId ):
        return self.listOfResponses[responseId]["value"]

    # --------------------------------------------------------------------------
    def getReportedFunctionReferenceValueOf( self, responseId ):
        return self.listOfResponses[responseId]["referenceValue"]

    # --------------------------------------------------------------------------    
    def getReportedGradientOf( self, responseId ):
        return self.listOfResponses[responseId]["gradient"]

    # --------------------------------------------------------------------------

# ==============================================================================
