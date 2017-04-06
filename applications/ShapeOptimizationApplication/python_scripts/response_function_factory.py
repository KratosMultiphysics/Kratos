# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
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
def CreateListOfResponseFunctions( inputModelPart, optimizationSettings ):
    listOfResponseFunctions = {}
    responseCreator = ResponseFunctionCreator( inputModelPart, optimizationSettings )
    responseCreator.AddSpecifiedKratosResponseFunctionsToList( listOfResponseFunctions )
    return listOfResponseFunctions

# ==============================================================================
class ResponseFunctionCreator: 

    # --------------------------------------------------------------------------
    def __init__( self, inputModelPart, optimizationSettings ):
        self.inputModelPart = inputModelPart
        self.optimizationSettings = optimizationSettings

     # --------------------------------------------------------------------------
    def AddSpecifiedKratosResponseFunctionsToList( self, listOfResponseFunctions ):        
        self.listOfResponseFunctions = listOfResponseFunctions
        self.addObjectivesToListOfResponseFunctions()
        self.addConstraintsToListOfResponseFunctions()
        
    # --------------------------------------------------------------------------
    def addObjectivesToListOfResponseFunctions( self ):

        numberOfObjectives = self.optimizationSettings["objectives"].size()

        for objectiveNumber in range(numberOfObjectives):

            objectiveId = self.optimizationSettings["objectives"][objectiveNumber]["identifier"].GetString()
            useKratos = self.optimizationSettings["objectives"][objectiveNumber]["use_kratos"].GetBool()

            if useKratos:
                self.checkIfGivenResponseFunctionIsAlreadyDefined( objectiveId )
                self.createAndAddGivenResponse( objectiveId, self.optimizationSettings["objectives"][objectiveNumber] )

        if not self.listOfResponseFunctions:
            raise ValueError("No objective function specified!")

    # --------------------------------------------------------------------------
    def addConstraintsToListOfResponseFunctions( self ):

        numberOfConstraints = self.optimizationSettings["constraints"].size()

        for constraintNumber in range(numberOfConstraints):

            constraintId = self.optimizationSettings["constraints"][constraintNumber]["identifier"].GetString()
            useKratos = self.optimizationSettings["constraints"][constraintNumber]["use_kratos"].GetBool()

            if useKratos:
                self.checkIfGivenResponseFunctionIsAlreadyDefined( constraintId )
                self.createAndAddGivenResponse( constraintId, self.optimizationSettings["constraints"][constraintNumber] )         

    # --------------------------------------------------------------------------
    def checkIfGivenResponseFunctionIsAlreadyDefined( self, responseId ):
        if responseId in self.listOfResponseFunctions.keys():
            raise NameError("There are multiple response functions with the following identifier: " + responseId)

    # --------------------------------------------------------------------------
    def createAndAddGivenResponse( self, responseId, solverSettings ):

        if responseId == "strain_energy":
            responseFunctionSolverIsNotImplemented = False
            self.inputModelPart.AddNodalSolutionStepVariable(STRAIN_ENERGY_SHAPE_GRADIENT)
            self.listOfResponseFunctions["strain_energy"] = StrainEnergyResponseFunction( self.inputModelPart, solverSettings )
        elif responseId == "mass":
            responseFunctionSolverIsNotImplemented = False
            self.inputModelPart.AddNodalSolutionStepVariable(MASS_SHAPE_GRADIENT)
            self.listOfResponseFunctions["mass"] = MassResponseFunction( self.inputModelPart, solverSettings )   
        else:
            raise NameError("The following response function is not specified: " + responseId)

# ==============================================================================
