# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
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
def CreateListOfResponseFunctions( OptimizationModelPart, OptimizationSettings ):
    listOfResponseFunctions = {}
    responseCreator = ResponseFunctionCreator( OptimizationModelPart, OptimizationSettings )
    responseCreator.AddSpecifiedKratosResponseFunctionsToList( listOfResponseFunctions )
    return listOfResponseFunctions

# ==============================================================================
class ResponseFunctionCreator: 

    # --------------------------------------------------------------------------
    def __init__( self, OptimizationModelPart, OptimizationSettings ):
        self.OptimizationModelPart = OptimizationModelPart
        self.OptimizationSettings = OptimizationSettings

     # --------------------------------------------------------------------------
    def AddSpecifiedKratosResponseFunctionsToList( self, listOfResponseFunctions ):        
        self.listOfResponseFunctions = listOfResponseFunctions
        self.__addObjectivesToListOfResponseFunctions()
        self.__addConstraintsToListOfResponseFunctions()
        
    # --------------------------------------------------------------------------
    def __addObjectivesToListOfResponseFunctions( self ):

        numberOfObjectives = self.OptimizationSettings["objectives"].size()

        for objectiveNumber in range(numberOfObjectives):

            objectiveId = self.OptimizationSettings["objectives"][objectiveNumber]["identifier"].GetString()
            useKratos = self.OptimizationSettings["objectives"][objectiveNumber]["use_kratos"].GetBool()

            if useKratos:
                self.__checkIfGivenResponseFunctionIsAlreadyDefined( objectiveId )
                self.__createAndAddGivenResponse( objectiveId, self.OptimizationSettings["objectives"][objectiveNumber] )

        if not self.listOfResponseFunctions:
            raise ValueError("No objective function specified!")

    # --------------------------------------------------------------------------
    def __addConstraintsToListOfResponseFunctions( self ):

        numberOfConstraints = self.OptimizationSettings["constraints"].size()

        for constraintNumber in range(numberOfConstraints):

            constraintId = self.OptimizationSettings["constraints"][constraintNumber]["identifier"].GetString()
            useKratos = self.OptimizationSettings["constraints"][constraintNumber]["use_kratos"].GetBool()

            if useKratos:
                self.__checkIfGivenResponseFunctionIsAlreadyDefined( constraintId )
                self.__createAndAddGivenResponse( constraintId, self.OptimizationSettings["constraints"][constraintNumber] )         

    # --------------------------------------------------------------------------
    def __checkIfGivenResponseFunctionIsAlreadyDefined( self, responseId ):
        if responseId in self.listOfResponseFunctions.keys():
            raise NameError("There are multiple response functions with the following identifier: " + responseId)

    # --------------------------------------------------------------------------
    def __createAndAddGivenResponse( self, responseId, solverSettings ):

        if responseId == "strain_energy":
            responseFunctionSolverIsNotImplemented = False
            self.OptimizationModelPart.AddNodalSolutionStepVariable(STRAIN_ENERGY_SHAPE_GRADIENT)
            self.listOfResponseFunctions["strain_energy"] = StrainEnergyResponseFunction( self.OptimizationModelPart, solverSettings )
        elif responseId == "mass":
            responseFunctionSolverIsNotImplemented = False
            self.OptimizationModelPart.AddNodalSolutionStepVariable(MASS_SHAPE_GRADIENT)
            self.listOfResponseFunctions["mass"] = MassResponseFunction( self.OptimizationModelPart, solverSettings )   
        else:
            raise NameError("The following response function is not specified: " + responseId)

# ==============================================================================
