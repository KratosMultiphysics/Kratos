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

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------
# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# ==============================================================================
def CreateSolver( inputModelPart, optimizationSettings ):
    listOfResponseFunctionSolver = {}
    ResponseFunctionSolverCreator( inputModelPart, optimizationSettings ).AddSpecifiedKratosResponseFunctionSolverToList( listOfResponseFunctionSolver )
    return listOfResponseFunctionSolver

# ==============================================================================
class ResponseFunctionSolverCreator: 

    # --------------------------------------------------------------------------
    def __init__( self, inputModelPart, optimizationSettings ):
        self.inputModelPart = inputModelPart
        self.optimizationSettings = optimizationSettings

     # --------------------------------------------------------------------------
    def AddSpecifiedKratosResponseFunctionSolverToList( self, listOfResponseFunctionSolver ):        
        self.listOfResponseFunctionSolver = listOfResponseFunctionSolver
        self.addObjectivesToListOfResponseFunctionSolver()
        self.addConstraintsToListOfResponseFunctionSolver()
        
    # --------------------------------------------------------------------------
    def addObjectivesToListOfResponseFunctionSolver( self ):

        numberOfObjectives = self.optimizationSettings["objectives"].size()

        for objectiveNumber in range(numberOfObjectives):

            objectiveId = self.optimizationSettings["objectives"][objectiveNumber]["identifier"].GetString()
            useKratosSolver = self.optimizationSettings["objectives"][objectiveNumber]["use_kratos_solver"].GetBool()

            if useKratosSolver:
                self.checkIfGivenResponseFunctionIsAlreadyDefined( objectiveId )
                self.createAndAddSolverForGivenResponse( objectiveId, self.optimizationSettings["objectives"][objectiveNumber] )

        if not self.listOfResponseFunctionSolver:
            raise ValueError("No objective function specified!")

    # --------------------------------------------------------------------------
    def addConstraintsToListOfResponseFunctionSolver( self ):

        numberOfConstraints = self.optimizationSettings["constraints"].size()

        for constraintNumber in range(numberOfConstraints):

            constraintId = self.optimizationSettings["constraints"][constraintNumber]["identifier"].GetString()
            useKratosSolver = self.optimizationSettings["constraints"][constraintNumber]["use_kratos_solver"].GetBool()

            if useKratosSolver:
                self.checkIfGivenResponseFunctionIsAlreadyDefined( constraintId )
                self.createAndAddSolverForGivenResponse( constraintId, self.optimizationSettings["constraints"][constraintNumber] )         

    # --------------------------------------------------------------------------
    def checkIfGivenResponseFunctionIsAlreadyDefined( self, responseId ):
        if responseId in self.listOfResponseFunctionSolver.keys():
            raise NameError("There are multiple response functions with the following identifier: " + responseId)

    # --------------------------------------------------------------------------
    def createAndAddSolverForGivenResponse( self, responseId, solverSettings ):

        responseFunctionSolverIsNotImplemented = True

        if responseId == "strain_energy":
            responseFunctionSolverIsNotImplemented = False
            self.inputModelPart.AddNodalSolutionStepVariable(STRAIN_ENERGY_SHAPE_GRADIENT)
            self.listOfResponseFunctionSolver["strain_energy"] = StrainEnergyResponseFunction( self.inputModelPart, solverSettings )

        if responseId == "mass":
            responseFunctionSolverIsNotImplemented = False
            self.inputModelPart.AddNodalSolutionStepVariable(MASS_SHAPE_GRADIENT)
            self.listOfResponseFunctionSolver["mass"] = MassResponseFunction( self.inputModelPart, solverSettings )   

        if responseFunctionSolverIsNotImplemented:
            raise NameError("The following response function is not specified: " + responseId)

# ==============================================================================
