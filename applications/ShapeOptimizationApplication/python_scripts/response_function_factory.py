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
def CreateListOfResponseFunctions( optimization_model_part, optimization_settings ):
    listOfResponseFunctions = {}
    responseCreator = ResponseFunctionCreator( optimization_model_part, optimization_settings )
    responseCreator.AddSpecifiedKratosResponseFunctionsToList( listOfResponseFunctions )
    return listOfResponseFunctions

# ==============================================================================
class ResponseFunctionCreator:

    # --------------------------------------------------------------------------
    def __init__( self, optimization_model_part, optimization_settings ):
        self.optimization_model_part = optimization_model_part
        self.optimization_settings = optimization_settings

     # --------------------------------------------------------------------------
    def AddSpecifiedKratosResponseFunctionsToList( self, listOfResponseFunctions ):
        self.listOfResponseFunctions = listOfResponseFunctions
        self.__addObjectivesToListOfResponseFunctions()
        self.__addConstraintsToListOfResponseFunctions()

    # --------------------------------------------------------------------------
    def __addObjectivesToListOfResponseFunctions( self ):
        for objective_number in range(self.optimization_settings["objectives"].size()):
            objective = self.optimization_settings["objectives"][objective_number]
            objective_id = objective["identifier"].GetString()
            if objective["use_kratos"].GetBool():
                self.__checkIfGivenResponseFunctionIsAlreadyDefined( objective_id )
                self.__createAndAddGivenResponse( objective_id, objective["kratos_response_settings"] )

        if not self.listOfResponseFunctions:
            raise ValueError("No objective function specified!")

    # --------------------------------------------------------------------------
    def __addConstraintsToListOfResponseFunctions( self ):
        for constraint_number in range(self.optimization_settings["constraints"].size()):
            constraint = self.optimization_settings["constraints"][constraint_number]
            constraint_id = constraint["identifier"].GetString()
            if constraint["use_kratos"].GetBool():
                self.__checkIfGivenResponseFunctionIsAlreadyDefined( constraint_id )
                self.__createAndAddGivenResponse( constraint_id, constraint["kratos_response_settings"] )

    # --------------------------------------------------------------------------
    def __checkIfGivenResponseFunctionIsAlreadyDefined( self, response_id ):
        if response_id in self.listOfResponseFunctions.keys():
            raise NameError("There are multiple response functions with the following identifier: " + response_id)

    # --------------------------------------------------------------------------
    def __createAndAddGivenResponse( self, response_id, response_settings ):
        if response_id == "strain_energy":
            self.optimization_model_part.AddNodalSolutionStepVariable(STRAIN_ENERGY_SHAPE_GRADIENT)
            self.listOfResponseFunctions["strain_energy"] = StrainEnergyResponseFunction( self.optimization_model_part, response_settings )

        elif response_id == "mass":
            self.optimization_model_part.AddNodalSolutionStepVariable(MASS_SHAPE_GRADIENT)
            self.listOfResponseFunctions["mass"] = MassResponseFunction( self.optimization_model_part, response_settings )

        elif response_id == "eigenfrequency":
            self.optimization_model_part.AddNodalSolutionStepVariable(EIGENFREQUENCY_SHAPE_GRADIENT)
            if not response_settings.Has("weighting_method") or response_settings["weighting_method"].GetString() == "none":
                self.listOfResponseFunctions["eigenfrequency"] = EigenfrequencyResponseFunction( self.optimization_model_part, response_settings )

            elif response_settings["weighting_method"].GetString() == "linear_scaling":
                self.listOfResponseFunctions["eigenfrequency"] = EigenfrequencyResponseFunctionLinScal( self.optimization_model_part, response_settings )

            else:
                raise NameError("The following weighting_method is not valid for eigenfrequency response: " + response_settings["weighting_method"].GetString())
        else:
            raise NameError("The following response function is not specified: " + response_id)

# ==============================================================================
