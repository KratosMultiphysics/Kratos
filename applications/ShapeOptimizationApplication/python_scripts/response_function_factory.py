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
import structural_response_function_factory
import time as timer

# ==============================================================================
def CreateListOfResponseFunctions( optimization_settings, optimization_model_part ):
    listOfResponseFunctions = {}
    responseCreator = ResponseFunctionCreator( optimization_settings, optimization_model_part )
    responseCreator.AddSpecifiedKratosResponseFunctionsToList( listOfResponseFunctions )
    return listOfResponseFunctions

# ==============================================================================
class ResponseFunctionCreator:
    # --------------------------------------------------------------------------
    def __init__( self, optimization_settings, optimization_model_part ):
        self.optimization_settings = optimization_settings
        self.optimization_model_part = optimization_model_part

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
        response_type = response_settings["response_type"].GetString()
        if response_type in ["strain_energy", "mass", "eigenfrequency", "adjoint_nodal_displacement", "adjoint_strain_energy"]:
            self.listOfResponseFunctions[response_id] = structural_response_function_factory.CreateResponseFunction(response_id, response_settings, self.optimization_model_part)
        else:
            raise NameError("The following response function is available for optimization: " + response_id)

# ==============================================================================
