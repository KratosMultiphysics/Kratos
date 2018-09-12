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

# ==============================================================================
def CreateListOfResponseFunctions( optimization_settings, optimization_model_part ):
    list_of_response_functions = {}
    response_creator = ResponseFunctionCreator( optimization_settings, optimization_model_part )
    response_creator.AddSpecifiedKratosResponseFunctionsToList( list_of_response_functions )
    return list_of_response_functions

# ==============================================================================
class ResponseFunctionCreator:
    # --------------------------------------------------------------------------
    def __init__( self, optimization_settings, optimization_model_part ):
        self.optimization_settings = optimization_settings
        self.optimization_model_part = optimization_model_part

     # --------------------------------------------------------------------------
    def AddSpecifiedKratosResponseFunctionsToList( self, list_of_response_functions ):
        self.list_of_response_functions = list_of_response_functions
        self.__AddObjectivesToListOfResponseFunctions()
        self.__AddConstraintsToListOfResponseFunctions()

    # --------------------------------------------------------------------------
    def __AddObjectivesToListOfResponseFunctions( self ):
        for objective_number in range(self.optimization_settings["objectives"].size()):
            objective = self.optimization_settings["objectives"][objective_number]
            objective_id = objective["identifier"].GetString()
            if objective["use_kratos"].GetBool():
                self.__CheckIfGivenResponseFunctionIsAlreadyDefined( objective_id )
                self.__CreateAndAddGivenResponse( objective_id, objective["kratos_response_settings"] )

        if not self.list_of_response_functions:
            raise ValueError("No objective function specified!")

    # --------------------------------------------------------------------------
    def __AddConstraintsToListOfResponseFunctions( self ):
        for constraint_number in range(self.optimization_settings["constraints"].size()):
            constraint = self.optimization_settings["constraints"][constraint_number]
            constraint_id = constraint["identifier"].GetString()
            if constraint["use_kratos"].GetBool():
                self.__CheckIfGivenResponseFunctionIsAlreadyDefined( constraint_id )
                self.__CreateAndAddGivenResponse( constraint_id, constraint["kratos_response_settings"] )

    # --------------------------------------------------------------------------
    def __CheckIfGivenResponseFunctionIsAlreadyDefined( self, response_id ):
        if response_id in self.list_of_response_functions.keys():
            raise NameError("There are multiple response functions with the following identifier: " + response_id)

    # --------------------------------------------------------------------------
    def __CreateAndAddGivenResponse( self, response_id, response_settings ):
        response_type = response_settings["response_type"].GetString()
        if response_type in ["strain_energy", "mass", "eigenfrequency"]:
            model = Model()
            model.AddModelPart(self.optimization_model_part)
            self.list_of_response_functions[response_id] = structural_response_function_factory.CreateResponseFunction(response_id, response_settings, model)
        else:
            raise NameError("The following response function is not available for optimization: " + response_id)

# ==============================================================================
