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

# Kratos Core and Apps
import KratosMultiphysics as km

# Additional imports
from KratosMultiphysics.StructuralMechanicsApplication import structural_response_function_factory as csm_response_factory
import time as timer

# ==============================================================================
class KratosInternalAnalyzer( (__import__("analyzer_base")).AnalyzerBaseClass ):
    # --------------------------------------------------------------------------
    def __init__( self, specified_responses, model_part_controller ):
        self.model_part_controller = model_part_controller
        self.response_functions = self.__CreateResponseFunctions(specified_responses, model_part_controller.GetModel())

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop( self ):
        for _, response in self.response_functions:
            response.Initialize()

    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):
        optimization_model_part = self.model_part_controller.GetOptimizationModelPart()

        time_before_analysis = optimization_model_part.ProcessInfo.GetValue(km.TIME)
        step_before_analysis = optimization_model_part.ProcessInfo.GetValue(km.STEP)
        delta_time_before_analysis = optimization_model_part.ProcessInfo.GetValue(km.DELTA_TIME)

        for identifier, response in self.response_functions:

            # Reset step/time iterators such that they match the optimization iteration after calling CalculateValue (which internally calls CloneTimeStep)
            optimization_model_part.ProcessInfo.SetValue(km.STEP, step_before_analysis-1)
            optimization_model_part.ProcessInfo.SetValue(km.TIME, time_before_analysis-1)
            optimization_model_part.ProcessInfo.SetValue(km.DELTA_TIME, 0)

            response.InitializeSolutionStep()

            # response values
            if communicator.isRequestingValueOf(identifier):
                response.CalculateValue()
                communicator.reportValue(identifier, response.GetValue())

            # response gradients
            if communicator.isRequestingGradientOf(identifier):
                response.CalculateGradient()
                communicator.reportGradient(identifier, response.GetShapeGradient())

            response.FinalizeSolutionStep()

            # Clear results or modifications on model part
            optimization_model_part.ProcessInfo.SetValue(km.STEP, step_before_analysis)
            optimization_model_part.ProcessInfo.SetValue(km.TIME, time_before_analysis)
            optimization_model_part.ProcessInfo.SetValue(km.DELTA_TIME, delta_time_before_analysis)

            self.model_part_controller.SetMeshToReferenceMesh()
            self.model_part_controller.SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop( self ):
        for _, response in self.response_functions:
            response.Finalize()

    # --------------------------------------------------------------------------
    @staticmethod
    def __CreateResponseFunctions( specified_responses, model ):
        available_csm_response_functions = ["strain_energy", "mass", "eigenfrequency"]

        response_functions = []
        response_ids = []

        for response_id, response_settings in specified_responses:
            if response_id in response_ids:
                raise NameError("There are multiple response functions with the following identifier: " + response_id)

            if response_settings["response_type"].GetString() not in available_csm_response_functions:
                raise NameError("The following structural response function is not available: " + response_id)

            response_ids.append(response_id)
            response = csm_response_factory.CreateResponseFunction(response_id, response_settings, model)
            response_functions.append((response_id, response))

        return response_functions

# ==============================================================================
