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
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *

# Additional imports
import response_function_factory
import time as timer

# ==============================================================================
class KratosInternalAnalyzer( (__import__("analyzer_base")).AnalyzerBaseClass ):
    # --------------------------------------------------------------------------
    def __init__( self, optimization_settings, model_part_controller ):
        self.model_part_controller = model_part_controller

        self.response_function_list = response_function_factory.CreateListOfResponseFunctions(optimization_settings, self.model_part_controller.GetModel())

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop( self ):
        for response in self.response_function_list.values():
            response.Initialize()
    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):
        optimization_model_part = self.model_part_controller.GetOptimizationModelPart()

        time_before_analysis = optimization_model_part.ProcessInfo.GetValue(TIME)
        step_before_analysis = optimization_model_part.ProcessInfo.GetValue(STEP)
        delta_time_before_analysis = optimization_model_part.ProcessInfo.GetValue(DELTA_TIME)

        for identifier, response in self.response_function_list.items():

            # Reset step/time iterators such that they match the optimization iteration after calling CalculateValue (which internally calls CloneTimeStep)
            optimization_model_part.ProcessInfo.SetValue(STEP, step_before_analysis-1)
            optimization_model_part.ProcessInfo.SetValue(TIME, time_before_analysis-1)
            optimization_model_part.ProcessInfo.SetValue(DELTA_TIME, 0)

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
            optimization_model_part.ProcessInfo.SetValue(STEP, step_before_analysis)
            optimization_model_part.ProcessInfo.SetValue(TIME, time_before_analysis)
            optimization_model_part.ProcessInfo.SetValue(DELTA_TIME, delta_time_before_analysis)

            self.model_part_controller.SetMeshToReferenceMesh()
            self.model_part_controller.SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop( self ):
        for response in self.response_function_list.values():
            response.Finalize()

# ==============================================================================
