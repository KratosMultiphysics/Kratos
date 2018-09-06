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
from KratosMultiphysics.StructuralMechanicsApplication import *

# Additional imports
import response_function_factory
import time as timer

# ==============================================================================
class KratosInternalAnalyzer( (__import__("analyzer_base")).AnalyzerBaseClass ):
    # --------------------------------------------------------------------------
    def __init__( self, optimization_settings, model_part_controller ):
        self.model_part_controller = model_part_controller

        analysis_mdpa = model_part_controller.GetOptimizationModelPart()
        self.response_function_list = response_function_factory.CreateListOfResponseFunctions(optimization_settings, analysis_mdpa)

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop( self ):
        for response in self.response_function_list.values():
            response.Initialize()
    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):

        for identifier, response in self.response_function_list.items():

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

            self.__ClearResultsFromModelPart()

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop( self ):
        for response in self.response_function_list.values():
            response.Finalize()

    # --------------------------------------------------------------------------
    def __ClearResultsFromModelPart( self ):
        self.model_part_controller.SetMeshToReferenceMesh()
        self.model_part_controller.SetDeformationVariablesToZero()

# ==============================================================================
