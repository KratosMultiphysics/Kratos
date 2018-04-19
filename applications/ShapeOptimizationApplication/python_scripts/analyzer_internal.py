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
    def __init__( self, project_parameters, model_part ):
        self.response_function_list = response_function_factory.CreateListOfResponseFunctions(project_parameters["optimization_settings"], model_part)

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop( self ):
        for response in self.response_function_list.values():
            response.Initialize()
    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):

        for identifier, response in self.response_function_list.items():
            # response values
            if communicator.isRequestingValueOf(identifier):
                startTime = timer.time()
                print("> Calculating response value of '" + identifier + "'...")
                value = response.CalculateValue()
                communicator.reportValue(identifier, value)
                print("> Time needed for calculating response value of '" + identifier + "' = ",round(timer.time() - startTime,2),"s")

            # response gradients
            if communicator.isRequestingGradientOf(identifier):
                startTime = timer.time()
                print("> Calculating response gradient of '" + identifier + "'...")
                response.CalculateGradient()
                communicator.reportGradient(identifier, response.GetShapeGradient())
                print("> Time needed for calculating response gradient of '" + identifier + "' = ",round(timer.time() - startTime,2),"s")

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop( self ):
        for response in self.response_function_list.values():
            response.Finalize()

# ==============================================================================
