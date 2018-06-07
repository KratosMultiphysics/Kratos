# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya, https://github.com/sunethwarna
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
import response_function_factory
import time as timer
import os
import shutil
import sys

# ==============================================================================
class KratosExternalTransientAnalyzer( (__import__("analyzer_base")).AnalyzerBaseClass ):
    # --------------------------------------------------------------------------
    def __init__( self, project_parameters, model_part ):
        self.output_directory = project_parameters["optimization_settings"]["output"]["output_directory"].GetString()
        sys.path.append(os.path.dirname(os.path.realpath(__file__)) + os.path.sep + "base_case")
        self.response_function_list = response_function_factory.CreateListOfResponseFunctions(project_parameters["optimization_settings"], model_part)

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop( self ):
        for response in self.response_function_list.values():
            response.Initialize()
    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):

        # creating the folder structure
        if not os.path.isdir("./%s/iterations" % self.output_directory):
            os.mkdir("./%s/iterations" % self.output_directory)
        
        current_iteration_dir = "./%s/iterations/%d" % (self.output_directory, optimizationIteration)
        
        if not os.path.isdir(current_iteration_dir):
            print("> Creating simulation iteration %d data..." % optimizationIteration)
            shutil.copytree("./base_case/", "%s/" % current_iteration_dir)

        # set the standard output to a file, redirecting both python and c++ outputs to local directory

        # changing the current working directory to iteration directory
        _main_working_dir = os.getcwd()
        os.chdir(current_iteration_dir)

        # writing the updated model part to the iteration directory

        for identifier, response in self.response_function_list.items():

            response.InitializeSolutionStep()

            # response values
            if communicator.isRequestingValueOf(identifier):
                response.CalculateValue()
                communicator.reportValue(identifier, response.GetValue())

            # response gradients
            if communicator.isRequestingGradientOf(identifier):
                communicator.reportGradient(identifier, response.GetShapeGradient())

            response.FinalizeSolutionStep()
        
        os.chdir(_main_working_dir)

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop( self ):
        for response in self.response_function_list.values():
            response.Finalize()
        
# ==============================================================================
