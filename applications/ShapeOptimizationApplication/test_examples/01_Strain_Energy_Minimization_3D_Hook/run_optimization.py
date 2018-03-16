from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# For time measures
import time as timer

# ======================================================================================================================================
# Model part & solver
# ======================================================================================================================================

#import define_output
parameter_file = open("optimization_project_parameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

#set echo level
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()

#defining the model_part
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

###TODO replace this "model" for real one once available in kratos core
Model = {ProjectParameters["problem_data"]["model_part_name"].GetString() : main_model_part}

# Create an optimizer
# Note that internally variables related to the optimizer are added to the model part
optimizerFactory = __import__("optimizer_factory")
optimizer = optimizerFactory.CreateOptimizer( main_model_part, ProjectParameters["optimization_settings"] )

# Create solver for all response functions specified in the optimization settings
# Note that internally variables related to the individual functions are added to the model part
responseFunctionFactory = __import__("response_function_factory")
response_function_list = responseFunctionFactory.CreateListOfResponseFunctions( main_model_part, ProjectParameters["optimization_settings"] )

# ======================================================================================================================================
# Analyzer
# ======================================================================================================================================

class kratosCSMAnalyzer( (__import__("analyzer_base")).analyzerBaseClass ):

    # --------------------------------------------------------------------------
    # NOTES:
    # Currently a separate modelpart is created for:
    #   - primal solver & optimizer (use the same)
    #   - adjoint solver
    # The nodal coordinates are synchronized at each analyzer call

    # --------------------------------------------------------------------------
    def initializeBeforeOptimizationLoop( self ):
        for response in response_function_list.values():
            response.Initialize()

    # --------------------------------------------------------------------------
    def analyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):

        # TODO initialize evaluation step

        # response values
        for identifier, response in response_function_list.items():
            if communicator.isRequestingValueOf(identifier):
                value = response.CalculateValue()
                communicator.reportValue(identifier, value)

        # add optional custom responses

        # response gradients
        for identifier, response in response_function_list.items():
            if communicator.isRequestingGradientOf(identifier):
                response.CalculateGradient()
                communicator.reportGradient(identifier, response.GetShapeGradient())

        # add optional custom gradients

        # TODO finalize evaluation step

    # --------------------------------------------------------------------------
    def finalizeAfterOptimizationLoop( self ):
        for response in response_function_list.values():
            response.Finalize()

    # --------------------------------------------------------------------------

structureAnalyzer = kratosCSMAnalyzer()

# ======================================================================================================================================
# Optimization
# ======================================================================================================================================

optimizer.importAnalyzer( structureAnalyzer )
optimizer.optimize()

# ======================================================================================================================================