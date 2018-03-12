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

# Read parameters
with open("ProjectParameters.json",'r') as parameter_file:
    ProjectParameters = Parameters(parameter_file.read())

# Defining the model_part
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

# Create an optimizer
# Note that internally variables related to the optimizer are added to the model part
optimizerFactory = __import__("optimizer_factory")
optimizer = optimizerFactory.CreateOptimizer(ProjectParameters["optimization_settings"], main_model_part)

# Create solver for all response functions specified in the optimization settings
# Note that internally variables related to the individual functions are added to the model part
responseFunctionFactory = __import__("response_function_factory")
listOfResponseFunctions = responseFunctionFactory.CreateListOfResponseFunctions(ProjectParameters["optimization_settings"], main_model_part)

# Create structural solver
# Note that internally variables related to the individual functions are added to the model part
csm_analysis = __import__("structural_mechanics_analysis").StructuralMechanicsAnalysis(ProjectParameters, main_model_part)

# ======================================================================================================================================
# Analyzer
# ======================================================================================================================================

class kratosCSMAnalyzer( (__import__("analyzer_base")).analyzerBaseClass ):

    # --------------------------------------------------------------------------
    def initializeBeforeOptimizationLoop( self ):
        csm_analysis.Initialize()

    # --------------------------------------------------------------------------
    def analyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):

        # Calculation of eigenfrequency
        if communicator.isRequestingValueOf("eigenfrequency"):

            print("\n> Starting StructuralMechanicsApplication to solve structure")
            startTime = timer.time()
            csm_analysis.InitializeTimeStep()
            csm_analysis.SolveTimeStep()
            csm_analysis.FinalizeTimeStep()
            print("> Time needed for solving the structure = ",round(timer.time() - startTime,2),"s")

            print("\n> Starting calculation of eigenfrequency")
            startTime = timer.time()
            listOfResponseFunctions["eigenfrequency"].CalculateValue()
            print("> Time needed for calculation of eigenfrequency = ",round(timer.time() - startTime,2),"s")

            communicator.reportValue("eigenfrequency", listOfResponseFunctions["eigenfrequency"].GetValue())

        # Calculation of gradient of eigenfrequency
        if communicator.isRequestingGradientOf("eigenfrequency"):

            print("\n> Starting calculation of gradients of eigenfrequency")
            startTime = timer.time()
            listOfResponseFunctions["eigenfrequency"].CalculateGradient()
            print("> Time needed for calculating gradients of eigenfrequency = ",round(timer.time() - startTime,2),"s")

            gradientForCompleteModelPart = listOfResponseFunctions["eigenfrequency"].GetGradient()
            communicator.reportGradient("eigenfrequency", gradientForCompleteModelPart)

        # Calculation of mass
        if communicator.isRequestingValueOf("mass"):

            print("\n> Starting calculation of mass")
            startTime = timer.time()
            listOfResponseFunctions["mass"].CalculateValue()
            constraintValue = listOfResponseFunctions["mass"].GetValue()
            print("> Time needed for calculation of mass = ",round(timer.time() - startTime,2),"s")

            communicator.reportValue("mass", constraintValue)

        # Calculation of gradient of mass
        if communicator.isRequestingGradientOf("mass"):

            print("\n> Starting calculation of gradient of mass")
            startTime = timer.time()
            listOfResponseFunctions["mass"].CalculateGradient()
            print("> Time needed for calculating gradient of mass = ",round(timer.time() - startTime,2),"s")

            gradientForCompleteModelPart = listOfResponseFunctions["mass"].GetGradient()
            communicator.reportGradient("mass", gradientForCompleteModelPart)

    # --------------------------------------------------------------------------
    def finalizeAfterOptimizationLoop( self ):
        csm_analysis.Finalize()

    # --------------------------------------------------------------------------

structureAnalyzer = kratosCSMAnalyzer()

# ======================================================================================================================================
# Optimization
# ======================================================================================================================================

optimizer.importAnalyzer( structureAnalyzer )
optimizer.optimize()

# ======================================================================================================================================