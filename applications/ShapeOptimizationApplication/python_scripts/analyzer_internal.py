from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# For time measures
import time as timer

# ==============================================================================
class KratosInternalAnalyzer( (__import__("analyzer_base")).AnalyzerBaseClass ):

    # --------------------------------------------------------------------------
    def __init__( self, project_parameters, model_part ):

        import response_function_factory
        self.listOfResponseFunctions = response_function_factory.CreateListOfResponseFunctions(project_parameters["optimization_settings"], model_part)

        import structural_mechanics_analysis
        self.csm_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(project_parameters, model_part)

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop( self ):
        self.csm_analysis.Initialize()

    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):

         # Calculation of value of strain energy
        if communicator.isRequestingValueOf("strain_energy"):

            print("\n> Starting StructuralMechanicsApplication to solve structure")
            startTime = timer.time()
            self.csm_analysis.InitializeTimeStep()
            self.csm_analysis.SolveTimeStep()
            self.csm_analysis.FinalizeTimeStep()
            print("> Time needed for solving the structure = ",round(timer.time() - startTime,2),"s")

            print("\n> Starting calculation of strain energy")
            startTime = timer.time()
            self.listOfResponseFunctions["strain_energy"].CalculateValue()
            print("> Time needed for calculation of strain energy = ",round(timer.time() - startTime,2),"s")

            communicator.reportValue("strain_energy", self.listOfResponseFunctions["strain_energy"].GetValue())

        # Calculation of value of mass
        if communicator.isRequestingValueOf("mass"):

            print("\n> Starting calculation of value of mass")
            startTime = timer.time()
            self.listOfResponseFunctions["mass"].CalculateValue()
            constraintValue = self.listOfResponseFunctions["mass"].GetValue()
            print("> Time needed for calculation of value of mass = ",round(timer.time() - startTime,2),"s")

            communicator.reportValue("mass", constraintValue)

        # Calculation of gradients of strain energy
        if communicator.isRequestingGradientOf("strain_energy"):

            print("\n> Starting calculation of gradient of strain energy")
            startTime = timer.time()
            self.listOfResponseFunctions["strain_energy"].CalculateGradient()
            print("> Time needed for calculating gradient of strain energy = ",round(timer.time() - startTime,2),"s")

            gradientForCompleteModelPart = self.listOfResponseFunctions["strain_energy"].GetGradient()
            communicator.reportGradient("strain_energy", gradientForCompleteModelPart)

        # Calculation of gradients of mass
        if communicator.isRequestingGradientOf("mass"):

            print("\n> Starting calculation of gradient of mass")
            startTime = timer.time()
            self.listOfResponseFunctions["mass"].CalculateGradient()
            print("> Time needed for calculating gradient of mass = ",round(timer.time() - startTime,2),"s")

            gradientForCompleteModelPart = self.listOfResponseFunctions["mass"].GetGradient()
            communicator.reportGradient("mass", gradientForCompleteModelPart)

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop( self ):
        self.csm_analysis.Finalize()

# ==============================================================================
