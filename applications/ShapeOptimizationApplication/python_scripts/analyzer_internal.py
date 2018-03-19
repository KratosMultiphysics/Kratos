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
from structural_mechanics_analysis import StructuralMechanicsAnalysis
import response_function_factory
import time as timer

# ==============================================================================
class KratosInternalAnalyzer( (__import__("analyzer_base")).AnalyzerBaseClass ):
    # --------------------------------------------------------------------------
    def __init__( self, project_parameters, model_part ):
        self.csm_analysis = StructuralMechanicsAnalysis(project_parameters, model_part)

        self.listOfResponseFunctions = response_function_factory.CreateListOfResponseFunctions(project_parameters["optimization_settings"], model_part)

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

        # Calculation of gradients of strain energy
        if communicator.isRequestingGradientOf("strain_energy"):

            print("\n> Starting calculation of gradient of strain energy")
            startTime = timer.time()
            self.listOfResponseFunctions["strain_energy"].CalculateGradient()
            print("> Time needed for calculating gradient of strain energy = ",round(timer.time() - startTime,2),"s")

            gradientForCompleteModelPart = self.listOfResponseFunctions["strain_energy"].GetGradient()
            communicator.reportGradient("strain_energy", gradientForCompleteModelPart)

        # Calculation of value of mass
        if communicator.isRequestingValueOf("mass"):

            print("\n> Starting calculation of value of mass")
            startTime = timer.time()
            self.listOfResponseFunctions["mass"].CalculateValue()
            constraintValue = self.listOfResponseFunctions["mass"].GetValue()
            print("> Time needed for calculation of value of mass = ",round(timer.time() - startTime,2),"s")

            communicator.reportValue("mass", constraintValue)

        # Calculation of gradients of mass
        if communicator.isRequestingGradientOf("mass"):

            print("\n> Starting calculation of gradient of mass")
            startTime = timer.time()
            self.listOfResponseFunctions["mass"].CalculateGradient()
            print("> Time needed for calculating gradient of mass = ",round(timer.time() - startTime,2),"s")

            gradientForCompleteModelPart = self.listOfResponseFunctions["mass"].GetGradient()
            communicator.reportGradient("mass", gradientForCompleteModelPart)

        # Calculation of eigenfrequency
        if communicator.isRequestingValueOf("eigenfrequency"):

            print("\n> Starting StructuralMechanicsApplication to solve structure")
            startTime = timer.time()
            self.csm_analysis.InitializeTimeStep()
            self.csm_analysis.SolveTimeStep()
            self.csm_analysis.FinalizeTimeStep()
            print("> Time needed for solving the structure = ",round(timer.time() - startTime,2),"s")

            print("\n> Starting calculation of eigenfrequency")
            startTime = timer.time()
            self.listOfResponseFunctions["eigenfrequency"].CalculateValue()
            print("> Time needed for calculation of eigenfrequency = ",round(timer.time() - startTime,2),"s")

            communicator.reportValue("eigenfrequency", self.listOfResponseFunctions["eigenfrequency"].GetValue())

        # Calculation of gradient of eigenfrequency
        if communicator.isRequestingGradientOf("eigenfrequency"):

            print("\n> Starting calculation of gradients of eigenfrequency")
            startTime = timer.time()
            self.listOfResponseFunctions["eigenfrequency"].CalculateGradient()
            print("> Time needed for calculating gradients of eigenfrequency = ",round(timer.time() - startTime,2),"s")

            gradientForCompleteModelPart = self.listOfResponseFunctions["eigenfrequency"].GetGradient()
            communicator.reportGradient("eigenfrequency", gradientForCompleteModelPart)

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop( self ):
        self.csm_analysis.Finalize()

# ==============================================================================
