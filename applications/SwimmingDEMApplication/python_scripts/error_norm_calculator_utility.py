from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import Vector
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import numpy as np
import sys

class ErrorNormCalculatorUtility:
    def __init__(self, model, parameters):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        """
        self.model_part = model

        self.u_characteristic = parameters["error_projection_parameters"]["u_characteristic"].GetDouble()

        self.rho = self.model_part.Elements.__iter__().__next__().Properties.GetValue(KratosMultiphysics.DENSITY)
        self.nu = self.model_part.Elements.__iter__().__next__().Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)/self.rho
        alpha_min = parameters["fluid_parameters"]["processes"]["initial_conditions_process_list"][0]["Parameters"]["benchmark_parameters"]["alpha_min"].GetDouble()
        if parameters["fluid_parameters"]["processes"]["initial_conditions_process_list"][0]["Parameters"]["benchmark_parameters"].Has("sigma"):
            sigma = parameters["fluid_parameters"]["processes"]["initial_conditions_process_list"][0]["Parameters"]["benchmark_parameters"]["sigma"].GetDouble()
        else:
            sigma = 0.0
        if parameters["fluid_parameters"]["processes"]["initial_conditions_process_list"][0]["Parameters"]["benchmark_parameters"].Has("omega"):
            omega = parameters["fluid_parameters"]["processes"]["initial_conditions_process_list"][0]["Parameters"]["benchmark_parameters"]["omega"].GetDouble()
        else:
            omega = 0.0


        Re = self.u_characteristic/self.nu
        Da = sigma/(self.nu * alpha_min)
        self.p_characteristic = self.u_characteristic * self.nu * (1 + omega/self.nu + Re + Da)

        self.model_part.AddNodalSolutionStepVariable(SDEM.VECTORIAL_ERROR)
        self.model_part.AddNodalSolutionStepVariable(SDEM.SCALAR_ERROR)

    def CalculateL2Norm(self, model_part):
        #self.ComputeDofsErrors()
        self.velocity_L2_error_norm = self.VectorL2ErrorNorm()
        self.pressure_L2_error_norm = self.ScalarL2ErrorNorm()

        return self.velocity_L2_error_norm/self.u_characteristic, self.pressure_L2_error_norm/self.p_characteristic, self.model_part

    def CalculateH1SemiNorm(self):

        self.velocity_H1_error_seminorm = self.VectorH1ErrorSemiNorm()
        self.pressure_H1_error_seminorm = self.ScalarH1ErrorSemiNorm()

        return self.velocity_H1_error_seminorm/self.u_characteristic, self.pressure_H1_error_seminorm/self.p_characteristic

    def ComputeDofsErrors(self):
        SDEM.ErrorNormCalculator().ComputeDofsErrors(self.model_part)

    def VectorL2ErrorNorm(self):
        return SDEM.ErrorNormCalculator().GetL2VectorErrorNorm(self.model_part, KratosMultiphysics.VELOCITY, KratosMultiphysics.EXACT_VELOCITY)

    def ScalarL2ErrorNorm(self):
        return SDEM.ErrorNormCalculator().GetL2ScalarErrorNorm(self.model_part, KratosMultiphysics.PRESSURE, SDEM.EXACT_PRESSURE)

    def ScalarH1ErrorSemiNorm(self):
        return SDEM.ErrorNormCalculator().GetH1ScalarErrorSemiNorm(self.model_part)

    def VectorH1ErrorSemiNorm(self):
        return SDEM.ErrorNormCalculator().GetH1VectorErrorSemiNorm(self.model_part)
