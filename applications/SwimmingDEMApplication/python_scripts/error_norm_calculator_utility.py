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

        self.p_characteristic = (1/2)*self.rho*self.u_characteristic**2

        self.model_part.AddNodalSolutionStepVariable(SDEM.VECTORIAL_ERROR)
        self.model_part.AddNodalSolutionStepVariable(SDEM.SCALAR_ERROR)

    def CalculateL2Norm(self):
        self.ComputeDofsErrors()
        self.velocity_L2_error_norm = self.VectorL2ErrorNorm()
        self.pressure_L2_error_norm = self.ScalarL2ErrorNorm()

        return self.velocity_L2_error_norm/self.u_characteristic, self.pressure_L2_error_norm, self.model_part

    def CalculateH1SemiNorm(self):

        self.velocity_H1_error_seminorm = self.VectorH1ErrorSemiNorm()
        self.pressure_H1_error_seminorm = self.ScalarH1ErrorSemiNorm()

        return self.velocity_H1_error_seminorm/self.u_characteristic, self.pressure_H1_error_seminorm/self.p_characteristic

    def ComputeDofsErrors(self):
        SDEM.ErrorNormCalculator().ComputeDofsErrors(self.model_part)

    def VectorL2ErrorNorm(self):
        return SDEM.ErrorNormCalculator().GetL2VectorErrorNorm(self.model_part, KratosMultiphysics.VELOCITY)

    def ScalarL2ErrorNorm(self):
        return SDEM.ErrorNormCalculator().GetL2ScalarErrorNorm(self.model_part, KratosMultiphysics.PRESSURE)

    def ScalarH1ErrorSemiNorm(self):
        return SDEM.ErrorNormCalculator().GetH1ScalarErrorSemiNorm(self.model_part)

    def VectorH1ErrorSemiNorm(self):
        return SDEM.ErrorNormCalculator().GetH1VectorErrorSemiNorm(self.model_part)
