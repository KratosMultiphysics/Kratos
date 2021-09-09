from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import Vector
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import sys

class L2ErrorCalculatorUtility:
    def __init__(self, model, parameters):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        """
        self.model_part = model

        self.u_characteristic = parameters["error_projection_parameters"]["u_characteristic"].GetDouble()

        for element in self.model_part.Elements:
            rho = element.Properties.GetValue(KratosMultiphysics.DENSITY)
            break

        self.p_characteristic = (1/2)*rho*self.u_characteristic**2

        self.model = KratosMultiphysics.Model()

        self.element_name = "Element3D4N"

        self.error_model_part = self.model.CreateModelPart("ErrorModelPart")

        self.error_model_part.AddNodalSolutionStepVariable(SDEM.VECTORIAL_ERROR)
        self.error_model_part.AddNodalSolutionStepVariable(SDEM.SCALAR_ERROR)
        self.error_model_part.AddNodalSolutionStepVariable(SDEM.ERROR_X)
        self.error_model_part.AddNodalSolutionStepVariable(SDEM.ERROR_Y)
        self.error_model_part.AddNodalSolutionStepVariable(SDEM.ERROR_Z)
        self.error_model_part.AddNodalSolutionStepVariable(SDEM.ERROR_P)

        model_part_cloner = KratosMultiphysics.ConnectivityPreserveModeler()
        model_part_cloner.GenerateModelPart(self.model_part, self.error_model_part, self.element_name)

        self.error_model_part.ProcessInfo = self.model_part.ProcessInfo

    def CalculateL2(self):
        self.ComputeDofsErrors(self.error_model_part)

        self.velocity_error_norm = self.VectorL2ErrorNorm(self.error_model_part)
        self.pressure_error_norm = self.ScalarL2ErrorNorm(self.error_model_part)

        return self.velocity_error_norm/self.u_characteristic, self.pressure_error_norm/self.p_characteristic, self.error_model_part

    def ComputeDofsErrors(self, error_model_part):
        SDEM.L2ErrorNormCalculator().ComputeDofsErrors(self.error_model_part)

    def VectorL2ErrorNorm(self, error_model_part):
        return SDEM.L2ErrorNormCalculator().GetL2VectorErrorNorm(self.error_model_part)

    def ScalarL2ErrorNorm(self, error_model_part):
        return SDEM.L2ErrorNormCalculator().GetL2ScalarErrorNorm(self.error_model_part)
