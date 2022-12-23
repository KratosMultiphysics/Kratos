# importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import Vector
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import numpy as np
import sys

class L2ErrorCalculatorUtility:
    def __init__(self, model, parameters):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        """
        self.model_part = model

        self.U = parameters["fluid_parameters"]["processes"]["initial_conditions_process_list"][0]["Parameters"]["benchmark_parameters"]["u_char"].GetDouble()

        for element in self.model_part.Elements:
            self.rho = element.Properties.GetValue(KratosMultiphysics.DENSITY)
            self.nu = element.Properties.GetValue(KratosMultiphysics.VISCOSITY)
            break

        porosity_field = [node.GetSolutionStepValue(KratosMultiphysics.FLUID_FRACTION) for node in self.model_part.Nodes]
        self.porosity_mean = np.mean(porosity_field)

        self.u_characteristic = self.U**2*0.01/self.porosity_mean

        self.p_characteristic = (1/2)*self.rho*self.u_characteristic**2

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
        self.ComputeDofsErrors()

        self.velocity_error_norm = self.VectorL2ErrorNorm()
        self.pressure_error_norm = self.ScalarL2ErrorNorm()

        self.reynolds_number = self.porosity_mean * self.u_characteristic / self.nu

        return self.velocity_error_norm/self.u_characteristic, self.pressure_error_norm/self.p_characteristic, self.error_model_part, self.reynolds_number, self.porosity_mean

    def ComputeDofsErrors(self):
        SDEM.L2ErrorNormCalculator().ComputeDofsErrors(self.error_model_part)

    def VectorL2ErrorNorm(self):
        return SDEM.L2ErrorNormCalculator().GetL2VectorErrorNorm(self.error_model_part)

    def ScalarL2ErrorNorm(self):
        return SDEM.L2ErrorNormCalculator().GetL2ScalarErrorNorm(self.error_model_part)
