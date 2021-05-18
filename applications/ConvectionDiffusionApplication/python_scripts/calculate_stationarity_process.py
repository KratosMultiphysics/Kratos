from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector
import sys
import KratosMultiphysics.ConvectionDiffusionApplication as CD
import KratosMultiphysics.ConvectionDiffusionApplication.WriteStationarityError as WriteStationarityError
from KratosMultiphysics import *
import sys
import os

def Factory(settings, Model):
    if(type(settings) != Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    # self.compute_stationarity = settings["Parameters"]["compute_stationarity"].GetBool()
    if settings["Parameters"]["compute_stationarity"].GetBool():
        return CalculateStationarityProcessUtility(Model,settings["Parameters"])

## All the processes python should be derived from "Process"
class CalculateStationarityProcessUtility(Kratos.Process):
    def __init__(self, Model, settings):
        super().__init__()
        default_settings = Kratos.Parameters("""
            {
                "file_name" : "PLEASE_SPECIFY_HDF5_FILENAME",
                "model_part_name" : "please_specify_model_part_name",
                "compute_stationarity" : true,
                "domain_size" : 2,
                "error_projection_parameters" : {
                	"u_characteristic"  : 1.0
            	}
            }
            """)

        settings.ValidateAndAssignDefaults(default_settings)
        self.parameters = settings
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.u_characteristic = settings["error_projection_parameters"]["u_characteristic"].GetDouble()
        self.problem_name = self.parameters["file_name"].GetString()
        self.velocity_error = []
        self.concentration_error = []
        self.step = []
        for element in self.fluid_model_part.Elements:
            rho = element.Properties.GetValue(Kratos.DENSITY)
            break
        self.p_characteristic = (1/2)*rho*self.u_characteristic**2

    def ExecuteFinalizeSolutionStep(self):
        step = self.fluid_model_part.ProcessInfo[Kratos.STEP]
        velocity_error, concentration_error, self.fluid_model_part=self.CalculateL2()
        self.velocity_error.append(velocity_error)
        self.concentration_error.append(concentration_error)
        self.step.append(step)
        WriteStationarityError.WriteStationarityErrorToHdf5(self.fluid_model_part, velocity_error, concentration_error, self.problem_name)

    def CalculateL2(self):
        self.ComputeDofsErrors(self.fluid_model_part)

        self.velocity_error_norm = self.VectorL2ErrorNorm(self.fluid_model_part)
        self.concentration_error_norm = self.ScalarL2ErrorNorm(self.fluid_model_part)
        return self.velocity_error_norm/self.u_characteristic, self.concentration_error_norm/self.p_characteristic, self.fluid_model_part

    def ComputeDofsErrors(self, fluid_model_part):
        CD.CalculateStationarityProcess().ComputeDofsErrorsCD(self.fluid_model_part)

    def VectorL2ErrorNorm(self, fluid_model_part):
        return CD.CalculateStationarityProcess().GetL2VectorErrorNormCD(self.fluid_model_part)

    def ScalarL2ErrorNorm(self, fluid_model_part):
        return CD.CalculateStationarityProcess().GetL2ScalarErrorNormCD(self.fluid_model_part)




