import KratosMultiphysics as Kratos
import KratosMultiphysics.SwimmingDEMApplication as SDEM
from KratosMultiphysics import Model, Parameters, Logger
import numpy as np

import os
import sys

file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
sys.path.insert(0, dir_path)

from swimming_DEM_analysis import SwimmingDEMAnalysis

import swimming_DEM_procedures as SDP

class FlowPastCylinderAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, iteration, Re, varying_parameters = Parameters("{}")):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        iteration -- The number of the fluid_model_part that is running
        model -- the container of the fluid model part.
        varying_parameters -- The whole project_parameters.
        """
        from KratosMultiphysics.SwimmingDEMApplication import hdf5_script
        super().__init__(model, varying_parameters)
        self.project_parameters = varying_parameters
        self.GetModelAttributes()
        self.max_iteration = self.project_parameters['fluid_parameters']['solver_settings']['maximum_iterations'].GetInt()
        self.reynolds_number = Re
        self.lowest_alpha = self.project_parameters["fluid_parameters"]["processes"]["initial_conditions_process_list"][0]["Parameters"]["benchmark_parameters"]["alpha_min"].GetDouble()
        self.u_characteristic = self.project_parameters["error_projection_parameters"]["u_characteristic"].GetDouble()
        self.pressure = hdf5_script.Pressure(iteration)
        # This model analysis is created to validate formulations so we have to make sure the fluid is computed in every time step

    def InitializeVariablesWithNonZeroValues(self):
        pass

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = 1)

    def _CreateSolver(self):
        import flow_past_cylinder_solver as sdem_solver
        return sdem_solver.FlowPastCylinderSolver(self.model,
                                                self.project_parameters,
                                                self.GetFieldUtility(),
                                                self._GetFluidAnalysis()._GetSolver(),
                                                self._GetDEMAnalysis()._GetSolver(),
                                                self.vars_man)

    def SetEmbeddedTools(self):
        pass

    def ComputePostProcessResults(self):
        pass

    def GetDerivativeRecoveryStrategy(self):
        pass

    def FinalizeSolutionStep(self):
        # printing if required
        if self._GetSolver().CannotIgnoreFluidNow():
            self._GetFluidAnalysis().FinalizeSolutionStep()

        self._GetDEMAnalysis().FinalizeSolutionStep()

        # coupling checks (debugging)
        if self.debug_info_counter.Tick():
            self.dem_volume_tool.UpdateDataAndPrint(
                self.project_parameters["fluid_domain_volume"].GetDouble())

        #super(SwimmingDEMAnalysis, self).FinalizeSolutionStep()

        self.pressure.WriteData(self.fluid_model_part, self.model_type, self.lowest_alpha, self.reynolds_number)

    def TransferBodyForceFromDisperseToFluid(self):
        pass

    def GetVolumeDebugTool(self):
        pass

    def SetInitialParticlePosition(self):
        pass

    def GetModelAttributes(self):
        if self.project_parameters["fluid_parameters"]["solver_settings"]["formulation"]["use_orthogonal_subscales"].GetBool() == True:
            self.projection_type = 'OSS'
        else:
            self.projection_type = 'ASGS'

        element_type = self.project_parameters["fluid_parameters"]["solver_settings"]["formulation"]["element_type"].GetString()

        if element_type == "advmsDEM":
            self.model_type = 'Drew model'
            self.subscale_type = 'dynamic'
        elif element_type == "aqsvmsDEM":
            self.model_type = 'Drew model'
            self.subscale_type = 'quasi-static'
        elif element_type == "qsvmsDEM":
            self.model_type = 'Jackson model'
            self.subscale_type = 'quasi-static'
        elif element_type == "dvmsDEM":
            self.model_type = 'Jackson model'
            self.subscale_type = 'dynamic'

if __name__ == "__main__":
    # Setting parameters

    with open('ProjectParameters.json','r') as parameter_file:
        parameters = Parameters(parameter_file.read())

    # Create Model
    model = Model()

    # To avoid too many prints
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)

    test = FluidFractionTestAnalysis(model, parameters)
    test.Run()