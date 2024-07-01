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

class FluidFractionTestAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, iteration, damkohler = 0, omega = 0, L = 1, varying_parameters = Parameters("{}")):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        iteration -- The number of the fluid_model_part that is running
        model -- the container of the fluid model part.
        varying_parameters -- The whole project_parameters.
        """
        from KratosMultiphysics.SwimmingDEMApplication import hdf5_script
        self.projector_post_process = hdf5_script.ErrorProjectionPostProcessTool(iteration)
        super().__init__(model, varying_parameters)
        self.project_parameters = varying_parameters
        self.GetModelAttributes()
        self.max_iteration = self.project_parameters['fluid_parameters']['solver_settings']['maximum_iterations'].GetInt()
        self.lowest_alpha = self.project_parameters["fluid_parameters"]["processes"]["initial_conditions_process_list"][0]["Parameters"]["benchmark_parameters"]["alpha_min"].GetDouble()
        self.u_characteristic = self.project_parameters["error_projection_parameters"]["u_characteristic"].GetDouble()
        self.damkohler_number = damkohler
        self.omega = omega
        self.length = L
        # This model analysis is created to validate formulations so we have to make sure the fluid is computed in every time step

    def InitializeVariablesWithNonZeroValues(self):
        pass

    def Initialize(self):
        super().Initialize()
        self._GetSolver().ConstructErrorNormCalculator()

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = 1)

    def _CreateSolver(self):
        import fluid_fraction_test_solver as sdem_solver
        return sdem_solver.FluidFractionTestSolver(self.model,
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

        super(SwimmingDEMAnalysis, self).FinalizeSolutionStep()
        self.n_iteration_number = self.fluid_model_part.ProcessInfo[Kratos.NL_ITERATION_NUMBER]
        self.relax_alpha = self.fluid_model_part.ProcessInfo[SDEM.RELAXATION_ALPHA]

        porosity_field = [node.GetSolutionStepValue(Kratos.FLUID_FRACTION) for node in self.fluid_model_part.Nodes]
        self.porosity_mean = np.mean(porosity_field)
        for node in self.fluid_model_part.Nodes:
            self.nu = node.GetSolutionStepValue(Kratos.VISCOSITY)
            break
        self.reynolds_number = self.length*self.u_characteristic/self.nu

        self.velocity_L2_error_norm, self.pressure_L2_error_norm, self.error_model_part = self._GetSolver().CalculateL2ErrorNorm()

        self.velocity_H1_error_seminorm, self.pressure_H1_error_seminorm = self._GetSolver().CalculateH1ErrorSemiNorm()

        self.projector_post_process.WriteData(self.error_model_part,
                                            self.velocity_L2_error_norm,
                                            self.pressure_L2_error_norm,
                                            self.velocity_H1_error_seminorm,
                                            self.pressure_H1_error_seminorm,
                                            self.projection_type,
                                            self.model_type,
                                            self.subscale_type,
                                            self.porosity_mean,
                                            self.n_iteration_number,
                                            self.max_iteration,
                                            self.relax_alpha,
                                            self.lowest_alpha,
                                            self.damkohler_number,
                                            self.omega,
                                            self.reynolds_number,
                                            2)

    def TransferBodyForceFromDisperseToFluid(self):
        pass

    def GetVolumeDebugTool(self):
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