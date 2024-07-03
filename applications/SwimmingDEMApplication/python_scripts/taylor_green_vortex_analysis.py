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

class TaylorGreenVortexAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, iteration, wave_number, maximum, difference, varying_parameters = Parameters("{}")):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        iteration -- The number of the fluid_model_part that is running
        model -- the container of the fluid model part.
        varying_parameters -- The whole project_parameters.
        """
        from KratosMultiphysics.SwimmingDEMApplication import hdf5_script
        self.kinetic_energy = hdf5_script.KineticEnergy(iteration)
        super().__init__(model, varying_parameters)
        self.project_parameters = varying_parameters
        self.GetModelAttributes()
        self.max_iteration = self.project_parameters['fluid_parameters']['solver_settings']['maximum_iterations'].GetInt()
        self.u_characteristic = self.project_parameters["error_projection_parameters"]["u_characteristic"].GetDouble()
        self.fix_pressure = self.project_parameters['fluid_parameters']["rve_settings"]["fix_pressure"].GetBool()
        self.max_deviation_from_mean = difference
        self.mean_alpha = maximum-self.max_deviation_from_mean
        self.k = wave_number
        # This model analysis is created to validate formulations so we have to make sure the fluid is computed in every time step

    def Initialize(self):
        super().Initialize()
        #self.kinetic_energy.WriteData(self.fluid_model_part, self.model_type, self.mean_alpha, self.time)

    def SetInitialParticlePosition(self):
        pass

    def InitializeVariablesWithNonZeroValues(self):
        pass

    def _GetFluidAnalysis(self):
        if not hasattr(self, '_fluid_phase_analysis'):
            import KratosMultiphysics.SwimmingDEMApplication.DEM_coupled_fluid_dynamics_analysis_periodicity as fluid_analysis
            self._fluid_phase_analysis = fluid_analysis.DEMCoupledFluidDynamicsAnalysisPeriodicity(self.model, self.project_parameters, self.vars_man)
            self._fluid_phase_analysis.main_path = self.main_path
        return self._fluid_phase_analysis

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = 1)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        #self.kinetic_energy.WriteData(self.fluid_model_part, self.model_type, self.mean_alpha, self.time)

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

    test = TaylorGreenVortexAnalysis(model, parameters)
    test.Run()