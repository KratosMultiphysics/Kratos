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

class GranularTemperatureAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, varying_parameters = Parameters("{}"), bc_conditions = Parameters("{}")):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the fluid model part.
        varying_parameters -- The whole project_parameters.
        """
        from KratosMultiphysics.SwimmingDEMApplication import hdf5_script
        self.averaging_variables = hdf5_script.AveragingVariablesPostProcessTool()

        self.bc = bc_conditions
        self.fix_top_velocity = self.bc["Boundary_at_top"]["fix_velocity"]
        self.top_velocity = self.bc["Boundary_at_top"]["top_velocity"].GetVector()
        self.top_pressure = self.bc["Boundary_at_top"]["top_pressure"].GetDouble()
        self.fix_bottom_velocity = self.bc["Boundary_at_bottom"]["fix_velocity"]
        self.bottom_velocity = self.bc["Boundary_at_bottom"]["bottom_velocity"].GetVector()
        self.bottom_pressure = self.bc["Boundary_at_bottom"]["bottom_pressure"].GetDouble()
        self.impose_top_pressure = self.bc["Boundary_at_top"]["fix_pressure"].GetBool()
        self.redefine_particle_boundary_interval = self.bc["bc_application_interval"].GetInt()
        self.bc_start_time = self.bc["bc_start_time"].GetDouble()
        self.bc_applied = False
        super().__init__(model, varying_parameters)
        # This model analysis is created to validate formulations so we have to make sure the fluid is computed in every time step

    def ModifyInputParametersForCoherence(self):
        # Making all time steps exactly commensurable
        output_time = self.project_parameters["output_interval"].GetDouble()
        self.output_time = int(output_time / self.time_step) * self.time_step
        self.project_parameters["output_interval"].SetDouble(self.output_time)

        self.project_parameters["dem_parameters"]["MaxTimeStep"].SetDouble(self.time_step)
        translational_scheme_name = self.project_parameters["custom_dem"]["translational_integration_scheme"].GetString()
        self.project_parameters["dem_parameters"]["TranslationalIntegrationScheme"].SetString(translational_scheme_name)

        self.SetDoSolveDEMVariable()

    def Initialize(self):
        super().Initialize()

    def InitializeSolutionStep(self):
        if (self.spheres_model_part.ProcessInfo[Kratos.TIME_STEPS] % self.redefine_particle_boundary_interval == 0.0 and self.spheres_model_part.ProcessInfo[Kratos.TIME] >= self.bc_start_time) or self._GetSolver().first_DEM_iteration:
            self.SetBoundaryConditions()
        super().InitializeSolutionStep()


    def SetBoundaryConditions(self):
        for elem in self.spheres_model_part.Elements:
            node = elem.GetGeometry()[0]
            #if node.Y > 0.0008*55:
            if node.Y > 0.0008*55:
                if self.fix_top_velocity[0].GetBool():
                    node.SetSolutionStepValue(Kratos.VELOCITY_X, self.top_velocity[0])
                    node.Fix(Kratos.VELOCITY_X)
                if self.fix_top_velocity[1].GetBool():
                    node.SetSolutionStepValue(Kratos.VELOCITY_Y, self.top_velocity[1])
                # if node.Y > 0.0008*59:
                #     node.SetSolutionStepValue(Kratos.VELOCITY_Y, 0.0)
                #     node.Fix(Kratos.VELOCITY_Y)
                if self.fix_top_velocity[2].GetBool():
                    node.SetSolutionStepValue(Kratos.VELOCITY_Z, self.top_velocity[2])
                    node.Fix(Kratos.VELOCITY_Z)
                if self.impose_top_pressure:
                    applied_force_y = -np.pi * node.GetSolutionStepValue(Kratos.RADIUS)**2 * self.top_pressure
                    node.SetSolutionStepValue(Kratos.EXTERNAL_APPLIED_FORCE_Y, applied_force_y)
                    node.Fix(Kratos.EXTERNAL_APPLIED_FORCE_Y)
            elif node.Y < 0.0008*5:
                if self.fix_bottom_velocity[0].GetBool():
                    node.SetSolutionStepValue(Kratos.VELOCITY_X, self.bottom_velocity[0])
                    node.Fix(Kratos.VELOCITY_X)
                if self.fix_bottom_velocity[1].GetBool():
                    node.SetSolutionStepValue(Kratos.VELOCITY_Y, self.bottom_velocity[1])
                    node.Fix(Kratos.VELOCITY_Y)
                if self.fix_bottom_velocity[2].GetBool():
                    node.SetSolutionStepValue(Kratos.VELOCITY_Z, self.bottom_velocity[2])
                    node.Fix(Kratos.VELOCITY_Z)
            else:
                node.Free(Kratos.VELOCITY_X)
                node.Free(Kratos.VELOCITY_Y)
                node.Free(Kratos.VELOCITY_Z)
                node.SetSolutionStepValue(Kratos.EXTERNAL_APPLIED_FORCE_Y, 0.0)
                node.Free(Kratos.EXTERNAL_APPLIED_FORCE_Y)

    def Finalize(self):
        import swimming_DEM_analysis as SDEM_script
        SDEM_script.Say('Finalizing simulation...\n')
        if self.do_print_results:
            self.swimming_DEM_gid_io.finalize_results()

        self.PerformFinalOperations(self.time)

        self.TellFinalSummary(self.time, self.step, self._GetSolver().fluid_step)

    def FluidInitialize(self):
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

    def OutputSolutionStep(self):
        # printing if required
        if self.print_counter.Tick():
            self._Print()
            self._disperse_phase_analysis.OutputSolutionStep()
            self.averaging_variables.WriteData(self.spheres_model_part)

        #Calling to AnalysisStage's method
        super(SwimmingDEMAnalysis, self).OutputSolutionStep()

    def _CreateSolver(self):
        import KratosMultiphysics.SwimmingDEMApplication.granular_temperature_solver as granular_temperature_solver
        return granular_temperature_solver.GranularTemperatureSolver(self.model,
                                                     self.project_parameters,
                                                     self.GetFieldUtility(),
                                                     self._GetFluidAnalysis()._GetSolver(),
                                                     self._GetDEMAnalysis()._GetSolver(),
                                                     self.vars_man)

    def GetBackwardCouplingCounter(self):
        return False

    def GetRecoveryCounter(self):
        return False

    def GetStationarityCounter(self):
        return False