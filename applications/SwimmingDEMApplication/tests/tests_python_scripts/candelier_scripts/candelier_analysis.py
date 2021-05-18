import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters, Vector
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import math
import os
import sys
file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
sys.path.insert(0, dir_path)
import candelier
import candelier_parameters as candelier_pp

def Cross(a, b):
    c0 = a[1]*b[2] - a[2]*b[1]
    c1 = a[2]*b[0] - a[0]*b[2]
    c2 = a[0]*b[1] - a[1]*b[0]
    return Vector([c0, c1, c2])

from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import Say
import KratosMultiphysics.SwimmingDEMApplication.parameters_tools as PT

class CandelierBenchmarkAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, varying_parameters=Parameters("{}")):
        super().__init__(model, varying_parameters)
        self._GetSolver().is_rotating_frame = self.project_parameters["frame_of_reference"]["frame_type"].GetInt()
        self._GetDEMAnalysis().mdpas_folder_path = os.path.join(self._GetDEMAnalysis().main_path, 'candelier_tests')
        candelier_pp.include_history_force = self.vars_man.do_include_history_force
        candelier_pp.include_lift = PT.RecursiveFindParametersWithCondition(self.project_parameters["properties"],
                                                                            'vorticity_induced_lift_parameters',
                                                                            condition=lambda value: not (value['name'].GetString()=='default'))

        Kratos.Logger.PrintInfo("SwimmingDEM", 'candelier_pp.include_lift ', candelier_pp.include_lift)
        candelier.sim = candelier.AnalyticSimulator(candelier_pp)
        self.project_parameters["custom_fluid"]["fluid_already_calculated"].SetBool(True)
        self.project_parameters.AddEmptyValue("load_derivatives").SetBool(False)

    def GetEmbeddedCounter(self):
        return SDP.Counter(is_dead=True)

    def GetBackwardCouplingCounter(self):
        return SDP.Counter(1, 4, 0)

    def GetDebugInfo(self):
        return SDP.Counter(self.project_parameters["debug_tool_cycle"].GetInt(), 1, is_dead = 1)

    def SetUpResultsDatabase(self):
        import candelier_scripts.candelier_hdf5 as candelier_hdf5
        self.results_database = candelier_hdf5.ResultsCandelier(self.project_parameters, self.main_path)

    def PerformZeroStepInitializations(self):
        # Impose initial velocity to be that of the fluid for the x/y-components
        # and the terminal velocity for the z-component
        candelier.sim.CalculateNonDimensionalVars()
        terminal_velocity_z = 2. / 9 * 9.81 * candelier_pp.a ** 2 / (candelier_pp.nu * candelier_pp.rho_f) * (candelier_pp.rho_f - candelier_pp.rho_p)

        for node in self.spheres_model_part.Nodes:
            r = Kratos.Vector([node.X, node.Y, node.Z])
            v0 = Kratos.Vector([candelier_pp.u0, candelier_pp.v0, terminal_velocity_z])

            if self._GetSolver().is_rotating_frame:
                v0 = self._GetSolver().GetVelocityRelativeToMovingFrame(r_rel = r, v_glob = v0)

            node.SetSolutionStepValue(Kratos.VELOCITY, v0)
            node.Fix(Kratos.VELOCITY_Z)
            node.SetSolutionStepValue(Kratos.VELOCITY_OLD, v0)
            node.Fix(Kratos.VELOCITY_OLD_Z)
            node.SetSolutionStepValue(Kratos.FLUID_VEL_PROJECTED, v0)

            if candelier_pp.include_lift:
                vorticity = Kratos.Vector([0.0, 0.0, 2.0 * candelier_pp.omega])
                node.SetSolutionStepValue(Kratos.FLUID_VORTICITY_PROJECTED, vorticity)

    def _CreateSolver(self):
        import tests_python_scripts.candelier_scripts.candelier_dem_solver as sdem_solver
        return sdem_solver.CandelierDEMSolver(self.model,
                                              self.project_parameters,
                                              self.GetFieldUtility(),
                                              self._GetFluidAnalysis()._GetSolver(),
                                              self._GetDEMAnalysis()._GetSolver(),
                                              self.vars_man)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        analytic_coors = [0., 0., 0.]
        candelier.sim.CalculatePosition(analytic_coors, self.time * candelier_pp.omega)
        node = [node for node in self.spheres_model_part.Nodes][0]
        Say('analytic', analytic_coors)
        Say('calculated', [value / candelier_pp.R for value in (node.X, node.Y, node.Z)])

        if self.project_parameters["do_print_results_option"].GetBool():
            for node in self.spheres_model_part.Nodes:
                r = Kratos.Vector([node.X, node.Y, node.Z])
                self.results_database.MakeReading(self.time, r)

                coor_calculated = [node.X, node.Y, node.Z]

                if self._GetSolver().is_rotating_frame:
                    sin = math.sin(self._GetSolver().omega * self.time)
                    cos = math.cos(self._GetSolver().omega * self.time)
                    coor_calculated[0] = node.X * cos - node.Y * sin
                    coor_calculated[1] = node.X * sin + node.Y * cos

                self.radial_error = self.results_database.CalculateError(self.time, coor_calculated)
                self.error_time = self.time

