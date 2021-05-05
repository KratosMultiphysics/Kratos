import sys
import os
import math

from KratosMultiphysics import Logger, Parameters

from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_solver import SwimmingDEMSolver
lib_path = os.path.abspath(os.path.join(__file__, '..', 'SwimmingDEMApplication', 'python_scripts','derivative_recovery'))
sys.path.append(lib_path)
import KratosMultiphysics.SwimmingDEMApplication.derivative_recovery.derivative_recovery_strategy as derivative_recoverer

import KratosMultiphysics.PlasmaDynamicsApplication.plasma_dynamics_procedures as PDP
import KratosMultiphysics.PlasmaDynamicsApplication.CFD_DEM_for_plasma_dynamics_coupling as fluid_DEM_coupling


def Say(*args):
    Logger.PrintInfo("PlasmaDynamics", *args)
    Logger.Flush()


class PlasmaDynamicsSolver(SwimmingDEMSolver):

    # TODO
    # def GetDefaultParameters(cls):
	#     pass


    def _ValidateSettings(self, project_parameters):
        default_processes_settings = Parameters("""{
                "python_module" : "calculate_nodal_area_process",
                "kratos_module" : "KratosMultiphysics.SwimmingDEMApplication",
                "process_name"  : "CalculateNodalAreaProcess",
                "Parameters"    : {
                    "model_part_name" : "FluidModelPart",
                    "domain_size" : 3,
                    "fixed_mesh": false
                }
            }

        """)

        if not project_parameters["processes"].Has('non_optional_solver_processes'):
            project_parameters["processes"].AddEmptyArray("non_optional_solver_processes")

        else: # reconstruct non_optional_solver_processes list making sure calculate_nodal_area_process is not added twice
            non_optional_processes_list = list(project_parameters["processes"]["non_optional_solver_processes"])
            project_parameters["processes"].Remove("non_optional_solver_processes")
            project_parameters["processes"].AddEmptyArray("non_optional_solver_processes")

            for process in non_optional_processes_list:
                if process["python_module"].GetString() != 'calculate_nodal_area_process':
                    project_parameters["processes"]["non_optional_solver_processes"].Append(process)

        non_optional_solver_processes = project_parameters["processes"]["non_optional_solver_processes"]
        non_optional_solver_processes.Append(default_processes_settings)
        nodal_area_process_parameters = non_optional_solver_processes[non_optional_solver_processes.size() -1]["Parameters"]
        nodal_area_process_parameters["model_part_name"].SetString(self.fluid_solver.main_model_part.Name)
        print('=========================')
        print(nodal_area_process_parameters["model_part_name"])
        print('=========================')
        nodal_area_process_parameters["domain_size"].SetInt(self.fluid_domain_dimension)
        the_mesh_moves = False
        if self.fluid_solver.settings.Has('move_mesh_flag'):
            the_mesh_moves = self.fluid_solver.settings["move_mesh_flag"].GetBool()
            nodal_area_process_parameters["fixed_mesh"].SetBool(not the_mesh_moves)
        elif self.fluid_solver.settings.Has('time_integration_settings'):
            the_mesh_moves = self.fluid_solver.settings["time_integration_settings"]["move_mesh_flag"].GetBool()
            nodal_area_process_parameters["fixed_mesh"].SetBool(not the_mesh_moves)
        elif self.fluid_solver.settings["solvers"][0]["Parameters"]["time_integration_settings"].Has('move_mesh_flag'):
            the_mesh_moves = self.fluid_solver.settings["solvers"][0]["Parameters"]["time_integration_settings"]["move_mesh_flag"].GetBool()
            nodal_area_process_parameters["fixed_mesh"].SetBool(not the_mesh_moves)
        self.move_mesh_flag = the_mesh_moves
        return project_parameters


    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        super().__init__(model, project_parameters, 
			 field_utility, 
			 fluid_solver, 
			 dem_solver, 
			 variables_manager)


    def _ConstructProjectionModule(self):
        self.h_min = 0.01
        n_balls = 1
        fluid_volume = 10
        self.project_parameters.AddEmptyValue("n_particles_in_depth").SetInt(int(math.sqrt(n_balls / fluid_volume)))

        # creating a projection module for the fluid-DEM coupling
        projection_module = fluid_DEM_coupling.ProjectionModule(
        	self.fluid_solver.main_model_part,
        	self.dem_solver.spheres_model_part,
        	self.dem_solver.all_model_parts.Get("RigidFacePart"),
        	self.project_parameters,
        	self.vars_man.coupling_dem_vars,
        	self.vars_man.coupling_fluid_vars,
        	self.vars_man.time_filtered_vars,
        	flow_field=self.field_utility,
        	domain_size=self.fluid_domain_dimension
        	)

        projection_module.UpdateDatabase(self.h_min)

        return projection_module


    def GetRecoveryCounter(self):
        there_is_something_to_recover = (
            self.project_parameters["coupling"]["coupling_level_type"].GetInt() or
            self.project_parameters["print_PRESSURE_GRADIENT_option"].GetBool())
        return PDP.Counter(1, 1, there_is_something_to_recover)


    def ConstructDerivativeRecoverer(self):
        self.derivative_recovery_counter = self.GetRecoveryCounter()

        self.recovery = derivative_recoverer.DerivativeRecoveryStrategy(
            self.project_parameters,
            self.fluid_solver.main_model_part,
            PDP.FunctionsCalculator(self.fluid_domain_dimension))


    def GetStationarityCounter(self):
        return PDP.Counter(
            steps_in_cycle=self.project_parameters["stationarity"]["time_steps_per_stationarity_step"].GetInt(),
            beginning_step=self.project_parameters["stationarity"]["time_steps_before_first_assessment"].GetInt(),
            is_active=self.project_parameters["stationarity"]["stationary_problem_option"].GetBool())    


    def GetHistoryForceQuadratureCounter(self):
        for prop in self.project_parameters["properties"].values():
            if prop["hydrodynamic_law_parameters"].Has("history_force_parameters"):
                history_force_parameters =  prop["hydrodynamic_law_parameters"]["history_force_parameters"]
                if history_force_parameters.Has("time_steps_per_quadrature_step"):
                    time_steps_per_quadrature_step = history_force_parameters["time_steps_per_quadrature_step"].GetInt()

                    return PDP.Counter(steps_in_cycle=time_steps_per_quadrature_step, beginning_step=1)

        return PDP.Counter(is_dead=True)
    

