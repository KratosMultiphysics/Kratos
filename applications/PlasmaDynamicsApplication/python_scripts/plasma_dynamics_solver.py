import sys
import os
import math

from KratosMultiphysics import Logger, Parameters

from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_solver import SwimmingDEMSolver
# lib_path = os.path.abspath(os.path.join(__file__, '..', 'SwimmingDEMApplication', 'python_scripts','derivative_recovery'))
# sys.path.append(lib_path)
import KratosMultiphysics.SwimmingDEMApplication.derivative_recovery.derivative_recovery_strategy as derivative_recoverer

import KratosMultiphysics.PlasmaDynamicsApplication.plasma_dynamics_procedures as PDP
import KratosMultiphysics.PlasmaDynamicsApplication.CFD_DEM_for_plasma_dynamics_coupling as fluid_DEM_coupling


def Say(*args):
    Logger.PrintInfo("PlasmaDynamics", *args)
    Logger.Flush()


class PlasmaDynamicsSolver(SwimmingDEMSolver):

    @classmethod
    def GetDefaultParameters(cls):
        # default settings string in json format
        default_settings = Parameters("""{
        "echo_level" : 1,

        "gravity_parameters" : {
            "modulus" : 9.81,
            "direction" : [0.0, 0.0, -1.0]
        },

        "time_stepping" : {
            "automatic_time_step" : true,
            "time_step" : 0.001
        },

        "problem_data"  : {
            "problem_name"  : "",
            "echo_level" : 1,
            "start_time" : 0.0,
            "end_time"   : 1,
            "parallel_type": "OpenMP",
            "number_of_threads": 1
        },

        "ElementType" : "SwimmingDEMElement",
        "body_force_per_unit_mass_variable_name" : "BODY_FORCE",
        "error_projection_parameters"   :{
            "u_characteristic"  : 1.0
        },
        "do_print_results_option" : true,
        "output_interval" : 0.5,

        "processes" : {},

        "coupling" : {
            "coupling_level_type" : 1,
            "coupling_weighing_type" : 2,
            "coupling_weighing_type_comment" : "{fluid_to_DEM, DEM_to_fluid, fluid_fraction} = {lin, lin, imposed} (-1), {lin, const, const} (0), {lin, lin, const} (1), {lin, lin, lin} (2), averaging method (3)",
            "interaction_start_time" : 0.0,

            "forward_coupling" : {},

            "backward_coupling" : {}
        },

        "frame_of_reference" : {},

        "non_newtonian_fluid" : {},

        "similarity" : {},

        "stationarity" : {},

        "debug_tool_cycle" : 10,
        "debug_tool_cycle_comment" : " number of 'ticks' per debug computations cycle",
        "print_debug_info_option" : false,
        "print_debug_info_option_comment" : " print a summary of global physical measures",
        "do_process_analytic_data" : true,
        "fluid_domain_volume" : 1.0,
        "fluid_domain_volume_comment" : "write down the volume you know it has, if available",

        "full_particle_history_watcher" : "Empty",


        "gradient_calculation_type" : 1,
        "gradient_calculation_type_comment" : "(Not calculated (0), volume-weighed average(1), Superconvergent recovery(2))",
        "material_acceleration_calculation_type" : 1,
        "laplacian_calculation_type" : 0,
        "laplacian_calculation_type_comment" : "(Not calculated (0), Finite element projection (1), Superconvergent recovery(2))",
        "vorticity_calculation_type" : 5,
        "store_full_gradient_option" : false,
        "add_each_hydro_force_option" : true,
        "add_each_hydro_force_option_comment" : " add each of the hydrodynamic forces (drag, lift and virtual mass)",
        "pressure_grad_recovery_type" : 0,
        "recovery_echo_level" : 1,
        "store_fluid_pressure_option" : false,

        "print_distance_option" : false,
        "print_steps_per_plot_step" : 1,
        "print_particles_results_option" : false,
        "make_results_directories_option" : true,
        "make_results_directories_option_comment": "results are written into a folder (../results) inside the problem folder",
        "print_particles_results_cycle" : 1,
        "print_particles_results_cycle_comment" : " number of 'ticks' per printing cycle",

        "drag_modifier_type" : 2,
        "drag_modifier_type_comment" : " Hayder (2), Chien (3)",

        "json_output_process" : [],
        "sdem_output_processes" : {},
        "properties": [{}],

        "fluid_phase_parameters" : {},

        "custom_fluid" : {},

        "dem_parameters" : {},

        "custom_dem" : {},

        "dem_nodal_results" : {},

        "fluid_nodal_results" : {}

        }""")

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings


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
        #nodal_area_process_parameters["model_part_name"].SetString(self.fluid_solver.fluid_solver.main_model_part.Name)
        
        # print('=========================')
        # print(self.fluid_solver.fluid_solver.main_model_part.Name)
        # print('=========================')
        
        nodal_area_process_parameters["domain_size"].SetInt(self.fluid_domain_dimension)
        the_mesh_moves = False
        
        # print('=========================')
        # print(self.fluid_solver.settings)
        # print('=========================')
        
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


    def __init__(self, model, project_parameters, field_utility, fluid_phase_solver, dem_solver, variables_manager):
#         super().__init__(model, project_parameters, 
# 			 field_utility, 
# 			 fluid_solver.fluid_solver, 
# 			 dem_solver, 
# 			 variables_manager)


        # Validate settings
        self.field_utility = field_utility
        self.vars_man = variables_manager
        self.fluid_domain_dimension = project_parameters["fluid_phase_parameters"]["solver_settings"]["domain_size"].GetInt()
        self.fluid_phase_solver = fluid_phase_solver
        
        self.fluid_solver = fluid_phase_solver.fluid_solver
        self.thermal_solver = fluid_phase_solver.thermal_solver
        self.dem_solver = dem_solver
        
        self.project_parameters = self._ValidateSettings(project_parameters)
        self.next_time_to_solve_fluid = project_parameters['problem_data']['start_time'].GetDouble()
        self.coupling_level_type = project_parameters["coupling"]["coupling_level_type"].GetInt()
        self.interaction_start_time = project_parameters["coupling"]["interaction_start_time"].GetDouble()
        self.integration_scheme = project_parameters["custom_dem"]["translational_integration_scheme"].GetString()
        
        self.fluid_dt = self.fluid_solver.settings["time_stepping"]["time_step"].GetDouble()
        
        self.do_solve_dem = project_parameters["custom_dem"]["do_solve_dem"].GetBool()
        self.solve_system = not self.project_parameters["custom_fluid"]["fluid_already_calculated"].GetBool()

        self.fluid_step = 0
        self.calculating_fluid_in_current_step = True
        self.first_DEM_iteration = True
        self.SetHistoryForceOptions()
        self.ConstructStationarityTool()
        self.ConstructDerivativeRecoverer()
        self.ConstructHistoryForceUtility()
    
        
        # Call the base Python solver constructor
        super(SwimmingDEMSolver, self).__init__(model, project_parameters)


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
    

