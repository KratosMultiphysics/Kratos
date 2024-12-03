import KratosMultiphysics

def GetDefaultInputParameters():

    default_settings = KratosMultiphysics.Parameters(
    """{
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

        "error_projection_parameters"   :{
            "u_characteristic"  : 1.0
        },

        "do_print_results_option" : true,
        "output_interval" : 0.5,
        "use_fluid_static" : 0,

        "processes" : {
        },

        "coupling" : {
            "coupling_level_type" : 1,
            "coupling_weighing_type" : 2,
            "coupling_weighing_type_comment" : "{fluid_to_DEM, DEM_to_fluid, fluid_fraction} = {lin, lin, imposed} (-1), {lin, const, const} (0), {lin, lin, const} (1), {lin, lin, lin} (2), averaging method (3)",
            "interaction_start_time" : 0.0,

            "forward_coupling" : {
                "time_averaging_type" : 0
            },

            "gentle_coupling_initiation": {
                "initiation_interval": 0
            },

            "backward_coupling" : {
                "backward_time_interval" : 1,
                "meso_scale_length" : 0.2,
                "meso_scale_length_comment" : " the radius of the support of the averaging function for homogenization (<=0 for automatic calculation)",
                "shape_factor" : 0.5,
                "filter_velocity_option" : false,
                "apply_time_filter_to_fluid_fraction_option" : false,
                "min_fluid_fraction" : 0.2,
                "fluid_fraction_grad_type" : 0,
                "calculate_diffusivity_option" : false,
                "viscosity_modification_type" : 0,
                "averaging_time_interval" : 1
            }
        },

        "frame_of_reference" : {
            "frame_type": 0,
            "angular_velocity_of_frame_X" : 0.0,
            "angular_velocity_of_frame_Y" : 0.0,
            "angular_velocity_of_frame_Z" : 0.0,
            "angular_velocity_of_frame_old_X" : 0.0,
            "angular_velocity_of_frame_old_Y" : 0.0,
            "angular_velocity_of_frame_old_Z" : 0.0,
            "acceleration_of_frame_origin_X" : 0.0,
            "acceleration_of_frame_origin_Y" : 0.0,
            "acceleration_of_frame_origin_Z" : 0.0,
            "angular_acceleration_of_frame_X" : 0.0,
            "angular_acceleration_of_frame_Y" : 0.0,
            "angular_acceleration_of_frame_Z" : 0.0,
            "frame_rotation_axis_initial_point" : [0.0,0.0,0.0],
            "frame_rotation_axis_final_point" : [0.0,0.0,1.0],
            "angular_velocity_magnitude" : 0.0
        },

        "non_newtonian_fluid" : {
            "non_newtonian_option" : false,
            "yield_stress" : 0.0,
            "regularization_coefficient" : 0.0,
            "power_law_tol" : 0.0,
            "power_law_k" : 0.0,
            "power_law_n" : 0.0
        },

        "similarity" : {
            "similarity_transformation_type" : 0,
            "similarity_transformation_type_comment" : " no transformation (0), Tsuji (1)",
            "model_over_real_diameter_factor" : 1.0,
            "model_over_real_diameter_factor_comment": " not active if similarity_transformation_type = 0"
        },

        "stationarity" : {
            "stationary_problem_option" : false,
            "stationary_problem_option_comment" : " stationary, stop calculating the fluid after it reaches the stationary state",
            "tolerance" : 1e-3,
            "tolerance_comment": " fraction of the historically-averaged, spatial-averaged time derivative of the pressure that is considered significant",
            "time_steps_per_stationarity_step" : 15,
            "time_steps_per_stationarity_step_comment" : " number of fluid time steps between consecutive assessment of stationarity steps",
            "time_steps_per_analytic_processing_step" : 1,
            "time_steps_before_first_assessment" : 4
        },

        "debug_tool_cycle" : 10,
        "debug_tool_cycle_comment" : " number of 'ticks' per debug computations cycle",
        "print_debug_info_option" : false,
        "print_debug_info_option_comment" : " print a summary of global physical measures",
        "do_process_analytic_data" : true,
        "fluid_domain_volume" : 1.0,
        "fluid_domain_volume_comment" : "write down the volume you know it has, if available",

        "full_particle_history_watcher" : "Empty",


            "gradient_calculation_type" : 1,
            "compute_exact_L2" : false,
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

        "fluid_parameters" : {},

        "custom_fluid" : {
            "fluid_already_calculated" : false,
            "embedded_option" : false,
            "embedded_option_comment" : "the embedded domain tools are to be used",
            "do_impose_flow_from_field_option" : false,
            "flow_in_porous_medium_option" : false,
            "flow_in_porous_medium_option_comment" : " the porosity is an imposed field",
            "flow_in_porous_DEM_medium_option" : false,
            "flow_in_porous_DEM_medium_option_comment" : "the DEM part is kept static",
            "body_force_on_fluid_option" : true,
            "ALE_option" : false,
            "fluid_model_type" : 1,
            "fluid_model_type_comment" : " untouched, velocity incremented by 1/fluid_fraction (0), modified mass conservation only (1)"
        },

        "dem_parameters" : {
            "seed" : 42
        },

        "custom_dem" : {
            "do_solve_dem" : true,
            "do_search_neighbours" : true,
            "do_search_dem_neighbours" : true,
            "do_search_fem_neighbours" : true,
            "type_of_dem_inlet" : "VelocityImposed",
            "translational_integration_scheme" : "Hybrid_Bashforth"
        },

        "dem_nodal_results" : {
            "REYNOLDS_NUMBER" : false,
            "SLIP_VELOCITY" : false,
            "RADIUS" : false,
            "ANGULAR_VELOCITY" : false,
            "ELASTIC_FORCES" : false,
            "CONTACT_FORCES" : false,
            "TOTAL_FORCES" : false,
            "EXTERNAL_APPLIED_FORCE" : false,
            "CATION_CONCENTRATION" : false,
            "PRESSURE" : false,
            "PRESSURE_GRAD_PROJECTED" : false,
            "HYDRODYNAMIC_FORCE" : false,
            "HYDRODYNAMIC_MOMENT" : false,
            "FLUID_VEL_PROJECTED" : false,
            "FLUID_VEL_PROJECTED_RATE" : false,
            "FLUID_VEL_LAPL_PROJECTED" : false,
            "FLUID_VEL_LAPL_RATE_PROJECTED" : false,
            "FLUID_ACCEL_PROJECTED" : false,
            "FLUID_ACCEL_FOLLOWING_PARTICLE_PROJECTED" : false,
            "FLUID_FRACTION_GRADIENT_PROJECTED" : false,
            "FLUID_VISCOSITY_PROJECTED" : false,
            "FLUID_FRACTION_PROJECTED" : false,
            "BUOYANCY" : false,
            "DRAG_FORCE" : false,
            "VIRTUAL_MASS_FORCE" : false,
            "BASSET_FORCE" : false,
            "LIFT_FORCE" : false,
            "IMPACT_WEAR" : false,
            "NON_DIMENSIONAL_VOLUME_WEAR" : false
        },

        "fluid_nodal_results" : {
            "MATERIAL_ACCELERATION" : false,
            "VELOCITY_GRADIENT" : false,
            "PRESSURE_GRADIENT" : false,
            "AVERAGED_FLUID_VELOCITY" : false,
            "FLUID_FRACTION" : false,
            "FLUID_FRACTION_OLD" : false,
            "DISPERSE_FRACTION" : false,
            "PARTICLE_VEL_FILTERED" : false,
            "TIME_AVERAGED_ARRAY_3" : false,
            "PHASE_FRACTION" : false,
            "FLUID_FRACTION_GRADIENT" : false,
            "FLUID_FRACTION_RATE" : false,
            "HYDRODYNAMIC_REACTION" : false,
            "MEAN_HYDRODYNAMIC_REACTION" : false,
            "POWER_LAW_N" : false,
            "POWER_LAW_K" : false,
            "YIELD_STRESS" : false,
            "GEL_STRENGTH" : false,
            "VISCOSITY" : false,
            "DISTANCE" : false,
            "SLIP_VELOCITY" : false,
            "VORTICITY" : false,
            "VELOCITY_LAPLACIAN" : false,
            "BODY_FORCE" : false,
            "CONDUCTIVITY" : false,
            "VECTORIAL_ERROR" : false,
            "VECTORIAL_ERROR_1" : false,
            "VELOCITY_LAPLACIAN_RATE" : false,
            "MESH_VELOCITY" : false
        }

        }""")

    return default_settings
