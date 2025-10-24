import KratosMultiphysics

def GetDefaultInputSettings():
    default_settings = KratosMultiphysics.Parameters("""
        {
            "thermal_integration_scheme"     : "forward_euler",
            "numerical_integration_method"   : "adaptive_simpson",
            "thermal_solve_frequency"        : 1,
            "voronoi_tesselation_frequency"  : 1000,
	        "porosity_update_frequency"      : 1000,
            "automatic_solve_frequency"      : false,
            "compute_forces"                 : true,
            "compute_motion"                 : true,
            "compute_direct_conduction"      : true,
            "compute_indirect_conduction"    : false,
            "compute_convection"             : false,
            "compute_radiation"              : false,
            "compute_heat_generation"        : false,
            "compute_adjusted_contact"       : false,
            "direct_conduction_model"        : "batchelor_obrien_simple",
            "indirect_conduction_model"      : "surrounding_layer",
            "nusselt_correlation"            : "sphere_hanz_marshall",
            "radiation_model"                : "continuum_zhou",
            "heat_generation_model"          : ["sliding_friction"],
            "adjusted_contact_model"         : "zhou",
            "voronoi_method"                 : "tesselation",
	        "porosity_method"                : "average_alpha_shape",
            "min_conduction_distance"        : 0.0000000275,
            "max_conduction_distance"        : 1.0,
            "conduction_radius"              : 1.0,
            "fluid_layer_thickness"          : 0.4,
            "isothermal_core_radius"         : 0.5,
            "max_radiation_distance"         : 2.0,
            "heat_generation_ratio"          : 1.0,
            "global_porosity"                : 0.0,
            "alpha_shape_parameter"          : 1.2,
            "integral_tolerance"             : 0.000001,
            "heat_map_corners"               : [[0,0,0],[1,1,1]],
            "heat_map_subdivisions"          : [10,10,10],
            "global_fluid_properties"        : {
                "fluid_density"              : 1.0,
                "fluid_viscosity"            : 1.0,
                "fluid_thermal_conductivity" : 1.0,
                "fluid_heat_capacity"        : 1.0,
                "fluid_temperature"          : 0.0,
                "fluid_velocity_X"           : 0.0,
                "fluid_velocity_Y"           : 0.0,
                "fluid_velocity_Z"           : 0.0
            }
        }""")
    return default_settings
