import KratosMultiphysics as Core
import KratosMultiphysics.GeoMechanicsApplication as Geo

def Factory(settings, model):
    if not isinstance(settings, Core.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")

    parameters = settings["Parameters"]

    default_parameters = Core.Parameters("""
    {
        "model_part_name": "",
        "increment_strategy": "fixed",
        "factor_increment": 0.1,
        "factor_list": [],
        "max_trials": 50,
        "monitored_node_ids": [],
        "failure_criterion": {
            "type": "non_convergence",
            "stop_stiffness_ratio": 0.02,
            "extrapolation_points": 3
        },
        "retry": {
            "max_cuts": 0,
            "cut_factor": 0.5,
            "bracket_tolerance": 0.0,
            "bracket_max_increment": 0.05,
            "bracket_max_stiffness_ratio": 0.3
        },
        "adaptive_increment": {
            "displacement_measure": "max_norm",
            "slope_averaging_steps": 2,
            "min_increment": 0.001,
            "max_increment": 0.2,
            "max_growth_factor": 1.5,
            "target_factor": 0.0,
            "target_approach_fraction": 0.25
        },
        "srf_indicators": {
            "enabled": false,
            "output_file_name": "srf_indicators.txt"
        }
    }
    """)

    parameters.ValidateAndAssignDefaults(default_parameters)
    parameters["srf_indicators"].ValidateAndAssignDefaults(default_parameters["srf_indicators"])
    parameters["failure_criterion"].ValidateAndAssignDefaults(default_parameters["failure_criterion"])
    parameters["retry"].ValidateAndAssignDefaults(default_parameters["retry"])
    parameters["adaptive_increment"].ValidateAndAssignDefaults(default_parameters["adaptive_increment"])

    return Geo.ApplyIncrementalCPhiReductionProcess(model, parameters)
