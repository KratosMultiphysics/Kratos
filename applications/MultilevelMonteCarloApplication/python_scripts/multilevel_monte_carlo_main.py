from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Importing the analysis stage classes of the different problems
from poisson_stochastic_analysis import SimulationScenario

# Import Monte Carlo library
import mc_utilities_new_analysis as mc
import cmlmc_utilities_new_analysis as mlmc


if __name__ == '__main__':

    '''
    |  V  |/ _|
    | \_/ | (_
    |_| |_|\__|
    '''

    '''set the ProjectParameters.json path'''
    project_parameters_path = "../tests/PoissonSquareTest/parameters_poisson_finer.json"
    '''customize setting parameters of the ML simulation'''
    settings_MC_simulation = KratosMultiphysics.Parameters("""
    {
        "tolerance" : 0.1,
        "cphi" : 5e-1,
        "batch_size" : 20,
        "convergence_criteria" : "MC_higher_moments_sequential_stopping_rule"
    }
    """)
    '''contruct MonteCarlo or MultilevelMonteCarlo class'''
    mc_manager = mc.MonteCarlo(settings_MC_simulation,project_parameters_path,SimulationScenario)
    '''execute algorithm'''
    mc_manager.Run()


    '''
    |  V  | | |  V  |/ _|
    | \_/ | |_| \_/ | (_
    |_| |_|___|_| |_|\__|
    '''

    '''set the ProjectParameters.json path'''
    project_parameters_path = "/home/kratos105b/Kratos/applications/MultilevelMonteCarloApplication/tests/PoissonSquareTest/parameters_poisson_coarse.json"
    '''customize setting parameters of the ML simulation'''
    settings_ML_simulation = KratosMultiphysics.Parameters("""
    {
        "tol0"                            : 0.25,
        "tolF"                            : 0.1,
        "cphi"                            : 1.0,
        "number_samples_screening"        : 25,
        "Lscreening"                      : 2,
        "Lmax"                            : 4,
        "initial_mesh_size"               : 0.5
    }
    """)
    '''customize setting parameters of the metric of the adaptive refinement utility'''
    custom_metric_refinement_parameters = KratosMultiphysics.Parameters("""
        {
            "hessian_strategy_parameters"           :{
                    "metric_variable"               : ["TEMPERATURE"],
                    "estimate_interpolation_error"  : false,
                    "interpolation_error"           : 0.004
            },
            "anisotropy_remeshing"                  : true,
            "anisotropy_parameters":{
                "reference_variable_name"           : "TEMPERATURE",
                "hmin_over_hmax_anisotropic_ratio"  : 0.15,
                "boundary_layer_max_distance"       : 1.0,
                "interpolation"                     : "Linear"
            },
            "local_gradient_variable"               : "TEMPERATURE"
        }
    """)
    '''customize setting parameters of the remesh of the adaptive refinement utility'''
    custom_remesh_refinement_settings = KratosMultiphysics.Parameters("""
        {
            "echo_level"                            : 0
        }
    """)
    '''contruct MultilevelMonteCarlo class'''
    mlmc_manager = mlmc.MultilevelMonteCarlo(settings_ML_simulation,project_parameters_path,custom_metric_refinement_parameters,custom_remesh_refinement_settings,SimulationScenario)
    mlmc_manager.RunScreening()
    mlmc_manager.Run()
