from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Importing the analysis stage classes of the different problems
from simulation_definition import SimulationScenario

# Import Monte Carlo library
import KratosMultiphysics.MultilevelMonteCarloApplication.mc_utilities as mc_utilities
import KratosMultiphysics.MultilevelMonteCarloApplication.mlmc_utilities as mlmc_utilities

# TODO: use json file instead of defining in the main file all the parameters

if __name__ == '__main__':

    """  __ __  __
        |  V  |/ _|
        | \_/ | (_
        |_| |_|\__|
    """

    # set the ProjectParameters.json path
    project_parameters_path = "problem_settings/parameters_poisson_square_2d_finer.json"
    # customize setting parameters of the MC simulation"""
    parameters_x_monte_carlo_path = "problem_settings/parameters_x_monte_carlo.json"
    # contruct MonteCarlo or MultilevelMonteCarlo class
    mc_manager = mc_utilities.MonteCarlo(parameters_x_monte_carlo_path,project_parameters_path,SimulationScenario)
    # execute algorithm
    mc_manager.Run()


    """ __ __ _   __ __  __
       |  V  | | |  V  |/ _|
       | \_/ | |_| \_/ | (_
       |_| |_|___|_| |_|\__|
    """

    # set the ProjectParameters.json path
    project_parameters_path = "problem_settings/parameters_poisson_square_2d_finer.json"
    # customize setting parameters of the MLMC simulation
    parameters_x_monte_carlo_path = "problem_settings/parameters_x_monte_carlo.json"
    # customize setting parameters of the metric of the adaptive refinement utility and setting parameters of the remesh of the adaptive refinement utility
    parameters_refinement_path = "problem_settings/parameters_refinement.json"
    # contruct MultilevelMonteCarlo class
    mlmc_manager = mlmc_utilities.MultilevelMonteCarlo(parameters_x_monte_carlo_path,project_parameters_path,parameters_refinement_path,SimulationScenario)
    mlmc_manager.Run()
