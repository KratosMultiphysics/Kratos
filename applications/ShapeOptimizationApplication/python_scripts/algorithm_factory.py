# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division


# ==============================================================================
def CreateOptimizationAlgorithm(optimization_settings, analyzer, communicator, model_part_controller):
    algorithm_name = optimization_settings["optimization_algorithm"]["name"].GetString()

    if algorithm_name == "steepest_descent":
        from .algorithm_steepest_descent import AlgorithmSteepestDescent
        return AlgorithmSteepestDescent(optimization_settings,
                                        analyzer,
                                        communicator,
                                        model_part_controller)
    elif algorithm_name == "gradient_projection":
        from .algorithm_gradient_projection import AlgorithmGradientProjection
        return AlgorithmGradientProjection(optimization_settings,
                                            analyzer,
                                            communicator,
                                            model_part_controller)
    elif algorithm_name == "penalized_projection":
        from .algorithm_penalized_projection import AlgorithmPenalizedProjection
        return AlgorithmPenalizedProjection(optimization_settings,
                                            analyzer,
                                            communicator,
                                            model_part_controller)
    elif algorithm_name == "trust_region":
        from .algorithm_trust_region import AlgorithmTrustRegion
        return AlgorithmTrustRegion(optimization_settings,
                                    analyzer,
                                    communicator,
                                    model_part_controller)
    elif algorithm_name == "bead_optimization":
        from .algorithm_bead_optimization import AlgorithmBeadOptimization
        return AlgorithmBeadOptimization(optimization_settings,
                                         analyzer,
                                         communicator,
                                         model_part_controller)
    elif algorithm_name == "relaxed_gradient_projection":
        from .algorithm_relaxed_gradient_projection import AlgorithmRelaxedGradientProjection
        return AlgorithmRelaxedGradientProjection(optimization_settings,
                                                  analyzer,
                                                  communicator,
                                                  model_part_controller)
    elif algorithm_name == "sequential_quadratic_programming":
        line_search = optimization_settings["optimization_algorithm"]["line_search"]["line_search_type"].GetString()
        if line_search == "quadratic_approximation":
            from .algorithm_sequential_quadratic_programming_with_line_search import AlgorithmSequentialQuadraticProgrammingWithLineSearch
            return AlgorithmSequentialQuadraticProgrammingWithLineSearch(optimization_settings,
                                                      analyzer,
                                                      communicator,
                                                      model_part_controller)
        else:
            from .algorithm_sequential_quadratic_programming import AlgorithmSequentialQuadraticProgramming
            return AlgorithmSequentialQuadraticProgramming(optimization_settings,
                                                      analyzer,
                                                      communicator,
                                                      model_part_controller)
    else:
        raise NameError("The following optimization algorithm is not supported by the algorithm factory: " + algorithm_name)

 # ==============================================================================
