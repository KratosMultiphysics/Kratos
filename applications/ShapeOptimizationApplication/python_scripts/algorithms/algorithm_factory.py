# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================



# ==============================================================================
def CreateOptimizationAlgorithm(optimization_settings, analyzer, communicator, model_part_controller):
    algorithm_name = optimization_settings["optimization_algorithm"]["name"].GetString()

    design_variables = optimization_settings["design_variables"]
    for design_variable in design_variables:
        if design_variable["type"].GetString() == "free_thickness":
            if algorithm_name == "gradient_projection":
                from .algorithm_free_thickness_optimization_v3 import AlgorithmFreeThicknessOptimizationv3
                return AlgorithmFreeThicknessOptimizationv3(optimization_settings,
                                                            analyzer,
                                                            communicator,
                                                            model_part_controller)

            elif algorithm_name == "relaxed_gradient_projection":
                from .algorithm_free_thickness_optimization_v3_rgp import AlgorithmFreeThicknessOptimizationv3RGP
                return AlgorithmFreeThicknessOptimizationv3RGP(optimization_settings,
                                                                analyzer,
                                                                communicator,
                                                                model_part_controller)
            else:
                RuntimeError("""Only Gradient Projection or Relaxed Gradient Projection algorithm
                             is available for Free Thickness optimization.""")

        elif design_variable["type"].GetString() in ["thickness_parameter", "free_thickness_original_vm"]:
            if design_variables.size() > 1:
                RuntimeError("""Thickness parameter optimization is not available with shape optimization simultaneously!
                             Use Free Thickness optimization instead.""")
            # TODO: remove free_thickness_original_vm after numerical experiments
            if design_variable["type"].GetString() == "free_thickness_original_vm":
                if algorithm_name == "gradient_projection":
                    from .algorithm_free_thickness_optimization import AlgorithmFreeThicknessOptimization
                    return AlgorithmFreeThicknessOptimization(optimization_settings,
                                                            analyzer,
                                                            communicator,
                                                            model_part_controller)

                elif algorithm_name == "relaxed_gradient_projection":
                    from .algorithm_free_thickness_rgp import AlgorithmFreeThicknessRelaxedGradientProjection
                    return AlgorithmFreeThicknessRelaxedGradientProjection(optimization_settings,
                                                                        analyzer,
                                                                        communicator,
                                                                        model_part_controller)
                else:
                    RuntimeError("""Only Gradient Projection or Relaxed Gradient Projection algorithm
                                 is available for Free Thickness optimization.""")

            from .algorithm_thickness_optimization import AlgorithmThicknessOptimization
            return AlgorithmThicknessOptimization(optimization_settings,
                                                  analyzer,
                                                  communicator,
                                                  model_part_controller)
    else:
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
        elif algorithm_name == "shape_fraction_optimization":
            from .algorithm_shape_fraction_optimization import AlgorithmShapeFractionOptimization
            return AlgorithmShapeFractionOptimization(optimization_settings,
                                            analyzer,
                                            communicator,
                                            model_part_controller)

        else:
            raise NameError("The following optimization algorithm is not supported by the algorithm factory: " + algorithm_name)

 # ==============================================================================
