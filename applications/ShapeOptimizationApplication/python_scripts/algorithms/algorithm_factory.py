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
            if design_variable["filtering_aproach"].GetString() == "initial":

                    if algorithm_name == "gradient_projection":
                        from .algorithm_shape_thickness_initial_filtering_gp import AlgorithmShapeThicknessInitialFilteringGradientProjection
                        return AlgorithmShapeThicknessInitialFilteringGradientProjection(optimization_settings,
                                                                                         analyzer,
                                                                                         communicator,
                                                                                         model_part_controller)

                    elif algorithm_name == "relaxed_gradient_projection":
                        from .algorithm_shape_thickness_initial_filtering_rgp import AlgorithmShapeThicknessInitialFilteringRelaxedGradientProjection
                        return AlgorithmShapeThicknessInitialFilteringRelaxedGradientProjection(optimization_settings,
                                                                                                analyzer,
                                                                                                communicator,
                                                                                                model_part_controller)
                    else:
                        RuntimeError("""Only Gradient Projection or Relaxed Gradient Projection algorithm
                                    is available for Free Thickness optimization.""")

            elif design_variable["filtering_aproach"].GetString() == "updated":
                    if design_variables.size() > 1:
                        RuntimeError("""Thickness optimization with updated filtering approach is not available with shape optimization simultaneously!
                                    Use initial filtering instead.""")

                    if algorithm_name == "gradient_projection":
                        from .algorithm_thickness_updated_filtering_gp import AlgorithmThicknessUpdatedFilteringGradientProjection
                        return AlgorithmThicknessUpdatedFilteringGradientProjection(optimization_settings,
                                                                                    analyzer,
                                                                                    communicator,
                                                                                    model_part_controller)

                    elif algorithm_name == "relaxed_gradient_projection":
                        from .algorithm_thickness_updated_filtering_rgp import AlgorithmThicknessUpdatedFilteringRelaxedGradientProjection
                        return AlgorithmThicknessUpdatedFilteringRelaxedGradientProjection(optimization_settings,
                                                                                           analyzer,
                                                                                           communicator,
                                                                                           model_part_controller)
                    else:
                        RuntimeError("""Only Gradient Projection or Relaxed Gradient Projection algorithm
                                    is available for Free Thickness optimization.""")
            else:
                NameError("The following filtering approach is not supported by the algorithm factory: " + design_variable["filtering_aproach"].GetString())

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
