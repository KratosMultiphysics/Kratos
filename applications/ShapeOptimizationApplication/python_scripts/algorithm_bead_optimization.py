# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from algorithm_base import OptimizationAlgorithm
import mapper_factory
import data_logger_factory
from custom_timer import Timer
from custom_variable_utilities import WriteDictionaryDataOnNodalVariable
import math, copy

# ==============================================================================
class AlgorithmBeadOptimization(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = Parameters("""
        {
            "name"                        : "bead_optimization",
            "bead_height"                 : 1.0,
            "bead_direction_mode"         : 2,
            "penalty_factor"              : 1000.0,
            "gradient_ratio"              : 0.1,
            "max_outer_iterations"        : 300,
            "max_inner_iterations"        : 50,
            "min_inner_iterations"        : 1,
            "inner_iteration_tolerance"   : 1e-3,
            "outer_iteration_tolerance"   : 1e-3,
            "line_search" : {
                "line_search_type"           : "manual_stepping",
                "normalize_search_direction" : true,
                "step_size"                  : 1.0
            }
        }""")
        self.algorithm_settings =  optimization_settings["optimization_algorithm"]
        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.objectives = optimization_settings["objectives"]

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surface = model_part_controller.GetDesignSurface()

        self.Mapper = mapper_factory.CreateMapper(self.design_surface, optimization_settings["design_variables"]["filter"])
        self.DataLogger = data_logger_factory.CreateDataLogger(model_part_controller, communicator, optimization_settings)

        self.OptimizationUtilities = OptimizationUtilities(self.design_surface, optimization_settings)

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("The augmented lagrange algorithm for bead optimization only supports one objective function!")
        if self.constraints.size() > 0:
            raise RuntimeError("The augmented lagrange algorithm for bead does not allow for any constraints!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.only_obj = self.objectives[0]

        self.bead_height = self.algorithm_settings["bead_height"].GetDouble()
        self.direction_mode = self.algorithm_settings["bead_direction_mode"].GetInt()
        self.penalty_factor = self.algorithm_settings["penalty_factor"].GetDouble()
        self.gradient_ratio = self.algorithm_settings["gradient_ratio"].GetDouble()
        self.max_outer_iterations = self.algorithm_settings["max_outer_iterations"].GetInt()
        self.max_inner_iterations = self.algorithm_settings["max_inner_iterations"].GetInt()
        self.min_inner_iterations = self.algorithm_settings["min_inner_iterations"].GetInt()
        self.inner_iteration_tolerance = self.algorithm_settings["inner_iteration_tolerance"].GetDouble()
        self.outer_iteration_tolerance = self.algorithm_settings["outer_iteration_tolerance"].GetDouble()
        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()

        self.alpha0 = {}
        self.lower_bounds = {}
        self.upper_bounds = {}
        if self.direction_mode == 1:
            for node in self.design_surface.Nodes:
                self.alpha0[node.Id] = 0.5
                self.lower_bounds[node.Id] = 0.0
                self.upper_bounds[node.Id] = 1.0
        elif self.direction_mode == -1:
            for node in self.design_surface.Nodes:
                self.alpha0[node.Id] = -0.5
                self.lower_bounds[node.Id] = -1.0
                self.upper_bounds[node.Id] = 0.0
        elif self.direction_mode == 2:
            for node in self.design_surface.Nodes:
                self.alpha0[node.Id] = 0.0001
                self.lower_bounds[node.Id] = -1.0
                self.upper_bounds[node.Id] = 1.0
        else:
            raise RuntimeError("Specified bead direction mode not supported!")

        self.lambda0 = 0.0
        self.constraint_scaling_0 = 1

        self.model_part_controller.InitializeMeshController()
        self.Mapper.InitializeMapping()
        self.analyzer.InitializeBeforeOptimizationLoop()
        self.DataLogger.InitializeDataLogging()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        alpha = self.alpha0
        current_lambda = self.lambda0
        constraint_scaling = self.constraint_scaling_0
        overall_iteration = 0
        previos_L = None

        # Compute normals once in the beginning
        self.model_part_controller.ComputeUnitSurfaceNormals()

        for outer_iteration in range(1,self.max_outer_iterations+1):
            for inner_iteration in range(1,self.max_inner_iterations+1):

                overall_iteration = overall_iteration+1
                timer.StartNewLap()

                print("\n>=====================================================================================")
                print("> ",timer.GetTimeStamp(),": Starting iteration ",outer_iteration,".",inner_iteration,".",overall_iteration,"(outer.inner.overall)")
                print(">=====================================================================================\n")

                # Initialize new shape
                self.model_part_controller.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
                self.model_part_controller.SetReferenceMeshToMesh()

                # Analyze shape
                self.communicator.initializeCommunication()
                self.communicator.requestValueOf(self.only_obj["identifier"].GetString())
                self.communicator.requestGradientOf(self.only_obj["identifier"].GetString())

                self.analyzer.AnalyzeDesignAndReportToCommunicator(self.design_surface, overall_iteration, self.communicator)

                objective_value = self.communicator.getStandardizedValue(self.only_obj["identifier"].GetString())
                objGradientDict = self.communicator.getStandardizedGradient(self.only_obj["identifier"].GetString())
                WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, DF1DX)

                # Compute sensitivities w.r.t. scalar design variable alpha
                dF1dalpha = {}
                for node in self.design_surface.Nodes:
                    raw_gradient = node.GetSolutionStepValue(DF1DX)
                    normal = node.GetSolutionStepValue(NORMALIZED_SURFACE_NORMAL)
                    dF1dalpha[node.Id] = self.bead_height*(raw_gradient[0]*normal[0] + raw_gradient[1]*normal[1] + raw_gradient[2]*normal[2])

                # Map gradient using temporarily an auxiliary variable
                aux_var = {}
                for node in self.design_surface.Nodes:
                    aux_var[node.Id] = [dF1dalpha[node.Id],0.0,0.0]
                WriteDictionaryDataOnNodalVariable(aux_var, self.design_surface, DF1DX_MAPPED)
                self.Mapper.MapToDesignSpace(DF1DX_MAPPED, DF1DX_MAPPED)

                dF1dalpha_mapped = {}
                for node in self.design_surface.Nodes:
                    dF1dalpha_mapped[node.Id] = node.GetSolutionStepValue(DF1DX_MAPPED_X)

                # Compute penalization term
                penalty_gradient = {}
                penalty_value = 0.0
                if self.direction_mode == 1:
                    for node in self.design_surface.Nodes:
                        alpha_i = alpha[node.Id]
                        penalty_value = penalty_value + alpha_i-alpha_i**2
                        penalty_gradient[node.Id] = 1-2*alpha_i

                elif self.direction_mode == -1:
                    for node in self.design_surface.Nodes:
                        alpha_i = alpha[node.Id]
                        penalty_value = penalty_value + -alpha_i-alpha_i**2
                        penalty_gradient[node.Id] = -1-2*alpha_i

                elif self.direction_mode == 2:
                    for node in self.design_surface.Nodes:
                        alpha_i = alpha[node.Id]
                        penalty_value = penalty_value + -alpha_i**2+1
                        penalty_gradient[node.Id] = -2*alpha_i

                # Compute Lagrange value
                L = objective_value + current_lambda*constraint_scaling*penalty_value + self.penalty_factor*(constraint_scaling*penalty_value)**2

                # Compute gradient of Lagrange function
                dLdalpha = {}
                for node in self.design_surface.Nodes:
                    dLdalpha[node.Id] = dF1dalpha_mapped[node.Id] + current_lambda*constraint_scaling*penalty_gradient[node.Id]

                # Normalization using infinity norm
                dLdalpha_for_normalization = copy.deepcopy(dLdalpha)
                dLdalpha_for_normalization.update({key: value**2 for key, value in dLdalpha_for_normalization.items()})

                for node in self.design_surface.Nodes:
                    if alpha[node.Id]==self.lower_bounds[node.Id] or alpha[node.Id]==self.upper_bounds[node.Id]:
                        dLdalpha_for_normalization[node.Id] = 0.0

                max_value = math.sqrt(max(dLdalpha_for_normalization.values()))
                if max_value != 0:
                    dLdalpha.update({key: value/max_value for key, value in dLdalpha.items()})

                # Compute updated design variable
                alpha_new = {}
                for node in self.design_surface.Nodes:
                    dalpha = -self.step_size*dLdalpha[node.Id]
                    alpha_new[node.Id] = alpha[node.Id] + dalpha

                # Enforce bounds
                for node in self.design_surface.Nodes:
                    alpha_new[node.Id] = max(alpha_new[node.Id], self.lower_bounds[node.Id])
                    alpha_new[node.Id] = min(alpha_new[node.Id], self.upper_bounds[node.Id])

                # Map design variables using temporarily an auxiliary variable
                aux_var = {}
                for node in self.design_surface.Nodes:
                    aux_var[node.Id] = [alpha_new[node.Id],0.0,0.0]
                WriteDictionaryDataOnNodalVariable(aux_var, self.design_surface, VECTOR_VARIABLE)
                self.Mapper.MapToGeometrySpace(VECTOR_VARIABLE, VECTOR_VARIABLE_MAPPED)

                # Compue actual shape update
                for node in self.design_surface.Nodes:
                    normal = node.GetSolutionStepValue(NORMALIZED_SURFACE_NORMAL)
                    original_alpha = node.GetSolutionStepValue(VECTOR_VARIABLE_X)
                    mapped_alpha = node.GetSolutionStepValue(VECTOR_VARIABLE_MAPPED_X)

                    vectorized_alpha = original_alpha * normal
                    shape_change = mapped_alpha * normal * self.bead_height

                    previos_shape_change = node.GetSolutionStepValue(SHAPE_CHANGE)
                    shape_update = shape_change-previos_shape_change

                    node.SetSolutionStepValue(CONTROL_POINT_CHANGE, vectorized_alpha)
                    node.SetSolutionStepValue(SHAPE_CHANGE, shape_change)
                    node.SetSolutionStepValue(SHAPE_UPDATE, shape_update)

                # Log current optimization step
                additional_values_to_log = {}
                additional_values_to_log["step_size"] = self.algorithm_settings["line_search"]["step_size"].GetDouble()
                additional_values_to_log["outer_iteration"] = outer_iteration
                additional_values_to_log["inner_iteration"] = inner_iteration
                additional_values_to_log["lagrange_value"] = lagrange_value
                additional_values_to_log["penalty_value"] = penalty_value
                additional_values_to_log["penalty_lambda"] = penalty_lambda

                self.DataLogger.LogCurrentValues(overall_iteration, additional_values_to_log)
                self.DataLogger.LogCurrentDesign(overall_iteration)

                # Iterate
                alpha = alpha_new
                previos_L = L

                # Check convergence
                if inner_iteration>=self.min_inner_iterations:
                    L_change_relative = 1-L/previos_L
                    if L_change_relative<self.inner_iteration_tolerance:
                        break

                print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
                print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

            # Update lambda
            current_lambda = current_lambda + self.penalty_factor*constraint_scaling*penalty_value

            if outer_iteration==1:
                norm_term_1 = 0.0
                norm_term_2 = 0.0

                for node in self.design_surface.Nodes:
                    temp_value = dF1dalpha_mapped[node.Id]
                    norm_term_1 = norm_term_1+temp_value**2
                    temp_value = current_lambda*penalty_gradient[node.Id]

                    norm_term_2 = norm_term_2+temp_value**2
                norm_term_1 = math.sqrt(norm_term_1)
                norm_term_2 = math.sqrt(norm_term_2)

                constraint_scaling = self.gradient_ratio*norm_term_1/norm_term_2

            print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
            print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

        # Check convergence of outer loop
        if overall_iteration > 1 :

            # Check if maximum iterations were reached
            if overall_iteration == self.max_outer_iterations:
                print("\n> Maximal iterations of optimization problem reached!")
                break

            # # Check for relative tolerance
            # relativeChangeOfObjectiveValue = self.DataLogger.GetValue("rel_change_obj", overall_iteration)
            # if abs(relativeChangeOfObjectiveValue) < self.relativeTolerance:
            #     print("\n> Optimization problem converged within a relative objective tolerance of ",self.relativeTolerance,"%.")
            #     return True

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.DataLogger.FinalizeDataLogging()
        self.analyzer.FinalizeAfterOptimizationLoop()

# ==============================================================================
