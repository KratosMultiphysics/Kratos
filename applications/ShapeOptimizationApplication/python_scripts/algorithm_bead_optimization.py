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
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
from .algorithm_base import OptimizationAlgorithm
from . import mapper_factory
from . import data_logger_factory
from .custom_timer import Timer
from .custom_variable_utilities import WriteDictionaryDataOnNodalVariable
import math

# ==============================================================================
class AlgorithmBeadOptimization(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = KM.Parameters("""
        {
            "name"                          : "bead_optimization",
            "bead_height"                   : 1.0,
            "bead_direction"                : [],
            "bead_side"                     : "both",
            "fix_boundaries"                : [],
            "estimated_lagrange_multiplier" : 1.0,
            "max_total_iterations"          : 10000,
            "max_outer_iterations"          : 10000,
            "max_inner_iterations"          : 30,
            "min_inner_iterations"          : 3,
            "inner_iteration_tolerance"     : 1e-3,
            "line_search" : {
                "line_search_type"           : "manual_stepping",
                "normalize_search_direction" : true,
                "step_size"                  : 1.0
            },
            "filter_penalty_term"           : false,
            "penalty_filter_radius"         : -1.0
        }""")
        self.algorithm_settings =  optimization_settings["optimization_algorithm"]
        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

        self.optimization_settings = optimization_settings
        self.mapper_settings = optimization_settings["design_variables"]["filter"]

        if self.algorithm_settings["filter_penalty_term"].GetBool():
            if self.algorithm_settings["penalty_filter_radius"].GetDouble() == -1.0:
                raise RuntimeError("The parameter `penalty_filter_radius` is missing in order to filter the penalty term!")

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.design_surface = None
        self.mapper = None
        self.penalty_filter = None
        self.data_logger = None
        self.optimization_utilities = None

        self.objectives = optimization_settings["objectives"]
        self.constraints = optimization_settings["constraints"]

        self.bead_height = self.algorithm_settings["bead_height"].GetDouble()
        self.bead_side = self.algorithm_settings["bead_side"].GetString()
        self.filter_penalty_term = self.algorithm_settings["filter_penalty_term"].GetBool()
        self.estimated_lagrange_multiplier = self.algorithm_settings["estimated_lagrange_multiplier"].GetDouble()
        self.max_total_iterations = self.algorithm_settings["max_total_iterations"].GetInt()
        self.max_outer_iterations = self.algorithm_settings["max_outer_iterations"].GetInt()
        self.max_inner_iterations = self.algorithm_settings["max_inner_iterations"].GetInt()
        self.min_inner_iterations = self.algorithm_settings["min_inner_iterations"].GetInt()
        self.inner_iteration_tolerance = self.algorithm_settings["inner_iteration_tolerance"].GetDouble()
        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()

        self.lower_bound = None
        self.upper_bound = None

        self.lambda0 = 0.0
        self.penalty_scaling_0 = 1.0
        self.penalty_factor_0 = 1.0

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.ALPHA)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.ALPHA_MAPPED)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.DF1DALPHA)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.DF1DALPHA_MAPPED)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.DPDALPHA)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.DPDALPHA_MAPPED)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.DLDALPHA)

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("The augmented lagrange algorithm for bead optimization only supports one objective function!")
        if self.constraints.size() > 0:
            raise RuntimeError("The augmented lagrange algorithm for bead does not allow for any constraints!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.model_part_controller.Initialize()
        self.model_part_controller.SetMinimalBufferSize(2)

        self.analyzer.InitializeBeforeOptimizationLoop()

        self.design_surface = self.model_part_controller.GetDesignSurface()

        self.mapper = mapper_factory.CreateMapper(self.design_surface, self.design_surface, self.mapper_settings)
        self.mapper.Initialize()

        if self.filter_penalty_term:
            penalty_filter_radius = self.algorithm_settings["penalty_filter_radius"].GetDouble()
            filter_radius = self.mapper_settings["filter_radius"].GetDouble()
            if abs(filter_radius - penalty_filter_radius) > 1e-9:
                penalty_filter_settings = self.mapper_settings.Clone()
                penalty_filter_settings["filter_radius"].SetDouble(self.algorithm_settings["penalty_filter_radius"].GetDouble())
                self.penalty_filter = mapper_factory.CreateMapper(self.design_surface, self.design_surface, penalty_filter_settings)
                self.penalty_filter.Initialize()
            else:
                self.penalty_filter = self.mapper

        self.data_logger = data_logger_factory.CreateDataLogger(self.model_part_controller, self.communicator, self.optimization_settings)
        self.data_logger.InitializeDataLogging()

        self.optimization_utilities = KSO.OptimizationUtilities(self.design_surface, self.optimization_settings)

        # Identify fixed design areas
        KM.VariableUtils().SetFlag(KM.BOUNDARY, False, self.optimization_model_part.Nodes)

        radius = self.mapper_settings["filter_radius"].GetDouble()
        search_based_functions = KSO.SearchBasedFunctions(self.design_surface)

        for itr in range(self.algorithm_settings["fix_boundaries"].size()):
            sub_model_part_name = self.algorithm_settings["fix_boundaries"][itr].GetString()
            node_set = self.optimization_model_part.GetSubModelPart(sub_model_part_name).Nodes
            search_based_functions.FlagNodesInRadius(node_set, KM.BOUNDARY, radius)

        # Specify bounds and assign starting values for ALPHA
        if self.bead_side == "positive":
            KM.VariableUtils().SetScalarVar(KSO.ALPHA, 0.5, self.design_surface.Nodes, KM.BOUNDARY, False)
            self.lower_bound = 0.0
            self.upper_bound = 1.0
        elif self.bead_side == "negative":
            KM.VariableUtils().SetScalarVar(KSO.ALPHA, -0.5, self.design_surface.Nodes, KM.BOUNDARY, False)
            self.lower_bound = -1.0
            self.upper_bound = 0.0
        elif self.bead_side == "both":
            KM.VariableUtils().SetScalarVar(KSO.ALPHA, 0.0, self.design_surface.Nodes, KM.BOUNDARY, False)
            self.lower_bound = -1.0
            self.upper_bound = 1.0
        else:
            raise RuntimeError("Specified bead direction mode not supported!")

        # Initialize ALPHA_MAPPED according to initial ALPHA values
        self.mapper.Map(KSO.ALPHA, KSO.ALPHA_MAPPED)

        # Specify bead direction
        bead_direction = self.algorithm_settings["bead_direction"].GetVector()
        if len(bead_direction) == 0:
            self.model_part_controller.ComputeUnitSurfaceNormals()
            for node in self.design_surface.Nodes:
                normalized_normal = node.GetSolutionStepValue(KSO.NORMALIZED_SURFACE_NORMAL)
                node.SetValue(KSO.BEAD_DIRECTION,normalized_normal)

        elif len(bead_direction) == 3:
            norm = math.sqrt(bead_direction[0]**2 + bead_direction[1]**2 + bead_direction[2]**2)
            normalized_bead_direction = [value/norm for value in bead_direction]
            KM.VariableUtils().SetNonHistoricalVectorVar(KSO.BEAD_DIRECTION, normalized_bead_direction, self.design_surface.Nodes)
        else:
            raise RuntimeError("Wrong definition of bead direction. Options are: 1) [] -> takes surface normal, 2) [x.x,x.x,x.x] -> takes specified vector.")

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        current_lambda = self.lambda0
        penalty_scaling = self.penalty_scaling_0
        penalty_factor = self.penalty_factor_0

        total_iteration = 0
        is_design_converged = False
        is_max_total_iterations_reached = False
        previos_L = None

        for outer_iteration in range(1,self.max_outer_iterations+1):
            for inner_iteration in range(1,self.max_inner_iterations+1):

                total_iteration += 1
                timer.StartNewLap()

                print("\n>=======================================================================================")
                print("> ",timer.GetTimeStamp(),": Starting iteration ",outer_iteration,".",inner_iteration,".",total_iteration,"(outer . inner . total)")
                print(">=======================================================================================\n")

                # Initialize new shape
                self.model_part_controller.UpdateTimeStep(total_iteration)

                for node in self.design_surface.Nodes:
                    new_shape_change = node.GetSolutionStepValue(KSO.ALPHA_MAPPED) * node.GetValue(KSO.BEAD_DIRECTION) * self.bead_height
                    node.SetSolutionStepValue(KSO.SHAPE_CHANGE, new_shape_change)

                self.model_part_controller.DampNodalVariableIfSpecified(KSO.SHAPE_CHANGE)

                for node in self.design_surface.Nodes:
                    shape_update = node.GetSolutionStepValue(KSO.SHAPE_CHANGE,0) - node.GetSolutionStepValue(KSO.SHAPE_CHANGE,1)
                    node.SetSolutionStepValue(KSO.SHAPE_UPDATE, shape_update)

                self.model_part_controller.UpdateMeshAccordingInputVariable(KSO.SHAPE_UPDATE)
                self.model_part_controller.SetReferenceMeshToMesh()

                # Analyze shape
                self.communicator.initializeCommunication()
                self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
                self.communicator.requestGradientOf(self.objectives[0]["identifier"].GetString())

                self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, total_iteration, self.communicator)

                objective_value = self.communicator.getStandardizedValue(self.objectives[0]["identifier"].GetString())
                objGradientDict = self.communicator.getStandardizedGradient(self.objectives[0]["identifier"].GetString())
                WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, KSO.DF1DX)

                self.model_part_controller.DampNodalVariableIfSpecified(KSO.DF1DX)

                # Compute sensitivities w.r.t. scalar design variable alpha
                for node in self.design_surface.Nodes:
                    raw_gradient = node.GetSolutionStepValue(KSO.DF1DX)
                    bead_dir = node.GetValue(KSO.BEAD_DIRECTION)

                    dF1dalpha_i = self.bead_height*(raw_gradient[0]*bead_dir[0] + raw_gradient[1]*bead_dir[1] + raw_gradient[2]*bead_dir[2])
                    node.SetSolutionStepValue(KSO.DF1DALPHA, dF1dalpha_i)

                # Map gradient of objective
                self.mapper.InverseMap(KSO.DF1DALPHA, KSO.DF1DALPHA_MAPPED)

                # Compute scaling
                max_norm_objective_gradient = self.optimization_utilities.ComputeMaxNormOfNodalVariable(KSO.DF1DALPHA_MAPPED)

                if outer_iteration == 1 and inner_iteration == min(3,self.max_inner_iterations):
                    if self.bead_side == "positive" or self.bead_side == "negative":
                        max_norm_penalty_gradient = 1.0
                    elif self.bead_side == "both":
                        max_norm_penalty_gradient = 2.0

                    penalty_scaling = max_norm_objective_gradient/max_norm_penalty_gradient

                # Compute penalization term
                penalty_value = 0.0
                if self.bead_side == "positive":
                    for node in self.design_surface.Nodes:
                        if not node.Is(KM.BOUNDARY):
                            alpha_i = node.GetSolutionStepValue(KSO.ALPHA)
                            penalty_value += penalty_scaling*(alpha_i-alpha_i**2)

                            penalty_gradient_i = penalty_scaling*(1-2*alpha_i)
                            node.SetSolutionStepValue(KSO.DPDALPHA, penalty_gradient_i)

                elif self.bead_side == "negative":
                    for node in self.design_surface.Nodes:
                        if not node.Is(KM.BOUNDARY):
                            alpha_i = node.GetSolutionStepValue(KSO.ALPHA)
                            penalty_value += penalty_scaling*(-alpha_i-alpha_i**2)

                            penalty_gradient_i = penalty_scaling*(-1-2*alpha_i)
                            node.SetSolutionStepValue(KSO.DPDALPHA, penalty_gradient_i)

                elif self.bead_side == "both":
                    for node in self.design_surface.Nodes:
                        if not node.Is(KM.BOUNDARY):
                            alpha_i = node.GetSolutionStepValue(KSO.ALPHA)
                            penalty_value += penalty_scaling*(-alpha_i**2+1)

                            penalty_gradient_i = penalty_scaling*(-2*alpha_i)
                            node.SetSolutionStepValue(KSO.DPDALPHA, penalty_gradient_i)

                # Filter penalty term if specified
                if self.filter_penalty_term:
                    self.penalty_filter.InverseMap(KSO.DPDALPHA, KSO.DPDALPHA_MAPPED)

                # Compute value of Lagrange function
                L = objective_value + current_lambda*penalty_value + 0.5*penalty_factor*penalty_value**2
                if inner_iteration == 1:
                    dL_relative = 0.0
                else:
                    dL_relative = 100*(L/previos_L-1)

                # Compute gradient of Lagrange function
                if self.filter_penalty_term:
                    penalty_gradient_variable = KSO.DPDALPHA_MAPPED
                else:
                    penalty_gradient_variable = KSO.DPDALPHA
                for node in self.design_surface.Nodes:
                    dLdalpha_i = node.GetSolutionStepValue(KSO.DF1DALPHA_MAPPED) + current_lambda*node.GetSolutionStepValue(penalty_gradient_variable)
                    node.SetSolutionStepValue(KSO.DLDALPHA, dLdalpha_i)

                # Normalization using infinity norm
                dLdalpha_for_normalization = {}
                for node in self.design_surface.Nodes:
                    nodal_alpha = node.GetSolutionStepValue(KSO.ALPHA)
                    if nodal_alpha==self.lower_bound or nodal_alpha==self.upper_bound or node.Is(KM.BOUNDARY):
                        dLdalpha_for_normalization[node.Id] = 0.0
                    else:
                        dLdalpha_for_normalization[node.Id] = node.GetSolutionStepValue(KSO.DLDALPHA)**2

                max_value = math.sqrt(max(dLdalpha_for_normalization.values()))
                if max_value == 0.0:
                    max_value = 1.0

                # Compute updated design variable
                for node in self.design_surface.Nodes:
                    dalpha = -self.step_size*node.GetSolutionStepValue(KSO.DLDALPHA)/max_value
                    alpha_new = node.GetSolutionStepValue(KSO.ALPHA) + dalpha

                    # Enforce bounds
                    alpha_new = max(alpha_new, self.lower_bound)
                    alpha_new = min(alpha_new, self.upper_bound)

                    # Enforce constraints
                    if node.Is(KM.BOUNDARY):
                        alpha_new = 0.0

                    node.SetSolutionStepValue(KSO.ALPHA,alpha_new)

                    alpha_new_vectorized = alpha_new * node.GetValue(KSO.BEAD_DIRECTION)
                    node.SetSolutionStepValue(KSO.CONTROL_POINT_CHANGE,alpha_new_vectorized)

                # Map design variables
                self.mapper.Map(KSO.ALPHA, KSO.ALPHA_MAPPED)

                # Log current optimization step and store values for next iteration
                additional_values_to_log = {}
                additional_values_to_log["step_size"] = self.algorithm_settings["line_search"]["step_size"].GetDouble()
                additional_values_to_log["outer_iteration"] = outer_iteration
                additional_values_to_log["inner_iteration"] = inner_iteration
                additional_values_to_log["lagrange_value"] = L
                additional_values_to_log["lagrange_value_relative_change"] = dL_relative
                additional_values_to_log["penalty_value"] = penalty_value
                additional_values_to_log["penalty_lambda"] = current_lambda
                additional_values_to_log["penalty_scaling"] = penalty_scaling
                additional_values_to_log["penalty_factor"] = penalty_factor
                additional_values_to_log["max_norm_objective_gradient"] = max_norm_objective_gradient

                self.data_logger.LogCurrentValues(total_iteration, additional_values_to_log)
                self.data_logger.LogCurrentDesign(total_iteration)

                previos_L = L

                # Convergence check of inner loop
                if total_iteration == self.max_total_iterations:
                    is_max_total_iterations_reached = True
                    break

                if inner_iteration >= self.min_inner_iterations and inner_iteration >1:
                    # In the first outer iteration, the constraint is not yet active and properly scaled. Therefore, the objective is used to check the relative improvement
                    if outer_iteration == 1:
                        if abs(self.data_logger.GetValues("rel_change_objective")[total_iteration]) < self.inner_iteration_tolerance:
                            break
                    else:
                        if abs(dL_relative) < self.inner_iteration_tolerance:
                            break

                if penalty_value == 0.0:
                    is_design_converged = True
                    break

                print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
                print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

            # Compute penalty factor such that estimated Lagrange multiplier is obtained
            if outer_iteration==1:
                penalty_factor = self.estimated_lagrange_multiplier/penalty_value

            # Update lambda
            current_lambda = current_lambda + penalty_factor*penalty_value

            print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
            print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

            # Check convergence of outer loop
            if outer_iteration == self.max_outer_iterations:
                print("\n> Maximal outer iterations of optimization problem reached!")
                break

            if is_max_total_iterations_reached:
                print("\n> Maximal total iterations of optimization problem reached!")
                break

            if is_design_converged:
                print("\n> Update of design variables is zero. Optimization converged!")
                break

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.data_logger.FinalizeDataLogging()
        self.analyzer.FinalizeAfterOptimizationLoop()

# ==============================================================================
