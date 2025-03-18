# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Schmoelz David, https://github.com/dschmoelz
#
#
# ==============================================================================


# Kratos Core and Apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.LinearSolversApplication import dense_linear_solver_factory

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication.algorithms.algorithm_base import OptimizationAlgorithm
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
from KratosMultiphysics.ShapeOptimizationApplication.loggers import data_logger_factory
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer

import math

# ==============================================================================
class AlgorithmFreeThicknessOptimization(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = KM.Parameters("""
        {
            "name"                    : "free_thickness_optimization",
            "max_correction_share"    : 0.75,
            "max_iterations"          : 100,
            "relative_tolerance"      : 1e-3,
            "line_search" : {
                "line_search_type"           : "manual_stepping",
                "normalize_search_direction" : true,
                "step_size"                  : 1.0
            }
        }""")
        self.algorithm_settings =  optimization_settings["optimization_algorithm"]
        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

        self.optimization_settings = optimization_settings
        self.mapper_settings = optimization_settings["design_variables"][0]["filter"]

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.design_surface = None
        self.mapper = None
        self.data_logger = None
        self.optimization_utilities = None

        self.objectives = optimization_settings["objectives"]
        self.constraints = optimization_settings["constraints"]
        self.constraint_gradient_variables = {}
        for itr, constraint in enumerate(self.constraints):
            self.constraint_gradient_variables.update({
                constraint["identifier"].GetString() : {
                    "gradient": KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT"),
                    "mapped_gradient": KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT_MAPPED")
                }
            })
        self.max_correction_share = self.algorithm_settings["max_correction_share"].GetDouble()

        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_SEARCH_DIRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CORRECTION)

        self.lower_bound = LowerBound(optimization_settings["design_variables"][0]["t_min"].GetDouble(), self.optimization_model_part)
        self.upper_bound = UpperBound(optimization_settings["design_variables"][0]["t_max"].GetDouble(), self.optimization_model_part)

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Gradient projection algorithm only supports one objective function!")
        if self.constraints.size() == 0:
            raise RuntimeError("Gradient projection algorithm requires definition of at least one constraint!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.model_part_controller.Initialize()

        self.analyzer.InitializeBeforeOptimizationLoop()

        self.model_part_controller.ModifyInitialProperties()

        self.design_surface = self.model_part_controller.GetDesignSurface()

        self.mapper = mapper_factory.CreateMapper(self.design_surface, self.design_surface, self.mapper_settings)
        self.mapper.Initialize()
        self.model_part_controller.InitializeDamping()

        self.data_logger = data_logger_factory.CreateDataLogger(self.model_part_controller, self.communicator, self.optimization_settings)
        self.data_logger.InitializeDataLogging()

        self.optimization_utilities = KSO.OptimizationUtilities

        self.__InitializeThicknessField()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        for self.optimization_iteration in range(1,self.max_iterations):
            KM.Logger.Print("")
            KM.Logger.Print("===============================================================================")
            KM.Logger.PrintInfo("ShapeOpt", timer.GetTimeStamp(), ": Starting optimization iteration ", self.optimization_iteration)
            KM.Logger.Print("===============================================================================\n")

            timer.StartNewLap()

            self.__initializeNewThickness()

            self.__analyzeThickness()

            self.__computeThicknessUpdate()

            self.__logCurrentOptimizationStep()

            KM.Logger.Print("")
            KM.Logger.PrintInfo("ShapeOpt", "Time needed for current optimization step = ", timer.GetLapTime(), "s")
            KM.Logger.PrintInfo("ShapeOpt", "Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

            if self.__isAlgorithmConverged():
                break
            else:
                self.__determineAbsoluteChanges()

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.data_logger.FinalizeDataLogging()
        self.analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeNewThickness(self):
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)
        self.model_part_controller.UpdateThicknessAccordingInputVariable(KSO.THICKNESS_UPDATE)

        self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyzeThickness(self):
        self.communicator.initializeCommunication()
        self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
        self.communicator.requestThicknessGradientOf(self.objectives[0]["identifier"].GetString())

        for constraint in self.constraints:
            con_id =  constraint["identifier"].GetString()
            self.communicator.requestValueOf(con_id)
            self.communicator.requestThicknessGradientOf(con_id)

        self.lower_bound.CalculateValueAndGradient()
        self.upper_bound.CalculateValueAndGradient()

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.optimization_iteration, self.communicator)

        # project and damp objective gradients
        objElementGradientDict = self.communicator.getStandardizedThicknessGradient(self.objectives[0]["identifier"].GetString())
        self.__mapElementGradientToNode(objElementGradientDict, KSO.DF1DT)

        self.model_part_controller.DampNodalSensitivityVariableIfSpecified(KSO.DF1DT)

        # project and damp constraint gradients
        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            conElementGradientDict = self.communicator.getStandardizedThicknessGradient(con_id)
            gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
            self.__mapElementGradientToNode(conElementGradientDict, gradient_variable)

            self.model_part_controller.DampNodalSensitivityVariableIfSpecified(gradient_variable)

        self.model_part_controller.DampNodalSensitivityVariableIfSpecified(self.lower_bound.gradient_variable)
        self.model_part_controller.DampNodalSensitivityVariableIfSpecified(self.upper_bound.gradient_variable)

    # --------------------------------------------------------------------------
    def __computeThicknessUpdate(self):
        self.mapper.Update()
        self.mapper.InverseMap(KSO.DF1DT, KSO.DF1DT_MAPPED)

        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
            mapped_gradient_variable = self.constraint_gradient_variables[con_id]["mapped_gradient"]
            self.mapper.InverseMap(gradient_variable, mapped_gradient_variable)

        self.mapper.InverseMap(self.lower_bound.gradient_variable, self.lower_bound.mapped_gradient_variable)
        self.mapper.InverseMap(self.upper_bound.gradient_variable, self.upper_bound.mapped_gradient_variable)

        self.__computeControlThicknessUpdate()

        self.mapper.Map(KSO.THICKNESS_CONTROL_UPDATE, KSO.THICKNESS_UPDATE)
        self.model_part_controller.DampNodalUpdateVariableIfSpecified(KSO.THICKNESS_UPDATE)

    # --------------------------------------------------------------------------
    def __computeControlThicknessUpdate(self):
        """adapted from https://msulaiman.org/onewebmedia/GradProj_2.pdf"""
        g_a, g_a_variables = self.__getActiveConstraints()

        KM.Logger.PrintInfo("ShapeOpt", "Assemble vector of objective gradient.")
        nabla_f = KM.Vector()
        s = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f, KSO.DF1DT_MAPPED)

        if len(g_a) == 0:
            KM.Logger.PrintInfo("ShapeOpt", "No constraints active, use negative objective gradient as search direction.")
            s = nabla_f * (-1.0)
            s *= self.step_size / s.norm_inf()
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_SEARCH_DIRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, [0.0]*len(s), KSO.THICKNESS_CORRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_CONTROL_UPDATE)
            return


        KM.Logger.PrintInfo("ShapeOpt", "Assemble matrix of constraint gradient.")
        N = KM.Matrix()
        self.optimization_utilities.AssembleMatrixScalarVariables(self.design_surface, N, g_a_variables)

        settings = KM.Parameters('{ "solver_type" : "LinearSolversApplication.dense_col_piv_householder_qr" }')
        solver = dense_linear_solver_factory.ConstructSolver(settings)

        KM.Logger.PrintInfo("ShapeOpt", "Calculate projected search direction and correction.")
        c = KM.Vector()
        self.optimization_utilities.CalculateProjectedSearchDirectionAndCorrection(
            nabla_f,
            N,
            g_a,
            solver,
            s,
            c)

        if c.norm_inf() != 0.0:
            if c.norm_inf() <= self.max_correction_share * self.step_size:
                delta = self.step_size - c.norm_inf()
                s *= delta/s.norm_inf()
            else:
                KM.Logger.PrintWarning("ShapeOpt", f"Correction is scaled down from {c.norm_inf()} to {self.max_correction_share * self.step_size}.")
                c *= self.max_correction_share * self.step_size / c.norm_inf()
                s *= (1.0 - self.max_correction_share) * self.step_size / s.norm_inf()
        else:
            s *= self.step_size / s.norm_inf()

        self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_SEARCH_DIRECTION)
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, c, KSO.THICKNESS_CORRECTION)
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, s+c, KSO.THICKNESS_CONTROL_UPDATE)

    # --------------------------------------------------------------------------
    def __getActiveConstraints(self):
        active_constraint_values = []
        active_constraint_variables = []

        for constraint in self.constraints:
            if self.__isConstraintActive(constraint):
                identifier = constraint["identifier"].GetString()
                constraint_value = self.communicator.getStandardizedValue(identifier)
                active_constraint_values.append(constraint_value)
                active_constraint_variables.append(
                    self.constraint_gradient_variables[identifier]["mapped_gradient"])

        if self.lower_bound.value > 0:
            active_constraint_values.append(self.lower_bound.value)
            active_constraint_variables.append(
                self.lower_bound.mapped_gradient_variable)

        if self.upper_bound.value > 0:
            active_constraint_values.append(self.upper_bound.value)
            active_constraint_variables.append(
                self.upper_bound.mapped_gradient_variable)

        return active_constraint_values, active_constraint_variables

    # --------------------------------------------------------------------------
    def __isConstraintActive(self, constraint):
        identifier = constraint["identifier"].GetString()
        constraint_value = self.communicator.getStandardizedValue(identifier)
        if constraint["type"].GetString() == "=" or constraint_value >= 0:
            gradient_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(
                self.design_surface, self.constraint_gradient_variables[identifier]["mapped_gradient"]
            )
            if math.isclose(gradient_norm, 0.0, abs_tol=1e-16):
                KM.Logger.PrintWarning("ShapeOpt", f"Gradient for constraint {identifier} is 0.0 - will not be considered!")
                return False
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.step_size
        additional_values_to_log["inf_norm_s"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_SEARCH_DIRECTION)
        additional_values_to_log["inf_norm_c"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_CORRECTION)
        # self.data_logger.LogSensitivityHeatmap(self.optimization_iteration, self.mapper)
        self.data_logger.LogCurrentValues(self.optimization_iteration, additional_values_to_log)
        self.data_logger.LogCurrentDesign(self.optimization_iteration)

    # --------------------------------------------------------------------------
    def __isAlgorithmConverged(self):

        if self.optimization_iteration > 1 :

            # Check if maximum iterations were reached
            if self.optimization_iteration == self.max_iterations:
                KM.Logger.Print("")
                KM.Logger.PrintInfo("ShapeOpt", "Maximal iterations of optimization problem reached!")
                return True

            # Check for relative tolerance
            relative_change_of_objective_value = self.data_logger.GetValues("rel_change_objective")[self.optimization_iteration]
            if abs(relative_change_of_objective_value) < self.relative_tolerance:
                KM.Logger.Print("")
                KM.Logger.PrintInfo("ShapeOpt", "Optimization problem converged within a relative objective tolerance of ",self.relative_tolerance,"%.")
                return True

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges(self):

        for node in self.optimization_model_part.Nodes:
            absolute_control_update = node.GetSolutionStepValue(KSO.THICKNESS_CONTROL_CHANGE)
            absolute_control_update += node.GetSolutionStepValue(KSO.THICKNESS_CONTROL_UPDATE)
            node.SetSolutionStepValue(KSO.THICKNESS_CONTROL_CHANGE, 0, absolute_control_update)

            absolute_update = node.GetSolutionStepValue(KSO.THICKNESS_CHANGE)
            absolute_update += node.GetSolutionStepValue(KSO.THICKNESS_UPDATE)
            node.SetSolutionStepValue(KSO.THICKNESS_CHANGE, 0, absolute_update)

            thickness = node.GetSolutionStepValue(KSO.THICKNESS)
            thickness += node.GetSolutionStepValue(KSO.THICKNESS_UPDATE)
            node.SetSolutionStepValue(KSO.THICKNESS, 0, thickness)

        return

    # --------------------------------------------------------------------------
    def __mapElementGradientToNode(self, element_gradient_dict, gradient_variable):

        # reset variables
        for node in self.design_surface.Nodes:
            node.SetSolutionStepValue(gradient_variable, 0, 0.0)

        for condition in self.design_surface.Conditions:
            condition.SetValue(gradient_variable, element_gradient_dict[condition.Id])

        total_node_areas = dict()
        for condition in self.design_surface.Conditions:
            df_dt = condition.GetValue(gradient_variable)
            for node in condition.GetNodes():
                if node.Id in total_node_areas:
                    total_node_areas[node.Id] += condition.GetGeometry().Area()
                else:
                    total_node_areas[node.Id] = condition.GetGeometry().Area()
                df_dt_node = node.GetSolutionStepValue(gradient_variable)
                df_dt_node += df_dt * condition.GetGeometry().Area()
                node.SetSolutionStepValue(gradient_variable, 0, df_dt_node)

        for node in self.design_surface.Nodes:
            total_node_area = total_node_areas[node.Id]
            df_dt_node = node.GetSolutionStepValue(gradient_variable)
            df_dt_node /= total_node_area
            node.SetSolutionStepValue(gradient_variable, 0, df_dt_node)

    # --------------------------------------------------------------------------
    def __InitializeThicknessField(self):

        element_thicknesses = dict()
        for condition in self.optimization_model_part.Conditions:
            element_thicknesses[condition.Id] = condition.Properties.GetValue(KM.THICKNESS)

        self.__mapElementGradientToNode(element_thicknesses, KSO.THICKNESS)

# ==============================================================================

class LowerBound(object):

    def __init__(self, t_min, model_part):
        self.t_min = t_min
        self.value = 0.0
        self.model_part = model_part

        self.gradient_variable = KM.KratosGlobals.GetVariable(f"DLBDT")
        self.mapped_gradient_variable = KM.KratosGlobals.GetVariable(f"DLBDT_MAPPED")

        self.model_part.AddNodalSolutionStepVariable(self.gradient_variable)
        self.model_part.AddNodalSolutionStepVariable(self.mapped_gradient_variable)

    def CalculateValueAndGradient(self):

        g = 0

        for node in self.model_part.Nodes:

            t_i = node.GetSolutionStepValue(KSO.THICKNESS)
            g_i = self.t_min - t_i
            dg_i_dt = 0

            if g_i > 0:
                g += g_i * g_i
                dg_i_dt = -1 * 2 * g_i

            node.SetSolutionStepValue(self.gradient_variable, 0, dg_i_dt)

        self.value = g

class UpperBound(object):

    def __init__(self, t_max, model_part):
        self.t_max = t_max
        self.value = 0.0
        self.model_part = model_part

        self.gradient_variable = KM.KratosGlobals.GetVariable(f"DUBDT")
        self.mapped_gradient_variable = KM.KratosGlobals.GetVariable(f"DUBDT_MAPPED")

        self.model_part.AddNodalSolutionStepVariable(self.gradient_variable)
        self.model_part.AddNodalSolutionStepVariable(self.mapped_gradient_variable)

    def CalculateValueAndGradient(self):

        g = 0

        for node in self.model_part.Nodes:

            t_i = node.GetSolutionStepValue(KSO.THICKNESS)
            g_i = t_i - self.t_max
            dg_i_dt = 0

            if g_i > 0:
                g += g_i * g_i
                dg_i_dt = 2 * g_i

            node.SetSolutionStepValue(self.gradient_variable, 0, dg_i_dt)

        self.value = g