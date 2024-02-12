# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    David Schmölz
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
import KratosMultiphysics as Kratos
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.LinearSolversApplication import dense_linear_solver_factory

# Additional imports
from .algorithm_base import OptimizationAlgorithm
from .. import mapper_factory
from KratosMultiphysics.ShapeOptimizationApplication.loggers import data_logger_factory
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer
import numpy as np
import sys


# ==============================================================================
class AlgorithmFreeThicknessRelaxedGradientProjection(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = Kratos.Parameters("""
        {
            "name"                    : "free_thickness_rgp",
            "max_iterations"          : 100,
            "max_inner_iter"          : 100,
            "relative_tolerance"      : 1e-3,
            "line_search" : {
                "line_search_type"           : "manual_stepping",
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
        self.constraint_buffer_variables = {}
        for itr, constraint in enumerate(self.constraints):
            constraint_id = constraint["identifier"].GetString()
            self.constraint_gradient_variables.update({
                constraint_id : {
                    "gradient": Kratos.KratosGlobals.GetVariable(f"DC{(itr+1)}DT"),
                    "mapped_gradient": Kratos.KratosGlobals.GetVariable(f"DC{(itr+1)}DT_MAPPED")
                }
            })
            self.constraint_buffer_variables.update({
                constraint_id : {
                    "buffer_value": 0.0,
                    "buffer_value-1": 0.0,
                    "buffer_size": 1e-12,
                    "buffer_size_factor": 2.0,
                    "central_buffer_value": 0.0,
                    "lower_buffer_value": - 1e-12,
                    "upper_buffer_value": 1e-12,
                    "g_i-1": 0.0,
                    "g_i-2": 0.0,
                    "g_i-3": 0.0,
                    "max_constraint_change": 0.0
                }
            })


        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.line_search_type = self.algorithm_settings["line_search"]["line_search_type"].GetString()
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()
        self.s_norm = 0.0

        # disable qnrgp
        self.max_inner_iter = self.algorithm_settings["max_inner_iter"].GetDouble()
        # self.max_inner_iter = 1
        self.buffer_coeff_update = 2.0 / self.max_inner_iter

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.INV_HESSIAN_THICKNESS)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_PROJECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CORRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_SEARCH_DIRECTION)

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
            Kratos.Logger.Print("")
            Kratos.Logger.Print("===============================================================================")
            Kratos.Logger.PrintInfo("ShapeOpt", timer.GetTimeStamp(), ": Starting optimization iteration ", self.optimization_iteration)
            Kratos.Logger.Print("===============================================================================\n")

            timer.StartNewLap()

            self.__initializeNewThickness()

            self.__analyzeThickness()

            self.__computeBufferValue()

            self.__computeThicknessUpdate()

            self.__logCurrentOptimizationStep()

            self.__updateBufferZone()

            Kratos.Logger.Print("")
            Kratos.Logger.PrintInfo("ShapeOpt", "Time needed for current optimization step = ", timer.GetLapTime(), "s")
            Kratos.Logger.PrintInfo("ShapeOpt", "Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

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
        self.lower_bound.AggregateNodalConstraints(self.optimization_utilities, self.design_surface)
        self.upper_bound.CalculateValueAndGradient()
        self.upper_bound.AggregateNodalConstraints(self.optimization_utilities, self.design_surface)

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
    def __computeBufferValue(self):
        # compute new buffer size and buffer values
        for constraint in self.constraints:
            identifier = constraint["identifier"].GetString()
            g_i = self.communicator.getStandardizedValue(identifier)
            g_i_m1 = self.constraint_buffer_variables[identifier]["g_i-1"]
            buffer_size_factor = self.constraint_buffer_variables[identifier]["buffer_size_factor"]

            self.constraint_buffer_variables[identifier]["buffer_value-1"] = self.constraint_buffer_variables[identifier]["buffer_value"]

            if self.optimization_iteration > 1:

                if abs(g_i - g_i_m1) > self.constraint_buffer_variables[identifier]["max_constraint_change"]:
                    self.constraint_buffer_variables[identifier]["max_constraint_change"] = abs(g_i - g_i_m1)

                max_constraint_change = self.constraint_buffer_variables[identifier]["max_constraint_change"]
                self.constraint_buffer_variables[identifier]["buffer_size"] = max(buffer_size_factor * max_constraint_change, 1e-12)

            buffer_size = self.constraint_buffer_variables[identifier]["buffer_size"]
            self.constraint_buffer_variables[identifier]["lower_buffer_value"] = self.constraint_buffer_variables[identifier]["central_buffer_value"] \
                - buffer_size
            self.constraint_buffer_variables[identifier]["upper_buffer_value"] = self.constraint_buffer_variables[identifier]["central_buffer_value"] \
                + buffer_size

            if self.__isConstraintActive(constraint):
                if constraint["type"].GetString() == "=":
                    self.constraint_buffer_variables[identifier]["buffer_value"] = min(1 - abs(g_i) / buffer_size, 2.0)
                else:
                    lower_buffer_value = self.constraint_buffer_variables[identifier]["lower_buffer_value"]
                    self.constraint_buffer_variables[identifier]["buffer_value"] = min( (g_i - lower_buffer_value) / buffer_size, 2.0 )
            else:
                self.constraint_buffer_variables[identifier]["buffer_value"] = 0.0

        self.lower_bound.ComputeBufferValue(self.optimization_iteration)
        self.upper_bound.ComputeBufferValue(self.optimization_iteration)

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

        self.inner_iter = 1

        while not self.__checkInnerConvergence():

            self.direction_has_changed = False

            Kratos.Logger.PrintInfo("ShapeOpt", "Inner Iteration to Find Shape Update = ", self.inner_iter)

            self.__computeControlThicknessUpdate()

            self.mapper.Map(KSO.THICKNESS_CONTROL_UPDATE, KSO.THICKNESS_UPDATE)
            self.model_part_controller.DampNodalUpdateVariableIfSpecified(KSO.THICKNESS_UPDATE)

            self.d_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_UPDATE)

            self.__checkConstraintValue()
            self.inner_iter += 1

        self.__saveLineSearchData()


    # --------------------------------------------------------------------------

    def __checkInnerConvergence(self):
        Kratos.Logger.PrintInfo("Check Convergence of the inner loop:")
        if self.inner_iter == 1:
            return False
        elif self.direction_has_changed and self.inner_iter <= self.max_inner_iter:
            return False
        else:
            return True

    def __checkConstraintValue(self):
        index = -1
        for constraint in self.constraints:
            if self.__isConstraintActive(constraint):
                index += 1
                identifier = constraint["identifier"].GetString()
                g_i = self.communicator.getStandardizedValue(identifier)
                g_a_variable = self.constraint_gradient_variables[identifier]["gradient"]
                shape_update = Kratos.Vector()
                gradient = Kratos.Vector()
                self.optimization_utilities.AssembleVector(self.design_surface, gradient, g_a_variable)
                self.optimization_utilities.AssembleVector(self.design_surface, shape_update, KSO.THICKNESS_UPDATE)
                new_g_i = g_i + np.dot(gradient, shape_update)
                Kratos.Logger.PrintInfo("Constraint ", identifier, "\n Linearized new value = ", new_g_i)
                if new_g_i > 0.0:
                    if self.relaxation_coefficients[index] < 1.0:
                        self.relaxation_coefficients[index] = min(self.relaxation_coefficients[index] + self.buffer_coeff_update, 1.0)
                        self.direction_has_changed = True
                    elif self.correction_coefficients[index] < 2.0:
                        self.correction_coefficients[index] = min (self.correction_coefficients[index] + self.buffer_coeff_update, 2.0)
                        self.direction_has_changed = True
                Kratos.Logger.PrintInfo("Constraint ", identifier, "\n W_R, W_C = ", self.relaxation_coefficients[index], self.correction_coefficients[index])

        direction_has_changed = self.lower_bound.checkConstraintValue(self.buffer_coeff_update, self.optimization_utilities, self.design_surface)
        if direction_has_changed:
            self.relaxation_coefficients[-2] = self.lower_bound.relaxation_coefficient
            self.correction_coefficients[-2] = self.lower_bound.correction_coefficient
            self.direction_has_changed = True

        direction_has_changed = self.upper_bound.checkConstraintValue(self.buffer_coeff_update, self.optimization_utilities, self.design_surface)
        if direction_has_changed:
            self.relaxation_coefficients[-1] = self.upper_bound.relaxation_coefficient
            self.correction_coefficients[-1] = self.upper_bound.correction_coefficient
            self.direction_has_changed = True

    # --------------------------------------------------------------------------
    def __LineSearch(self):
        Kratos.Logger.PrintInfo("Line Search ...")
        if self.line_search_type == "manual_stepping":
            self.__manualStep()
        elif self.line_search_type == "QNBB_method":
            if self.optimization_iteration == 1:
                self.max_step_size = self.step_size
                # Do initial small step
                self.step_size /= 5
                self.__manualStep()
            else:
                self.__QNBBStep()
        elif self.line_search_type == "BB_method":
            if self.optimization_iteration == 1:
                self.max_step_size = self.step_size
                # Do initial small step
                self.step_size /= 5
                self.__manualStep()
            else:
                self.__BBStep()


    def __manualStep(self):
        step_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_CONTROL_UPDATE)
        if abs(step_norm) > 1e-10:
            step = Kratos.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, step, KSO.THICKNESS_CONTROL_UPDATE)
            step *= 1.0 / step_norm
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, step, KSO.THICKNESS_SEARCH_DIRECTION)
            step *= self.step_size
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, step, KSO.THICKNESS_CONTROL_UPDATE)

    def __QNBBStep(self):
        self.s_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_SEARCH_DIRECTION)
        if abs(self.s_norm) > 1e-10:
            s = Kratos.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, s, KSO.THICKNESS_SEARCH_DIRECTION)
            s /= self.s_norm
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_SEARCH_DIRECTION)

            for index, node in enumerate(self.design_surface.Nodes):
                y_i = self.prev_s[index] - s[index]
                d_i = self.d[index]
                if (y_i * y_i) < 1e-9:
                    step_i = self.max_step_size
                else:
                    step_i = abs((d_i * y_i) / (y_i * y_i))
                if step_i > self.max_step_size:
                    step_i = self.max_step_size
                node.SetSolutionStepValue(KSO.INV_HESSIAN_THICKNESS, step_i)
                s[index] = s[index] * step_i
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_CONTROL_UPDATE)

    def __BBStep(self):
        self.s_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_SEARCH_DIRECTION)
        if abs(self.s_norm) > 1e-10:
            s = Kratos.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, s, KSO.THICKNESS_SEARCH_DIRECTION)
            s /= self.s_norm
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_SEARCH_DIRECTION)
            y = self.prev_s - s
            if np.dot(y, y) < 1e-9:
                step = self.max_step_size
            else:
                step = abs(np.dot(y, self.d) / np.dot(y, y))
            if step > self.max_step_size:
                step = self.max_step_size
            s = s * step
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_CONTROL_UPDATE)

    def __saveLineSearchData(self):
        self.prev_s = Kratos.Vector()
        self.d = Kratos.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, self.d, KSO.THICKNESS_CONTROL_UPDATE)
        self.optimization_utilities.AssembleVector(self.design_surface, self.prev_s, KSO.THICKNESS_SEARCH_DIRECTION)

    # --------------------------------------------------------------------------
    def __computeControlThicknessUpdate(self):
        """adapted from https://msulaiman.org/onewebmedia/GradProj_2.pdf"""
        if self.inner_iter == 1:
            self.g_a, self.g_a_variables, self.relaxation_coefficients, self.correction_coefficients = self.__getActiveConstraints()

        Kratos.Logger.PrintInfo("ShapeOpt", "Assemble vector of objective gradient.")
        nabla_f = Kratos.Vector()
        p = Kratos.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f, KSO.DF1DT_MAPPED)
        f_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.DF1DT_MAPPED)

        if abs(f_norm) > 1e-10:
            nabla_f *= 1.0/f_norm

        if len(self.g_a) == 0:
            Kratos.Logger.PrintInfo("ShapeOpt", "No constraints active, use negative objective gradient as search direction.")
            p = nabla_f * (-1.0)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, p, KSO.THICKNESS_SEARCH_DIRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, p, KSO.THICKNESS_PROJECTION)
            Kratos.VariableUtils().SetHistoricalVariableToZero(KSO.THICKNESS_CORRECTION, self.design_surface.Nodes)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, p, KSO.THICKNESS_CONTROL_UPDATE)
            self.__LineSearch()
            return

        omega_r = Kratos.Matrix()
        self.optimization_utilities.AssembleBufferMatrix(omega_r, self.relaxation_coefficients)
        omega_c = Kratos.Vector(self.correction_coefficients)

        Kratos.Logger.PrintInfo("ShapeOpt", "Assemble matrix of constraint gradient.")
        N = Kratos.Matrix()
        self.optimization_utilities.AssembleMatrixScalarVariables(self.design_surface, N, self.g_a_variables)

        settings = Kratos.Parameters('{ "solver_type" : "LinearSolversApplication.dense_col_piv_householder_qr" }')
        solver = dense_linear_solver_factory.ConstructSolver(settings)

        c = Kratos.Vector()

        Kratos.Logger.PrintInfo("ShapeOpt", "Calculate projected search direction and correction.")
        self.optimization_utilities.CalculateRelaxedProjectedSearchDirectionAndCorrection(
            nabla_f,
            N,
            omega_r,
            omega_c,
            solver,
            p,
            c)

        # additional normalization step

        self.optimization_utilities.AssignVectorToVariable(self.design_surface, p, KSO.THICKNESS_PROJECTION)
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, c, KSO.THICKNESS_CORRECTION)
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, p+c, KSO.THICKNESS_SEARCH_DIRECTION)
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, p+c, KSO.THICKNESS_CONTROL_UPDATE)
        self.__LineSearch()

    # --------------------------------------------------------------------------
    def __getActiveConstraints(self):
        active_constraint_values = []
        active_constraint_variables = []
        active_relaxation_coefficient = []
        active_correction_coefficient = []

        for constraint in self.constraints:
            if self.__isConstraintActive(constraint):
                identifier = constraint["identifier"].GetString()
                g_i = self.communicator.getStandardizedValue(identifier)
                buffer_value = self.constraint_buffer_variables[identifier]["buffer_value"]

                active_constraint_values.append(g_i)
                g_a_variable = self.constraint_gradient_variables[identifier]["mapped_gradient"]
                g_a_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, g_a_variable)
                g_a_variable_vector = Kratos.Vector()
                self.optimization_utilities.AssembleVector(self.design_surface, g_a_variable_vector, g_a_variable)
                if abs(g_a_norm) > 1e-10:
                    g_a_variable_vector /= g_a_norm

                self.optimization_utilities.AssignVectorToVariable(self.design_surface, g_a_variable_vector, g_a_variable)

                active_constraint_variables.append(g_a_variable)
                active_relaxation_coefficient.append(min(buffer_value,1.0))

                max_buffer = 2.0
                if buffer_value > 1.0:
                    if buffer_value < max_buffer:
                        active_correction_coefficient.append(2*(buffer_value - 1))
                    else:
                        active_correction_coefficient.append(2*(max_buffer-1))
                else:
                    active_correction_coefficient.append(0.0)

        if self.lower_bound.isActive():
            self.lower_bound.getActive(self.optimization_utilities, self.design_surface)
            active_constraint_values.append(self.lower_bound.value)
            active_constraint_variables.append(self.lower_bound.mapped_gradient_variable)
            active_relaxation_coefficient.append(self.lower_bound.relaxation_coefficient)
            active_correction_coefficient.append(self.lower_bound.correction_coefficient)

        if self.upper_bound.isActive():
            self.upper_bound.getActive(self.optimization_utilities, self.design_surface)
            active_constraint_values.append(self.upper_bound.value)
            active_constraint_variables.append(self.upper_bound.mapped_gradient_variable)
            active_relaxation_coefficient.append(self.upper_bound.relaxation_coefficient)
            active_correction_coefficient.append(self.upper_bound.correction_coefficient)

        return active_constraint_values, active_constraint_variables, active_relaxation_coefficient, active_correction_coefficient

    # --------------------------------------------------------------------------
    def __isConstraintActive(self, constraint):
        identifier = constraint["identifier"].GetString()
        g_i = self.communicator.getStandardizedValue(identifier)
        if constraint["type"].GetString() == "=":
            return True
        elif g_i >= self.constraint_buffer_variables[identifier]["lower_buffer_value"]:
            return True
        else:
            return False
    # --------------------------------------------------------------------------
    def __updateBufferZone(self):
        # adapt the buffer zones for zig-zagging, too much or too little correction
        for constraint in self.constraints:

            identifier = constraint["identifier"].GetString()
            g_i = self.communicator.getStandardizedValue(identifier)
            g_i_m1 = self.constraint_buffer_variables[identifier]["g_i-1"]
            g_i_m2 = self.constraint_buffer_variables[identifier]["g_i-2"]
            g_i_m3 = self.constraint_buffer_variables[identifier]["g_i-3"]

            buffer_value = self.constraint_buffer_variables[identifier]["buffer_value"]
            buffer_value_m1 = self.constraint_buffer_variables[identifier]["buffer_value-1"]

            if self.optimization_iteration > 3:
                delta_g_1 = g_i - g_i_m1
                delta_g_2 = g_i_m1 -g_i_m2
                delta_g_3 = g_i_m2 -g_i_m3

                if delta_g_1*delta_g_2 < 0 and delta_g_2*delta_g_3 < 0:
                    self.constraint_buffer_variables[identifier]["buffer_size_factor"] += abs(buffer_value - buffer_value_m1)

            if self.optimization_iteration > 1:
                delta_g = g_i - g_i_m1
                if delta_g >= 0.0 and g_i_m1 > 0:
                    self.constraint_buffer_variables[identifier]["central_buffer_value"] -= g_i_m1
                elif delta_g <= 0.0 and g_i_m1 < 0:
                    self.constraint_buffer_variables[identifier]["central_buffer_value"] += g_i_m1
                    self.constraint_buffer_variables[identifier]["central_buffer_value"] = \
                        max(self.constraint_buffer_variables[identifier]["central_buffer_value"], 0.0)

            self.constraint_buffer_variables[identifier]["g_i-3"] = g_i_m2
            self.constraint_buffer_variables[identifier]["g_i-2"] = g_i_m1
            self.constraint_buffer_variables[identifier]["g_i-1"] = g_i

        self.lower_bound.updateBufferZone(self.optimization_iteration)
        self.upper_bound.updateBufferZone(self.optimization_iteration)

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.d_norm
        additional_values_to_log["inf_norm_p"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_PROJECTION)
        additional_values_to_log["inf_norm_c"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_CORRECTION)
        additional_values_to_log["projection_norm"] = self.s_norm

        additional_values_to_log["t_lower_bound"] = self.lower_bound.value
        additional_values_to_log["t_upper_bound"] = self.upper_bound.value
        additional_values_to_log["t_lower_bound_buffer_value"] = self.lower_bound.buffer_variables["buffer_value"]
        additional_values_to_log["t_lower_bound_buffer_size"] = self.lower_bound.buffer_variables["buffer_size"]
        additional_values_to_log["t_lower_bound_buffer_size_factor"] = self.lower_bound.buffer_variables["buffer_size_factor"]
        additional_values_to_log["t_lower_bound_central_buffer_value"] = self.lower_bound.buffer_variables["central_buffer_value"]
        additional_values_to_log["t_lower_bound_lower_buffer_value"] = self.lower_bound.buffer_variables["lower_buffer_value"]
        additional_values_to_log["t_lower_bound_upper_buffer_value"] = self.lower_bound.buffer_variables["upper_buffer_value"]
        additional_values_to_log["t_upper_bound_buffer_value"] = self.upper_bound.buffer_variables["buffer_value"]
        additional_values_to_log["t_upper_bound_buffer_size"] = self.upper_bound.buffer_variables["buffer_size"]
        additional_values_to_log["t_upper_bound_buffer_size_factor"] = self.upper_bound.buffer_variables["buffer_size_factor"]
        additional_values_to_log["t_upper_bound_central_buffer_value"] = self.upper_bound.buffer_variables["central_buffer_value"]
        additional_values_to_log["t_upper_bound_lower_buffer_value"] = self.upper_bound.buffer_variables["lower_buffer_value"]
        additional_values_to_log["t_upper_bound_upper_buffer_value"] = self.upper_bound.buffer_variables["upper_buffer_value"]

        itr = 0
        for constraint in self.constraints:
            identifier = constraint["identifier"].GetString()
            additional_values_to_log["c"+str(itr+1)+"_buffer_value"] = self.constraint_buffer_variables[identifier]["buffer_value"]
            additional_values_to_log["c"+str(itr+1)+"_buffer_size"] = self.constraint_buffer_variables[identifier]["buffer_size"]
            additional_values_to_log["c"+str(itr+1)+"_buffer_size_factor"] = self.constraint_buffer_variables[identifier]["buffer_size_factor"]
            additional_values_to_log["c"+str(itr+1)+"_central_buffer_value"] = self.constraint_buffer_variables[identifier]["central_buffer_value"]
            additional_values_to_log["c"+str(itr+1)+"_lower_buffer_value"] = self.constraint_buffer_variables[identifier]["lower_buffer_value"]
            additional_values_to_log["c"+str(itr+1)+"_upper_buffer_value"] = self.constraint_buffer_variables[identifier]["upper_buffer_value"]
            itr += 1
        self.data_logger.LogCurrentValues(self.optimization_iteration, additional_values_to_log)
        self.data_logger.LogCurrentDesign(self.optimization_iteration)

    # --------------------------------------------------------------------------
    def __isAlgorithmConverged(self):

        if self.optimization_iteration > 1 :

            # Check if maximum iterations were reached
            if self.optimization_iteration == self.max_iterations:
                Kratos.Logger.Print("")
                Kratos.Logger.PrintInfo("ShapeOpt", "Maximal iterations of optimization problem reached!")
                return True

            # # Check for relative tolerance
            # relative_change_of_objective_value = self.data_logger.GetValues("rel_change_objective")[self.optimization_iteration]
            # if abs(relative_change_of_objective_value) < self.relative_tolerance:
            #     Kratos.Logger.Print("")
            #     Kratos.Logger.PrintInfo("ShapeOpt", "Optimization problem converged within a relative objective tolerance of ",self.relative_tolerance,"%.")
            #     return True

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges(self):
        # TODO: move to c++
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
            element_thicknesses[condition.Id] = condition.Properties.GetValue(Kratos.THICKNESS)

        self.__mapElementGradientToNode(element_thicknesses, KSO.THICKNESS)

# ==============================================================================

class LowerBound(object):

    def __init__(self, t_min, model_part):
        self.t_min = t_min
        self.value = 0.0
        self.model_part = model_part

        self.identifier = "lower thickness limit"

        self.nodal_values_variable = Kratos.KratosGlobals.GetVariable(f"LB_values")
        self.buffer_coeff_variable = Kratos.KratosGlobals.GetVariable(f"LB_buffer_coeff")
        self.gradient_variable = Kratos.KratosGlobals.GetVariable(f"DLBDT")
        self.mapped_gradient_variable = Kratos.KratosGlobals.GetVariable(f"DLBDT_MAPPED")

        self.model_part.AddNodalSolutionStepVariable(self.nodal_values_variable)
        self.model_part.AddNodalSolutionStepVariable(self.buffer_coeff_variable)
        self.model_part.AddNodalSolutionStepVariable(self.gradient_variable)
        self.model_part.AddNodalSolutionStepVariable(self.mapped_gradient_variable)

        self.buffer_variables = {
            "buffer_value": 0.0,
            "buffer_value-1": 0.0,
            "buffer_size": 1e-12,
            "buffer_size_factor": 2.0,
            "central_buffer_value": 0.0,
            "lower_buffer_value": - 1e-12,
            "upper_buffer_value": 1e-12,
            "g_i-1": 0.0,
            "g_i-2": 0.0,
            "g_i-3": 0.0,
            "max_constraint_change": 0.0
        }

        self.relaxation_coefficient = 0.0
        self.correction_coefficient = 0.0

    def CalculateValueAndGradient(self):

        g = - sys.float_info.max

        for node in self.model_part.Nodes:

            t_i = node.GetSolutionStepValue(KSO.THICKNESS)
            g_i = self.t_min - t_i

            # dg_i_dt = 0

            # max aggregation
            g = max(g, g_i)
            dg_i_dt = - 1

            # quadratic aggregation
            # if g_i > 0:
            #     g += g_i * g_i
            #     dg_i_dt = -1 * 2 * g_i

            node.SetSolutionStepValue(self.gradient_variable, 0, dg_i_dt)
            node.SetSolutionStepValue(self.nodal_values_variable, 0, g_i)

        self.value = g

    def ComputeBufferValue(self, optimization_iteration):

        g_i_m1 = self.buffer_variables["g_i-1"]
        buffer_size_factor = self.buffer_variables["buffer_size_factor"]

        self.buffer_variables["buffer_value-1"] = self.buffer_variables["buffer_value"]

        if optimization_iteration > 1:

            if abs(self.value - g_i_m1) > self.buffer_variables["max_constraint_change"]:
                self.buffer_variables["max_constraint_change"] = abs(self.value - g_i_m1)

            max_constraint_change = self.buffer_variables["max_constraint_change"]
            self.buffer_variables["buffer_size"] = max(buffer_size_factor * max_constraint_change, 1e-12)

        buffer_size = self.buffer_variables["buffer_size"]
        self.buffer_variables["lower_buffer_value"] = self.buffer_variables["central_buffer_value"] \
            - buffer_size
        self.buffer_variables["upper_buffer_value"] = self.buffer_variables["central_buffer_value"] \
            + buffer_size

        if self.isActive():
            lower_buffer_value = self.buffer_variables["lower_buffer_value"]
            self.buffer_variables["buffer_value"] = min( (self.value - lower_buffer_value) / buffer_size, 2.0 )
        else:
            self.buffer_variables["buffer_value"] = 0.0

    def getActive(self, optimization_utilities, design_surface):

        buffer_value = self.buffer_variables["buffer_value"]

        # normalize gradient
        g_a_variable = self.mapped_gradient_variable
        g_a_norm = optimization_utilities.ComputeMaxNormOfNodalVariable(design_surface, g_a_variable)
        g_a_variable_vector = Kratos.Vector()
        optimization_utilities.AssembleVector(design_surface, g_a_variable_vector, g_a_variable)
        if abs(g_a_norm) > 1e-10:
            g_a_variable_vector /= g_a_norm
        optimization_utilities.AssignVectorToVariable(design_surface, g_a_variable_vector, g_a_variable)

        active_relaxation_coefficient = min(buffer_value,1.0)
        self.relaxation_coefficient = active_relaxation_coefficient

        active_correction_coefficient = 0.0

        max_buffer = 2.0
        if buffer_value > 1.0:
            if buffer_value < max_buffer:
                active_correction_coefficient = 2*(buffer_value - 1)
            else:
                active_correction_coefficient = 2*(max_buffer-1)

        self.correction_coefficient = active_correction_coefficient

    def checkConstraintValue(self, buffer_coeff_update, optimization_utilities, design_surface):

        direction_has_changed = False

        if self.isActive():

            new_g_i = self._PredictResponseNodalValues(optimization_utilities, design_surface)
            new_g = np.max(new_g_i)

            Kratos.Logger.PrintInfo("Constraint ", self.identifier, "\n Linearized new value = ", new_g)
            if new_g > 0.0:
                if self.relaxation_coefficient < 1.0:
                    self.relaxation_coefficient = min(self.relaxation_coefficient + buffer_coeff_update, 1.0)
                    direction_has_changed = True
                elif self.correction_coefficient < 2.0:
                    self.correction_coefficient = min (self.correction_coefficient + buffer_coeff_update, 2.0)
                    direction_has_changed = True
            Kratos.Logger.PrintInfo("Constraint ", self.identifier, "\n W_R, W_C = ", self.relaxation_coefficient, self.correction_coefficient)

            self._ApplyNodeBasedCorrection(new_g_i, buffer_coeff_update, optimization_utilities, design_surface)

        return direction_has_changed

    def _PredictResponseNodalValues(self, optimization_utilities, design_surface):
        g_i = Kratos.Vector()
        optimization_utilities.AssembleVector(design_surface, g_i, self.nodal_values_variable)

        thickness_update = Kratos.Vector()
        optimization_utilities.AssembleVector(design_surface, thickness_update, KSO.THICKNESS_UPDATE)

        new_g_i = Kratos.Vector(len(g_i))
        for index, node in enumerate(design_surface.Nodes):
            if self.identifier == "lower thickness limit":
                new_g_i[index] = g_i[index] - thickness_update[index]
            elif self.identifier == "upper thickness limit":
                new_g_i[index] = g_i[index] + thickness_update[index]

        return new_g_i

    def _ApplyNodeBasedCorrection(self, new_g_i, buffer_coeff_update, optimization_utilities, design_surface):

        gradient = Kratos.Vector()
        optimization_utilities.AssembleVector(design_surface, gradient, self.mapped_gradient_variable)

        for index, node in enumerate(design_surface.Nodes):
            if new_g_i[index] < 0.0:
                gradient[index] /= 1 + buffer_coeff_update

        optimization_utilities.AssignVectorToVariable(design_surface, gradient, self.mapped_gradient_variable)

    def updateBufferZone(self, optimization_iteration):

        # update value history
        g_i = self.value
        g_i_m1 = self.buffer_variables["g_i-1"]
        g_i_m2 = self.buffer_variables["g_i-2"]

        self.buffer_variables["g_i-3"] = g_i_m2
        self.buffer_variables["g_i-2"] = g_i_m1
        self.buffer_variables["g_i-1"] = g_i

        # g_i = self.value
        # g_i_m1 = self.buffer_variables["g_i-1"]
        # g_i_m2 = self.buffer_variables["g_i-2"]
        # g_i_m3 = self.buffer_variables["g_i-3"]

        # buffer_value = self.buffer_variables["buffer_value"]
        # buffer_value_m1 = self.buffer_variables["buffer_value-1"]

        # if optimization_iteration > 3:
        #     delta_g_1 = g_i - g_i_m1
        #     delta_g_2 = g_i_m1 -g_i_m2
        #     delta_g_3 = g_i_m2 -g_i_m3

        #     if delta_g_1*delta_g_2 < 0 and delta_g_2*delta_g_3 < 0:
        #         self.buffer_variables["buffer_size_factor"] += abs(buffer_value - buffer_value_m1)

        # if optimization_iteration > 1:
        #     delta_g = g_i - g_i_m1
        #     if delta_g >= 0.0 and g_i_m1 > 0:
        #         self.buffer_variables["central_buffer_value"] -= g_i_m1
        #     elif delta_g <= 0.0 and g_i_m1 < 0:
        #         self.buffer_variables["central_buffer_value"] += g_i_m1
        #         self.buffer_variables["central_buffer_value"] = \
        #             max(self.buffer_variables["central_buffer_value"], 0.0)

        # self.buffer_variables["g_i-3"] = g_i_m2
        # self.buffer_variables["g_i-2"] = g_i_m1
        # self.buffer_variables["g_i-1"] = g_i

    def ComputeBufferCoefficients(self, optimization_utilities, design_surface):

        g_i = Kratos.Vector()
        optimization_utilities.AssembleVector(design_surface, g_i, self.nodal_values_variable)

        buffer_coeff = Kratos.Vector()
        optimization_utilities.AssembleVector(design_surface, buffer_coeff, self.buffer_coeff_variable)

        for index, node in enumerate(design_surface.Nodes):
            if g_i[index] > self.buffer_variables["lower_buffer_value"]:
                buffer_coeff[index] = (g_i[index] - self.buffer_variables["lower_buffer_value"]) / (
                    self.buffer_variables["upper_buffer_value"] - self.buffer_variables["lower_buffer_value"]) * 2.0
            else:
                buffer_coeff[index] = 0.0

        optimization_utilities.AssignVectorToVariable(design_surface, buffer_coeff, self.buffer_coeff_variable)

        return buffer_coeff

    def AggregateNodalConstraints(self, optimization_utilities, design_surface):
        buffer_coefficients = self.ComputeBufferCoefficients(optimization_utilities, design_surface)

        gradient = Kratos.Vector()
        optimization_utilities.AssembleVector(design_surface, gradient, self.gradient_variable)
        for index, node in enumerate(design_surface.Nodes):
            gradient[index] *= buffer_coefficients[index]

        optimization_utilities.AssignVectorToVariable(design_surface, gradient, self.gradient_variable)

    def isActive(self):

        if self.value >= self.buffer_variables["lower_buffer_value"]:
            return True
        else:
            return False


class UpperBound(LowerBound):

    def __init__(self, t_max, model_part):
        self.t_max = t_max
        self.value = 0.0
        self.model_part = model_part

        self.nodal_values_variable = Kratos.KratosGlobals.GetVariable(f"UB_values")
        self.buffer_coeff_variable = Kratos.KratosGlobals.GetVariable(f"UB_buffer_coeff")
        self.gradient_variable = Kratos.KratosGlobals.GetVariable(f"DUBDT")
        self.mapped_gradient_variable = Kratos.KratosGlobals.GetVariable(f"DUBDT_MAPPED")

        self.model_part.AddNodalSolutionStepVariable(self.nodal_values_variable)
        self.model_part.AddNodalSolutionStepVariable(self.buffer_coeff_variable)
        self.model_part.AddNodalSolutionStepVariable(self.gradient_variable)
        self.model_part.AddNodalSolutionStepVariable(self.mapped_gradient_variable)

        self.buffer_variables = {
            "buffer_value": 0.0,
            "buffer_value-1": 0.0,
            "buffer_size": 1e-12,
            "buffer_size_factor": 2.0,
            "central_buffer_value": 0.0,
            "lower_buffer_value": - 1e-12,
            "upper_buffer_value": 1e-12,
            "g_i-1": 0.0,
            "g_i-2": 0.0,
            "g_i-3": 0.0,
            "max_constraint_change": 0.0
        }

        self.identifier = "upper thickness limit"

    def CalculateValueAndGradient(self):

        g = - sys.float_info.max

        for node in self.model_part.Nodes:

            t_i = node.GetSolutionStepValue(KSO.THICKNESS)
            g_i = t_i - self.t_max
            # dg_i_dt = 0

            # max aggregation
            g = max(g, g_i)
            dg_i_dt = 1

            # quadratic aggregation
            # if g_i > 0:
            #     g += g_i * g_i
            #     dg_i_dt = 2 * g_i

            node.SetSolutionStepValue(self.gradient_variable, 0, dg_i_dt)
            node.SetSolutionStepValue(self.nodal_values_variable, 0, g_i)

        self.value = g