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
# from KratosMultiphysics.ShapeOptimizationApplication.algorithms import algorithm_gradient_projection
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
from KratosMultiphysics.ShapeOptimizationApplication.loggers import data_logger_factory
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer
from KratosMultiphysics.ShapeOptimizationApplication import model_part_controller_factory
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_variable_utilities import WriteDictionaryDataOnNodalVariable

import math
import numpy as np

# ==============================================================================
class AlgorithmFreeThicknessOptimizationv2RGP(OptimizationAlgorithm):
    """
        Algorithm for free thickness optimization using a filtering operation which is
        defined on the initial geometry (thus not changing over the course of the optimization)
        and filtering the total thickness update t = T + dt = T + A * dt_control.
        NOT READY YET: Shape optimization is conducted seperately and not simultaneously.
    """
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
        self.shape_opt = False

        self.algorithm_settings =  optimization_settings["optimization_algorithm"]
        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

        self.optimization_settings = optimization_settings
        self.mapper_settings = optimization_settings["design_variables"]["filter"]

        if optimization_settings["design_variables"].Has("projection"):
            self.projection = optimization_settings["design_variables"]["projection"]
        else:
            self.projection = False

        if self.projection:
            self.thickness_targets = optimization_settings["design_variables"]["projection_settings"]["available_thicknesses"].GetVector()
            self.beta = optimization_settings["design_variables"]["projection_settings"]["initial_beta"].GetDouble()
            self.q = optimization_settings["design_variables"]["projection_settings"]["increase_beta_factor"].GetDouble()
            self.deactivated_design_variables_indices = []

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.design_surface = None
        self.mapper = None
        self.data_logger = None
        self.optimization_utilities = None
        self.variable_utils = None

        self.objectives = optimization_settings["objectives"]
        self.constraints = optimization_settings["constraints"]

        if self.shape_opt:
            self.shape_constraint_gradient_variables = {}

        self.constraint_gradient_variables = {}
        self.constraint_buffer_variables = {}
        for itr, constraint in enumerate(self.constraints):
            self.constraint_gradient_variables.update({
                constraint["identifier"].GetString() : {
                    "gradient": KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT"),
                    "mapped_gradient": KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT_MAPPED"),
                    "projected_gradient": KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT_PROJECTED")
                }
            })
            self.constraint_buffer_variables.update({
                constraint["identifier"].GetString() : {
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
            if self.shape_opt:
                self.shape_constraint_gradient_variables.update({
                    constraint["identifier"].GetString() : {
                        "gradient": KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DX"),
                        "mapped_gradient": KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DX_MAPPED")
                    }
                })


        # create initial model and model_controller
        initial_model = KM.Model()
        initial_model_settings = optimization_settings["model_settings"].Clone()
        initial_model_settings["model_part_name"].SetString(optimization_settings["model_settings"]["model_part_name"].GetString() + "_initial")
        self.initial_model_part_controller = model_part_controller_factory.CreateController(initial_model_settings, initial_model)

        initial_model_part = self.initial_model_part_controller.GetOptimizationModelPart()

        KM.Logger.PrintInfo("ShapeOpt", "Model Part: ", self.model_part_controller.GetOptimizationModelPart().Name)
        KM.Logger.PrintInfo("ShapeOpt", "Initial Model Part: ", self.initial_model_part_controller.GetOptimizationModelPart().Name)

        nodal_variable = KM.KratosGlobals.GetVariable("DF1DT")
        initial_model_part.AddNodalSolutionStepVariable(nodal_variable)
        nodal_variable = KM.KratosGlobals.GetVariable("DF1DT_MAPPED")
        initial_model_part.AddNodalSolutionStepVariable(nodal_variable)

        for itr, constraint in enumerate(self.constraints):
            nodal_variable = KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT")
            initial_model_part.AddNodalSolutionStepVariable(nodal_variable)
            nodal_variable = KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT_MAPPED")
            initial_model_part.AddNodalSolutionStepVariable(nodal_variable)

        initial_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CHANGE_CONTROL_PROJECTED)
        initial_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CHANGE_CONTROL)
        initial_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CHANGE)

        self.max_correction_share = self.algorithm_settings["max_correction_share"].GetDouble()

        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.line_search_type = self.algorithm_settings["line_search"]["line_search_type"].GetString()
        if self.shape_opt:
            self.shape_step_size = 0.5
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_SEARCH_DIRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CORRECTION)

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Gradient projection algorithm only supports one objective function!")
        if self.constraints.size() == 0:
            raise RuntimeError("Gradient projection algorithm requires definition of at least one constraint!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.initial_model_part_controller.Initialize()
        self.model_part_controller.Initialize()

        self.analyzer.InitializeBeforeOptimizationLoop()

        self.initial_design_surface = self.initial_model_part_controller.GetDesignSurface()
        self.design_surface = self.model_part_controller.GetDesignSurface()

        self.initial_mapper = mapper_factory.CreateMapper(self.initial_design_surface, self.initial_design_surface, self.mapper_settings)
        self.initial_mapper.Initialize()
        self.initial_model_part_controller.InitializeDamping()

        self.mapper = mapper_factory.CreateMapper(self.design_surface, self.design_surface, self.mapper_settings)
        self.mapper.Initialize()
        self.model_part_controller.InitializeDamping()

        optimization_settings = self.optimization_settings.Clone()
        self.data_logger = data_logger_factory.CreateDataLogger(self.model_part_controller, self.communicator, self.optimization_settings)
        self.data_logger.InitializeDataLogging()

        optimization_settings["output"]["design_output_mode"].SetString("WriteOptimizationModelPart")
        self.initial_data_logger = data_logger_factory.CreateDataLogger(self.initial_model_part_controller, self.communicator, optimization_settings)
        self.initial_data_logger.InitializeDataLogging()

        self.optimization_utilities = KSO.OptimizationUtilities
        self.variable_utils = KM.VariableUtils()

        self.__InitializeThicknessField()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        for self.optimization_iteration in range(1,self.max_iterations):
            KM.Logger.Print("")
            KM.Logger.Print("===============================================================================")
            KM.Logger.PrintInfo("TicknessOpt", timer.GetTimeStamp(), ": Starting optimization iteration ", self.optimization_iteration)
            KM.Logger.Print("===============================================================================\n")

            timer.StartNewLap()

            self.__initializeNewThickness()

            if self.shape_opt:
                self.__initializeNewShape()

            self.__analyzeThickness()

            self.__computeBufferValue()

            if self.shape_opt:
                self.__analyzeShape()

            self.__computeThicknessUpdate()

            if self.shape_opt:
                self.__computeShapeUpdate()

            self.__logCurrentOptimizationStep()

            self.__updateBufferZone()

            KM.Logger.Print("")
            KM.Logger.PrintInfo("TicknessOpt", "Time needed for current optimization step = ", timer.GetLapTime(), "s")
            KM.Logger.PrintInfo("TicknessOpt", "Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

            if self.__isAlgorithmConverged():
                break
            else:
                self.__updateBeta()
                if self.shape_opt:
                    self.__determineAbsoluteShapeChanges()

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.data_logger.FinalizeDataLogging()
        self.analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeNewThickness(self):
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)
        self.model_part_controller.UpdateThicknessAccordingInitialAndInputVariable(KSO.THICKNESS_CHANGE)

        # update solution step variables
        for node in self.optimization_model_part.Nodes:
            initial_thickness = node.GetSolutionStepValue(KSO.THICKNESS_INITIAL)
            thickness_change = node.GetSolutionStepValue(KSO.THICKNESS_CHANGE)

            new_thickness = initial_thickness + thickness_change
            node.SetSolutionStepValue(KSO.THICKNESS, 0, new_thickness)

        self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyzeThickness(self):
        self.communicator.initializeCommunication()
        self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
        self.communicator.requestThicknessGradientOf(self.objectives[0]["identifier"].GetString())
        if self.shape_opt:
            self.communicator.requestGradientOf(self.objectives[0]["identifier"].GetString())

        for constraint in self.constraints:
            con_id =  constraint["identifier"].GetString()
            self.communicator.requestValueOf(con_id)
            self.communicator.requestThicknessGradientOf(con_id)
            if self.shape_opt:
                self.communicator.requestGradientOf(con_id)

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.optimization_iteration, self.communicator)

        # project and damp objective gradients
        objGradientDict = self.communicator.getStandardizedThicknessGradient(self.objectives[0]["identifier"].GetString())
        objElementGradientDict = dict()
        self.__mapPropertyGradientToElement(objGradientDict, objElementGradientDict)
        self.__mapElementGradientToNode(objElementGradientDict, KSO.DF1DT)

        self.variable_utils.CopyModelPartNodalVar(KSO.DF1DT, self.model_part_controller.GetOptimizationModelPart(),
                                                  self.initial_model_part_controller.GetOptimizationModelPart(), 0)

        # self.initial_model_part_controller.DampNodalSensitivityVariableIfSpecified(KSO.DF1DT)
        # self.model_part_controller.DampNodalSensitivityVariableIfSpecified(KSO.DF1DT)

        # project and damp constraint gradients
        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            conGradientDict = self.communicator.getStandardizedThicknessGradient(con_id)
            gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
            conElementGradientDict = dict()
            self.__mapPropertyGradientToElement(conGradientDict, conElementGradientDict)
            self.__mapElementGradientToNode(conElementGradientDict, gradient_variable)

            self.variable_utils.CopyModelPartNodalVar(gradient_variable, self.model_part_controller.GetOptimizationModelPart(),
                                                      self.initial_model_part_controller.GetOptimizationModelPart(), 0)

            # self.initial_model_part_controller.DampNodalSensitivityVariableIfSpecified(gradient_variable)
            # self.model_part_controller.DampNodalSensitivityVariableIfSpecified(gradient_variable)

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

    # --------------------------------------------------------------------------
    def __computeThicknessUpdate(self):

        ### 1. Filter gradients
        self.initial_mapper.Update()
        self.initial_mapper.InverseMap(KSO.DF1DT, KSO.DF1DT_MAPPED)
        self.variable_utils.CopyModelPartNodalVar(KSO.DF1DT_MAPPED, self.initial_model_part_controller.GetOptimizationModelPart(),
                                                  self.model_part_controller.GetOptimizationModelPart(), 0)

        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
            mapped_gradient_variable = self.constraint_gradient_variables[con_id]["mapped_gradient"]
            self.initial_mapper.InverseMap(gradient_variable, mapped_gradient_variable)
            self.variable_utils.CopyModelPartNodalVar(mapped_gradient_variable, self.initial_model_part_controller.GetOptimizationModelPart(),
                                                      self.model_part_controller.GetOptimizationModelPart(), 0)


        ### 2. Project gradients to control thickness
        self.__ProjectGradients(KSO.DF1DT_MAPPED, KSO.DF1DT_PROJECTED)

        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            mapped_gradient_variable = self.constraint_gradient_variables[con_id]["mapped_gradient"]
            projected_gradient_variable = self.constraint_gradient_variables[con_id]["projected_gradient"]
            self.__ProjectGradients(mapped_gradient_variable, projected_gradient_variable)

        ### 3. Compute control thickness update via gradient projection
        self.__computeControlThicknessUpdate()

        ### 4. Project control thickness to obtain filtered thickness
        self.__ProjectThickness()

        ### 5. Filter thickness to obtain physical thickness
        if self.projection:
            self.variable_utils.CopyModelPartNodalVar(KSO.THICKNESS_CHANGE_CONTROL_PROJECTED, self.model_part_controller.GetOptimizationModelPart(),
                                                      self.initial_model_part_controller.GetOptimizationModelPart(), 0)
            self.initial_mapper.Map(KSO.THICKNESS_CHANGE_CONTROL_PROJECTED, KSO.THICKNESS_CHANGE)
            # self.mapper.Map(KSO.THICKNESS_CHANGE_CONTROL_PROJECTED, KSO.THICKNESS_CHANGE)

        else:
            self.variable_utils.CopyModelPartNodalVar(KSO.THICKNESS_CHANGE_CONTROL, self.model_part_controller.GetOptimizationModelPart(),
                                                      self.initial_model_part_controller.GetOptimizationModelPart(), 0)
            self.initial_mapper.Map(KSO.THICKNESS_CHANGE_CONTROL, KSO.THICKNESS_CHANGE)
            # self.mapper.Map(KSO.THICKNESS_CHANGE_CONTROL, KSO.THICKNESS_CHANGE)

        # self.initial_model_part_controller.DampNodalUpdateVariableIfSpecified(KSO.THICKNESS_CHANGE)
        # self.model_part_controller.DampNodalUpdateVariableIfSpecified(KSO.THICKNESS_CHANGE)

        self.variable_utils.CopyModelPartNodalVar(KSO.THICKNESS_CHANGE, self.initial_model_part_controller.GetOptimizationModelPart(),
                                                  self.model_part_controller.GetOptimizationModelPart(), 0)

        self.__saveLineSearchData()

    # --------------------------------------------------------------------------
    def __LineSearch(self):
        KM.Logger.PrintInfo("Line Search ...")
        if self.line_search_type == "manual_stepping":
            self.__manualStep()
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
            step = KM.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, step, KSO.THICKNESS_CONTROL_UPDATE)
            step *= 1.0 / step_norm
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, step, KSO.THICKNESS_SEARCH_DIRECTION)
            step *= self.step_size
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, step, KSO.THICKNESS_CONTROL_UPDATE)

    def __BBStep(self):
        self.s_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_SEARCH_DIRECTION)
        if abs(self.s_norm) > 1e-10:
            s = KM.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, s, KSO.THICKNESS_SEARCH_DIRECTION)
            s /= self.s_norm
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_SEARCH_DIRECTION)
            y = self.prev_s - s


            print(f"self.deactivated_design_variables_indices: {self.deactivated_design_variables_indices}")
            if self.projection:
                for deact_design_var_index in self.deactivated_design_variables_indices:
                    self.d[deact_design_var_index] = 0.0
                    y[deact_design_var_index] = 0.0

            if np.dot(y, y) < 1e-9:
                step = self.max_step_size
            else:
                step = abs(np.dot(y, self.d) / np.dot(y, y))
            if step > self.max_step_size:
                step = self.max_step_size
            s = s * step
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_CONTROL_UPDATE)

    def __saveLineSearchData(self):
        self.prev_s = KM.Vector()
        self.d = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, self.d, KSO.THICKNESS_CONTROL_UPDATE)
        self.optimization_utilities.AssembleVector(self.design_surface, self.prev_s, KSO.THICKNESS_SEARCH_DIRECTION)

    # --------------------------------------------------------------------------
    def __computeControlThicknessUpdate(self):
        """adapted from https://msulaiman.org/onewebmedia/GradProj_2.pdf"""
        g_a, g_a_variables, self.relaxation_coefficients, self.correction_coefficients = self.__getActiveConstraints()

        prev_s = KM.Vector(len(self.design_surface.Nodes))

        max_inner_iter = 10
        for inner_iter in range(max_inner_iter):

            KM.Logger.PrintInfo("ThicknessOpt", f"Inner Iteration: {inner_iter+1}")
            KM.Logger.PrintInfo("ThicknessOpt", "Assemble vector of objective gradient.")
            nabla_f = KM.Vector()
            if self.projection:
                self.optimization_utilities.AssembleVector(self.design_surface, nabla_f, KSO.DF1DT_PROJECTED)
            else:
                self.optimization_utilities.AssembleVector(self.design_surface, nabla_f, KSO.DF1DT_MAPPED)

            # normalize objective gradient
            if abs(nabla_f.norm_inf()) > 1e-10:
                nabla_f *= 1.0/nabla_f.norm_inf()

            s = KM.Vector()
            p = KM.Vector()
            print(f"Norm of objective sens: {nabla_f.norm_inf()}")

            # normalize constraint gradients
            for g_a_variable in g_a_variables:
                constraint_gradient = KM.Vector()
                self.optimization_utilities.AssembleVector(self.design_surface, constraint_gradient, g_a_variable)
                if abs(constraint_gradient.norm_inf()) > 1e-10:
                    constraint_gradient *= 1.0/constraint_gradient.norm_inf()
                self.optimization_utilities.AssignVectorToVariable(self.design_surface, constraint_gradient, g_a_variable)
                print(f"Norm of constraint sens: {constraint_gradient.norm_inf()}")

            delta_t_control = KM.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, delta_t_control, KSO.THICKNESS_CHANGE_CONTROL)
            delta_t_control_update = KM.Vector()

            if len(g_a) == 0:
                KM.Logger.PrintInfo("ThicknessOpt", "No constraints active, use negative objective gradient as search direction.")
                s = nabla_f * (-1.0)

                if self.projection:
                    s = self.__ProjectSearchDirectionAndGradients(s, g_a_variables)

                if not self.projection or (s - prev_s).norm_inf() == 0.0 or inner_iter == max_inner_iter - 1:
                    self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_PROJECTION)
                    self.optimization_utilities.AssignVectorToVariable(self.design_surface, [0.0]*len(s), KSO.THICKNESS_CORRECTION)
                    self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_SEARCH_DIRECTION)
                    self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_CONTROL_UPDATE)
                    self.__LineSearch()
                    self.optimization_utilities.AssembleVector(self.design_surface, delta_t_control_update, KSO.THICKNESS_CONTROL_UPDATE)
                    self.optimization_utilities.AssignVectorToVariable(self.design_surface, delta_t_control+delta_t_control_update, KSO.THICKNESS_CHANGE_CONTROL)
                    return
                else:
                    prev_s = s
                    continue

            omega_r = KM.Matrix()
            self.optimization_utilities.AssembleBufferMatrix(omega_r, self.relaxation_coefficients)
            omega_c = KM.Vector(self.correction_coefficients)

            KM.Logger.PrintInfo("ThicknessOpt", "Assemble matrix of constraint gradient.")
            N = KM.Matrix()
            self.optimization_utilities.AssembleMatrixScalarVariables(self.design_surface, N, g_a_variables)

            settings = KM.Parameters('{ "solver_type" : "LinearSolversApplication.dense_col_piv_householder_qr" }')
            solver = dense_linear_solver_factory.ConstructSolver(settings)

            KM.Logger.PrintInfo("ThicknessOpt", "Calculate projected search direction and correction.")
            c = KM.Vector()
            self.optimization_utilities.CalculateRelaxedProjectedSearchDirectionAndCorrection(
                nabla_f,
                N,
                omega_r,
                omega_c,
                solver,
                p,
                c)

            s = p+c
            print(f"Norm of projection: {p.norm_inf()}")
            print(f"Norm of correction: {c.norm_inf()}")
            print(f"Norm of search direction: {s.norm_inf()}")

            if self.projection:
                s = self.__ProjectSearchDirectionAndGradients(s, g_a_variables)

            if not self.projection or (s - prev_s).norm_inf() == 0.0 or inner_iter == max_inner_iter - 1:
                self.optimization_utilities.AssignVectorToVariable(self.design_surface, p, KSO.THICKNESS_PROJECTION)
                self.optimization_utilities.AssignVectorToVariable(self.design_surface, c, KSO.THICKNESS_CORRECTION)
                self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_SEARCH_DIRECTION)
                self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.THICKNESS_CONTROL_UPDATE)
                self.__LineSearch()
                self.optimization_utilities.AssembleVector(self.design_surface, delta_t_control_update, KSO.THICKNESS_CONTROL_UPDATE)
                self.optimization_utilities.AssignVectorToVariable(self.design_surface, delta_t_control+delta_t_control_update, KSO.THICKNESS_CHANGE_CONTROL)
                return
            else:
                prev_s = s
                continue

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
                if self.projection:
                    g_a_variable = self.constraint_gradient_variables[identifier]["projected_gradient"]
                else:
                    g_a_variable = self.constraint_gradient_variables[identifier]["mapped_gradient"]

                g_a_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, g_a_variable)
                g_a_variable_vector = KM.Vector()
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

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):

        for node in self.optimization_model_part.Nodes:
            initial_thickness = node.GetSolutionStepValue(KSO.THICKNESS_INITIAL)
            thickness_change = node.GetSolutionStepValue(KSO.THICKNESS_CHANGE)
            old_thickness = node.GetSolutionStepValue(KSO.THICKNESS)

            new_thickness = initial_thickness + thickness_change
            thickness_update = new_thickness - old_thickness

            node.SetSolutionStepValue(KSO.THICKNESS_UPDATE, 0, thickness_update)

        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_UPDATE)
        additional_values_to_log["inf_norm_p"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_PROJECTION)
        additional_values_to_log["inf_norm_c"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_CORRECTION)
        additional_values_to_log["projection_norm"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_SEARCH_DIRECTION)

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

        # self.data_logger.LogSensitivityHeatmap(self.optimization_iteration, self.mapper)
        self.data_logger.LogCurrentValues(self.optimization_iteration, additional_values_to_log)
        self.data_logger.LogCurrentDesign(self.optimization_iteration)
        # self.initial_data_logger.LogCurrentDesign(self.optimization_iteration)

    # --------------------------------------------------------------------------
    def __isAlgorithmConverged(self):

        if self.optimization_iteration > 1 :

            # Check if maximum iterations were reached
            if self.optimization_iteration == self.max_iterations:
                KM.Logger.Print("")
                KM.Logger.PrintInfo("ThicknessOpt", "Maximal iterations of optimization problem reached!")
                return True

            # # # Check for relative tolerance
            # relative_change_of_objective_value = self.data_logger.GetValues("rel_change_objective")[self.optimization_iteration]
            # if abs(relative_change_of_objective_value) < self.relative_tolerance:
            #     KM.Logger.Print("")
            #     KM.Logger.PrintInfo("ThicknessOpt", "Optimization problem converged within a relative objective tolerance of ",self.relative_tolerance,"%.")
            #     return True

    # --------------------------------------------------------------------------
    def __updateBeta(self):

        if self.projection:
            self.beta *= self.q

    # --------------------------------------------------------------------------
    def __determineAbsoluteShapeChanges(self):
        self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.CONTROL_POINT_UPDATE, KSO.CONTROL_POINT_CHANGE)
        self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.SHAPE_UPDATE, KSO.SHAPE_CHANGE)

    # --------------------------------------------------------------------------
    def __mapPropertyGradientToElement(self, property_gradient_dict, element_gradient_dict):

        for condition in self.optimization_model_part.Conditions:
            element_gradient_dict[condition.Id] = property_gradient_dict[condition.Properties.Id]

    # --------------------------------------------------------------------------
    def __mapElementGradientToNode(self, element_gradient_dict, gradient_variable):

        # reset variables
        for node in self.optimization_model_part.Nodes:
            node.SetSolutionStepValue(gradient_variable, 0, 0.0)

        for condition in self.optimization_model_part.Conditions:
            condition.SetValue(gradient_variable, element_gradient_dict[condition.Id])

        total_node_areas = dict()
        for condition in self.optimization_model_part.Conditions:
            df_dt = condition.GetValue(gradient_variable)
            for node in condition.GetNodes():
                if node.Id in total_node_areas:
                    total_node_areas[node.Id] += condition.GetGeometry().Area()
                else:
                    total_node_areas[node.Id] = condition.GetGeometry().Area()
                df_dt_node = node.GetSolutionStepValue(gradient_variable)
                df_dt_node += df_dt * condition.GetGeometry().Area()
                node.SetSolutionStepValue(gradient_variable, 0, df_dt_node)

        for node in self.optimization_model_part.Nodes:
            total_node_area = total_node_areas[node.Id]
            df_dt_node = node.GetSolutionStepValue(gradient_variable)
            df_dt_node /= total_node_area
            node.SetSolutionStepValue(gradient_variable, 0, df_dt_node)

    # --------------------------------------------------------------------------
    def __InitializeThicknessField(self):

        element_thicknesses = dict()
        for condition in self.optimization_model_part.Conditions:
            element_thicknesses[condition.Id] = condition.Properties.GetValue(KM.THICKNESS)
            condition.Properties.SetValue(KSO.THICKNESS_INITIAL, element_thicknesses[condition.Id])

        self.__mapElementGradientToNode(element_thicknesses, KSO.THICKNESS)
        self.__mapElementGradientToNode(element_thicknesses, KSO.THICKNESS_INITIAL)

    # --------------------------------------------------------------------------
    def __ProjectGradients(self, mapped_gradient_variable, projected_gradient_variable):

        if not self.projection:
            return

        for node in self.optimization_model_part.Nodes:
            delta_t_control = node.GetSolutionStepValue(KSO.THICKNESS_CHANGE_CONTROL)
            mapped_gradient = node.GetSolutionStepValue(mapped_gradient_variable)
            # print(f"t: {t}")

            delta_t_m = self.___GetInterval(node, delta_t_control)
            # print(f"t_m: {t_m}")

            w = (delta_t_control - delta_t_m[0]) / (delta_t_m[1] - delta_t_m[0])
            tanh_numerator = 2 * np.tanh(self.beta * 0.5)
            dw_dt_p = 1 / (delta_t_m[1] - delta_t_m[0])
            dt_p_dw = ((delta_t_m[1] - delta_t_m[0]) / tanh_numerator) \
                * (1 - np.tanh(self.beta * (w - 0.5))) * self.beta
            projected_gradient = mapped_gradient * dw_dt_p * dt_p_dw

            node.SetSolutionStepValue(projected_gradient_variable, 0, projected_gradient)

    def __ProjectSearchDirectionAndGradients(self, search_direction, g_a_variables):

        if not self.projection:
            return search_direction

        delta_t_control = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, delta_t_control, KSO.THICKNESS_CHANGE_CONTROL)

        nabla_f = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f, KSO.DF1DT_PROJECTED)

        g_a = {}
        for g_a_variable in g_a_variables:
            constraint_gradient = KM.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, constraint_gradient, g_a_variable)
            g_a[g_a_variable] = constraint_gradient

        t_init = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, t_init, KSO.THICKNESS_INITIAL)

        for i in range(len(self.design_surface.Nodes)):
            t = delta_t_control[i]
            delta_t_targets = np.zeros(len(self.thickness_targets))
            for j in range(len(self.thickness_targets)):
                delta_t_targets[j] = self.thickness_targets[j] - t_init[i]

            if t <= delta_t_targets[0] and search_direction[i] < 0.0:
                self.deactivated_design_variables_indices.append(i)
                search_direction[i] = 0.0
                nabla_f[i] = 0.0
                for constraint_gradient in g_a.values():
                    constraint_gradient[i] = 0.0
            elif t >= delta_t_targets[-1] and search_direction[i] > 0.0:
                self.deactivated_design_variables_indices.append(i)
                search_direction[i] = 0.0
                nabla_f[i] = 0.0
                for constraint_gradient in g_a.values():
                    constraint_gradient[i] = 0.0

        print(f"ProjectSearchDirection: self.deactivated_design_variables_indices: {self.deactivated_design_variables_indices}")

        self.optimization_utilities.AssignVectorToVariable(self.design_surface, nabla_f, KSO.DF1DT_PROJECTED)
        for g_a_variable, constraint_gradient in g_a.items():
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, constraint_gradient, g_a_variable)

        return search_direction

    def __ProjectThickness(self):

        if not self.projection:
            return

        for node in self.optimization_model_part.Nodes:
            delta_t_control = node.GetSolutionStepValue(KSO.THICKNESS_CHANGE_CONTROL)

            bound = self.__CheckBounds(node, delta_t_control)
            if bound:
                delta_t_projected = bound
            else:
                delta_t_m = self.___GetInterval(node, delta_t_control)

                w = (delta_t_control - delta_t_m[0]) / (delta_t_m[1] - delta_t_m[0])
                tanh_denominator = np.tanh(self.beta * 0.5) + np.tanh(self.beta * (w - 0.5))
                tanh_numerator = 2 * np.tanh(self.beta * 0.5)
                tanh_term = tanh_denominator / tanh_numerator
                delta_t_projected = delta_t_m[0] + (delta_t_m[1] - delta_t_m[0]) * tanh_term

            node.SetSolutionStepValue(KSO.THICKNESS_CHANGE_CONTROL_PROJECTED, 0, delta_t_projected)


    def ___GetInterval(self, node, delta_t_projected):

        t_init = node.GetSolutionStepValue(KSO.THICKNESS_INITIAL)

        delta_t_targets = KM.Vector(len(self.thickness_targets))

        for i in range(len(self.thickness_targets)):
            delta_t_targets[i] = self.thickness_targets[i] - t_init

        if delta_t_projected <= delta_t_targets[0]:
            return np.array([delta_t_targets[0], delta_t_targets[1]])
        if delta_t_projected >= delta_t_targets[len(delta_t_targets)-1]:
            return np.array([delta_t_targets[len(delta_t_targets)-2], delta_t_targets[len(delta_t_targets)-1]])

        for i in range(len(delta_t_targets)-1):
            if delta_t_projected >= delta_t_targets[i] and delta_t_projected <= delta_t_targets[i+1]:
                return np.array([delta_t_targets[i], delta_t_targets[i+1]])

    def __CheckBounds(self, node, delta_t_projected):

        t_init = node.GetSolutionStepValue(KSO.THICKNESS_INITIAL)

        delta_t_targets = KM.Vector(len(self.thickness_targets))

        for i in range(len(self.thickness_targets)):
            delta_t_targets[i] = self.thickness_targets[i] - t_init

        if delta_t_projected <= delta_t_targets[0]:
            return delta_t_targets[0]
        elif delta_t_projected >= delta_t_targets[len(delta_t_targets)-1]:
            return delta_t_targets[len(delta_t_targets)-1]

        return None


## Kopie aus gradient projection für shape optimization
    # --------------------------------------------------------------------------
    def __initializeNewShape(self):
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)
        self.model_part_controller.UpdateMeshAccordingInputVariable(KSO.SHAPE_UPDATE)
        self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyzeShape(self):
        # self.communicator.initializeCommunication()
        # self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
        # self.communicator.requestGradientOf(self.objectives[0]["identifier"].GetString())

        # for constraint in self.constraints:
        #     con_id =  constraint["identifier"].GetString()
        #     self.communicator.requestValueOf(con_id)
        #     self.communicator.requestGradientOf(con_id)

        # self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.optimization_iteration, self.communicator)

        # compute normals only if required
        surface_normals_required = self.objectives[0]["project_gradient_on_surface_normals"].GetBool()
        for constraint in self.constraints:
            if constraint["project_gradient_on_surface_normals"].GetBool():
                surface_normals_required = True

        if surface_normals_required:
            self.model_part_controller.ComputeUnitSurfaceNormals()

        # project and damp objective gradients
        objGradientDict = self.communicator.getStandardizedGradient(self.objectives[0]["identifier"].GetString())
        # print(f"shape objGradientDict: {objGradientDict}")
        WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, KSO.DF1DX)
        nabla_f_raw = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f_raw, KSO.DF1DX)
        print(f"Max objective gradient raw : nabla_f_raw : {nabla_f_raw.norm_inf()}")

        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(KSO.DF1DX)

        self.model_part_controller.DampNodalSensitivityVariableIfSpecified(KSO.DF1DX)

        nabla_f = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f, KSO.DF1DX)
        print(f"Max objective gradient : nabla_f : {nabla_f.norm_inf()}")

        # project and damp constraint gradients
        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            conGradientDict = self.communicator.getStandardizedGradient(con_id)
            gradient_variable = self.shape_constraint_gradient_variables[con_id]["gradient"]
            WriteDictionaryDataOnNodalVariable(conGradientDict, self.optimization_model_part, gradient_variable)

            if constraint["project_gradient_on_surface_normals"].GetBool():
                self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(gradient_variable)

            self.model_part_controller.DampNodalSensitivityVariableIfSpecified(gradient_variable)

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):
        self.mapper.Update()
        self.mapper.InverseMap(KSO.DF1DX, KSO.DF1DX_MAPPED)

        nabla_f_mapped = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f_mapped, KSO.DF1DX_MAPPED)
        print(f"Max objective gradient mapped : nabla_f_mapped : {nabla_f_mapped.norm_inf()}")


        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            gradient_variable = self.shape_constraint_gradient_variables[con_id]["gradient"]
            mapped_gradient_variable = self.shape_constraint_gradient_variables[con_id]["mapped_gradient"]
            self.mapper.InverseMap(gradient_variable, mapped_gradient_variable)

        self.__computeControlPointUpdate()

        self.mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
        self.model_part_controller.DampNodalUpdateVariableIfSpecified(KSO.SHAPE_UPDATE)

    # --------------------------------------------------------------------------
    def __computeControlPointUpdate(self):
        """adapted from https://msulaiman.org/onewebmedia/GradProj_2.pdf"""
        g_a, g_a_variables = self.__getActiveShapeConstraints()

        KM.Logger.PrintInfo("ShapeOpt", "Assemble vector of objective gradient.")
        nabla_f = KM.Vector()
        s = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f, KSO.DF1DX_MAPPED)

        if len(g_a) == 0:
            KM.Logger.PrintInfo("ShapeOpt", "No constraints active, use negative objective gradient as search direction.")
            s = nabla_f * (-1.0)
            s *= self.shape_step_size / s.norm_inf()
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.SEARCH_DIRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, [0.0]*len(s), KSO.CORRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.CONTROL_POINT_UPDATE)
            return


        KM.Logger.PrintInfo("ShapeOpt", "Assemble matrix of constraint gradient.")
        N = KM.Matrix()
        self.optimization_utilities.AssembleMatrix(self.design_surface, N, g_a_variables)

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
            if c.norm_inf() <= self.max_correction_share * self.shape_step_size:
                delta = self.shape_step_size - c.norm_inf()
                s *= delta/s.norm_inf()
            else:
                KM.Logger.PrintWarning("ShapeOpt", f"Correction is scaled down from {c.norm_inf()} to {self.max_correction_share * self.step_size}.")
                c *= self.max_correction_share * self.shape_step_size / c.norm_inf()
                s *= (1.0 - self.max_correction_share) * self.shape_step_size / s.norm_inf()
        else:
            s *= self.shape_step_size / s.norm_inf()

        self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.SEARCH_DIRECTION)
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, c, KSO.CORRECTION)
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, s+c, KSO.CONTROL_POINT_UPDATE)

    # --------------------------------------------------------------------------
    def __getActiveShapeConstraints(self):
        active_constraint_values = []
        active_constraint_variables = []

        for constraint in self.constraints:
            if self.__isConstraintActive(constraint):
                identifier = constraint["identifier"].GetString()
                constraint_value = self.communicator.getStandardizedValue(identifier)
                active_constraint_values.append(constraint_value)
                active_constraint_variables.append(
                    self.shape_constraint_gradient_variables[identifier]["mapped_gradient"])

        return active_constraint_values, active_constraint_variables
# ==============================================================================
