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
class AlgorithmFreeThicknessOptimizationv3(OptimizationAlgorithm):
    """
        Algorithm for free thickness optimization using a filtering operation which is
        defined on the initial geometry (thus not changing over the course of the optimization)
        and filtering the total thickness update t = T + dt = T + A * dt_control.
        Shape optimization is conducted simultaneously!
    """
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = KM.Parameters("""
        {
            "name"                    : "gradient_projection",
            "max_correction_share"    : 0.75,
            "max_iterations"          : 100,
            "relative_tolerance"      : 1e-3,
            "shape_scaling_divisor"   : 4.0,
            "thickness_scaling_divisor": 10.0,
            "line_search" : {
                "line_search_type"           : "manual_stepping",
                "normalize_search_direction" : true,
                "step_size"                  : 1.0
            }
        }""")
        self.shape_opt = False
        free_thickness_settings = None
        for design_variable in optimization_settings["design_variables"].values():
            if design_variable["type"].GetString() == "vertex_morphing":
                self.shape_opt = True
                self.shape_mapper_settings = design_variable["filter"].Clone()
            elif design_variable["type"].GetString() == "free_thickness":
                free_thickness_settings = design_variable.Clone()

        self.algorithm_settings =  optimization_settings["optimization_algorithm"]
        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

        self.optimization_settings = optimization_settings
        self.initial_mapper_settings = free_thickness_settings["filter"]

        if free_thickness_settings.Has("projection"):
            self.projection = free_thickness_settings["projection"]
        else:
            self.projection = False

        if self.projection:
            self.thickness_targets = free_thickness_settings["projection_settings"]["available_thicknesses"].GetVector()
            self.beta = free_thickness_settings["projection_settings"]["initial_beta"].GetDouble()
            self.q = free_thickness_settings["projection_settings"]["increase_beta_factor"].GetDouble()

        if free_thickness_settings.Has("create_property_copies"):
            self.create_property_copies = free_thickness_settings["create_property_copies"].GetBool()
        else:
            self.create_property_copies = True

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.design_surface = None
        self.shape_mapper = None
        self.data_logger = None
        self.optimization_utilities = None
        self.variable_utils = None

        self.objectives = optimization_settings["objectives"]
        self.constraints = optimization_settings["constraints"].values()

        self.shape_constraint_gradient_variables = {}
        self.constraint_gradient_variables = {}
        for itr, constraint in enumerate(self.constraints):
            self.constraint_gradient_variables.update({
                constraint["identifier"].GetString() : {
                    "gradient": KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT"),
                    "mapped_gradient": KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT_MAPPED"),
                    "projected_gradient": KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT_PROJECTED")
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
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_SEARCH_DIRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CORRECTION)

        self.variable_scaling = False
        if self.shape_opt:
            self.optimization_model_part.AddNodalSolutionStepVariable(KSO.PROJECTION)
            self.optimization_model_part.AddNodalSolutionStepVariable(KSO.CORRECTION)
            self.optimization_model_part.AddNodalSolutionStepVariable(KSO.SEARCH_DIRECTION)
            self.variable_scaling = True

        self.shape_scaling_divisor = self.algorithm_settings["shape_scaling_divisor"].GetDouble()
        self.thickness_scaling_divisor = self.algorithm_settings["thickness_scaling_divisor"].GetDouble()

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Gradient projection algorithm only supports one objective function!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.initial_model_part_controller.Initialize()
        self.model_part_controller.Initialize()

        self.analyzer.InitializeBeforeOptimizationLoop()

        if self.create_property_copies:
            self.model_part_controller.ModifyInitialProperties()

        self.initial_design_surface = self.initial_model_part_controller.GetDesignSurface()
        self.design_surface = self.model_part_controller.GetDesignSurface()

        if self.variable_scaling:
            step_t = ( self.thickness_targets[len(self.thickness_targets)-1] - self.thickness_targets[0] ) / self.thickness_scaling_divisor
            h = KSO.GeometryUtilities(self.design_surface).CalculateAverageElementSize()
            step_x = h / self.shape_scaling_divisor
            self.shape_scaling_factor = step_x
            self.thickness_scaling_factor = step_t
            KM.Logger.PrintInfo(f"Opt", f"Average element size {h}")
            KM.Logger.PrintInfo(f"Opt", f"Shape scaling factor {self.shape_scaling_factor}")
            KM.Logger.PrintInfo(f"Opt", f"Thickness scaling factor {self.thickness_scaling_factor}")

        self.initial_mapper = mapper_factory.CreateMapper(self.initial_design_surface, self.initial_design_surface, self.initial_mapper_settings)
        self.initial_mapper.Initialize()
        self.initial_model_part_controller.InitializeDamping()

        if self.shape_opt:
            self.shape_mapper = mapper_factory.CreateMapper(self.design_surface, self.design_surface, self.shape_mapper_settings)
            self.shape_mapper.Initialize()
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

            self.__updatePhysicalVariables()

            self.__analyze()

            self.__mapGradients()

            self.__computeControlUpdate()

            self.__mapControlUpdate()

            self.__logCurrentOptimizationStep()

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
    def __updatePhysicalVariables(self):

        self.__updateThickness()
        if self.shape_opt:
            self.__updateShape()

    # --------------------------------------------------------------------------
    def __updateThickness(self):
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)
        self.model_part_controller.UpdateThicknessAccordingInitialAndInputVariable(KSO.THICKNESS_CHANGE)

        # update solution step variables
        initial_thickness = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, initial_thickness, KSO.THICKNESS_INITIAL)
        thickness_change = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, thickness_change, KSO.THICKNESS_CHANGE)
        new_thickness = initial_thickness + thickness_change
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, new_thickness, KSO.THICKNESS)

        # visualization on mesh
        for condition in self.optimization_model_part.Conditions:
            condition.SetValue(KSO.THICKNESS_ELEMENTAL, condition.Properties.GetValue(KM.THICKNESS))

    # --------------------------------------------------------------------------
    def __updateShape(self):
        self.model_part_controller.UpdateMeshAccordingInputVariable(KSO.SHAPE_UPDATE)
        self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyze(self):

        self.communicator.initializeCommunication()
        self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
        if self.shape_opt:
            self.communicator.requestGradientOf(self.objectives[0]["identifier"].GetString())
        self.communicator.requestThicknessGradientOf(self.objectives[0]["identifier"].GetString())

        for constraint in self.constraints:
            con_id =  constraint["identifier"].GetString()
            self.communicator.requestValueOf(con_id)
            if self.shape_opt:
                self.communicator.requestGradientOf(con_id)
            self.communicator.requestThicknessGradientOf(con_id)

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.optimization_iteration, self.communicator)

        self.__analyzeThickness()
        if self.shape_opt:
            self.__analyzeShape()

    # --------------------------------------------------------------------------
    def __analyzeThickness(self):

        # project and damp objective gradients
        objElementGradientDict = self.communicator.getStandardizedThicknessGradient(self.objectives[0]["identifier"].GetString())
        self.__mapElementGradientToNode(objElementGradientDict, KSO.DF1DT)

        self.model_part_controller.DampNodalThicknessVariableIfSpecified(KSO.DF1DT)

        self.variable_utils.CopyModelPartNodalVar(KSO.DF1DT, self.model_part_controller.GetOptimizationModelPart(),
                                                  self.initial_model_part_controller.GetOptimizationModelPart(), 0)

        # project and damp constraint gradients
        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            conElementGradientDict = self.communicator.getStandardizedThicknessGradient(con_id)
            gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
            self.__mapElementGradientToNode(conElementGradientDict, gradient_variable)

            self.model_part_controller.DampNodalThicknessVariableIfSpecified(gradient_variable)

            self.variable_utils.CopyModelPartNodalVar(gradient_variable, self.model_part_controller.GetOptimizationModelPart(),
                                                      self.initial_model_part_controller.GetOptimizationModelPart(), 0)

    # --------------------------------------------------------------------------
    def __analyzeShape(self):

        # compute normals only if required
        surface_normals_required = self.objectives[0]["project_gradient_on_surface_normals"].GetBool()
        for constraint in self.constraints:
            if constraint["project_gradient_on_surface_normals"].GetBool():
                surface_normals_required = True

        if surface_normals_required:
            self.model_part_controller.ComputeUnitSurfaceNormals()

        # project and damp objective gradients
        objGradientDict = self.communicator.getStandardizedGradient(self.objectives[0]["identifier"].GetString())
        WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, KSO.DF1DX)
        nabla_f_raw = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f_raw, KSO.DF1DX)

        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(KSO.DF1DX)

        self.model_part_controller.DampNodalSensitivityVariableIfSpecified(KSO.DF1DX)

        nabla_f = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f, KSO.DF1DX)

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
    def __mapGradients(self):

        self.__mapThicknessGradients()
        if self.shape_opt:
            self.__mapShapeGradients()

    # --------------------------------------------------------------------------
    def __mapThicknessGradients(self):

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

    def __mapShapeGradients(self):

        self.shape_mapper.Update()
        self.shape_mapper.InverseMap(KSO.DF1DX, KSO.DF1DX_MAPPED)

        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            gradient_variable = self.shape_constraint_gradient_variables[con_id]["gradient"]
            mapped_gradient_variable = self.shape_constraint_gradient_variables[con_id]["mapped_gradient"]
            self.shape_mapper.InverseMap(gradient_variable, mapped_gradient_variable)

    # --------------------------------------------------------------------------
    def __computeControlUpdate(self):

        number_of_control_variables = len(self.design_surface.Nodes)
        if self.shape_opt:
            number_of_control_variables *= 4

        prev_s = KM.Vector(number_of_control_variables)

        max_inner_iter = 10
        for inner_iter in range(max_inner_iter):

            nabla_f = self.__assembleObjectiveGradient()

            g_a, g_a_t_variables, g_a_x_variables = self.__getActiveConstraints()

            s = KM.Vector()
            corr = KM.Vector()

            if len(g_a) == 0:
                KM.Logger.PrintInfo("ShapeOpt", "No constraints active, use negative objective gradient as search direction.")
                s = nabla_f * (-1.0)

                if self.projection:
                    s = self.__ProjectSearchDirectionAndGradients(s, g_a_t_variables)

                if not self.projection or (s - prev_s).norm_inf() == 0.0 or inner_iter == max_inner_iter - 1:
                    if s.norm_inf() > 0:
                        s *= self.step_size / s.norm_inf()
                    corr = 0.0 * s
                    self.__saveAlgorithmResultsToVariables(s, corr)
                    return
                else:
                    prev_s = s
                    continue

            KM.Logger.PrintInfo("ShapeOpt", "Assemble matrix of constraint gradient.")
            g_a_gradients = self.__assembleConstraintGradients(g_a_t_variables, g_a_x_variables)

            N = KM.Matrix()
            self.optimization_utilities.AssembleMatrixFromGradientVectors(self.design_surface, N, g_a_gradients)

            settings = KM.Parameters('{ "solver_type" : "LinearSolversApplication.dense_col_piv_householder_qr" }')
            solver = dense_linear_solver_factory.ConstructSolver(settings)

            KM.Logger.PrintInfo("ShapeOpt", "Calculate projected search direction and correction.")
            self.optimization_utilities.CalculateProjectedSearchDirectionAndCorrection(
                nabla_f,
                N,
                g_a,
                solver,
                s,
                corr)

            if self.projection:
                s = self.__ProjectSearchDirectionAndGradients(s, g_a_t_variables)

            if not self.projection or (s - prev_s).norm_inf() == 0.0 or inner_iter == max_inner_iter - 1:
                if corr.norm_inf() != 0.0:
                    if corr.norm_inf() <= self.max_correction_share * self.step_size:
                        delta = self.step_size - corr.norm_inf()
                        s *= delta/s.norm_inf()
                    else:
                        KM.Logger.PrintWarning("ShapeOpt", f"Correction is scaled down from {corr.norm_inf()} to {self.max_correction_share * self.step_size}.")
                        corr *= self.max_correction_share * self.step_size / corr.norm_inf()
                        s *= (1.0 - self.max_correction_share) * self.step_size / s.norm_inf()
                else:
                    if s.norm_inf() > 0:
                        s *= self.step_size / s.norm_inf()

                self.__saveAlgorithmResultsToVariables(s, corr)
                return
            else:
                prev_s = s
                continue

    # --------------------------------------------------------------------------
    def __assembleObjectiveGradient(self):

        df_dt = KM.Vector()

        if self.projection:
            self.optimization_utilities.AssembleVector(self.design_surface, df_dt, KSO.DF1DT_PROJECTED)
        else:
            self.optimization_utilities.AssembleVector(self.design_surface, df_dt, KSO.DF1DT_MAPPED)

        if self.shape_opt:
            df_dx = KM.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, df_dx, KSO.DF1DX_MAPPED)

            df_dc = KM.Vector(df_dt.Size()+df_dx.Size())

            if self.variable_scaling:
                df_dt *= self.thickness_scaling_factor
                df_dx *= self.shape_scaling_factor

            for i in range(df_dt.Size()):
                df_dc[i] = df_dt[i]

            shift = df_dt.Size()
            for i in range(df_dx.Size()):
                df_dc[shift+i] = df_dx[i]

            return df_dc

        else:
            return df_dt

    # --------------------------------------------------------------------------
    def __assembleConstraintGradients(self, g_a_t_variables, g_a_x_variables):

        constraint_gradients = []
        for iter, g_a_t_variable in enumerate(g_a_t_variables):
            dg_dt = KM.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, dg_dt, g_a_t_variable)

            if self.shape_opt:
                g_a_x_variable = g_a_x_variables[iter]
                dg_dx = KM.Vector()
                self.optimization_utilities.AssembleVector(self.design_surface, dg_dx, g_a_x_variable)
                dg_dc = KM.Vector(dg_dt.Size()+dg_dx.Size())

                if self.variable_scaling:
                    dg_dt *= self.thickness_scaling_factor
                    dg_dx *= self.shape_scaling_factor

                for i in range(dg_dt.Size()):
                    dg_dc[i] = dg_dt[i]

                shift = dg_dt.Size()
                for i in range(dg_dx.Size()):
                    dg_dc[shift+i] = dg_dx[i]

                constraint_gradients.append(dg_dc)
            else:
                constraint_gradients.append(dg_dt)

        return constraint_gradients

    # --------------------------------------------------------------------------
    def __saveAlgorithmResultsToVariables(self, search_direction, correction):

        s_t = KM.Vector(len(self.design_surface.Nodes))
        c_t = KM.Vector(len(self.design_surface.Nodes))

        for i in range(len(self.design_surface.Nodes)):
            s_t[i] = search_direction[i]
            c_t[i] = correction[i]

        delta_t_control = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, delta_t_control, KSO.THICKNESS_CHANGE_CONTROL)

        self.optimization_utilities.AssignVectorToVariable(self.design_surface, s_t, KSO.THICKNESS_SEARCH_DIRECTION)
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, c_t, KSO.THICKNESS_CORRECTION)
        delta_t_control_update = s_t+c_t
        if self.variable_scaling:
            delta_t_control_update *= self.thickness_scaling_factor
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, delta_t_control_update, KSO.THICKNESS_CONTROL_UPDATE)
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, delta_t_control+delta_t_control_update, KSO.THICKNESS_CHANGE_CONTROL)

        if self.shape_opt:
            shift = len(self.design_surface.Nodes)
            s_x = KM.Vector(3*len(self.design_surface.Nodes))
            c_x = KM.Vector(3*len(self.design_surface.Nodes))

            for i in range(3*len(self.design_surface.Nodes)):
                s_x[i] = search_direction[shift+i]
                c_x[i] = correction[shift+i]

            control_x_update = s_x+c_x
            if self.variable_scaling:
                control_x_update *= self.shape_scaling_factor

            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s_x, KSO.SEARCH_DIRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, c_x, KSO.CORRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, control_x_update, KSO.CONTROL_POINT_UPDATE)

    # --------------------------------------------------------------------------
    def __mapControlUpdate(self):

        self.__mapControlThicknessUpdate()
        if self.shape_opt:
            self.__mapControlShapeUpdate()

    # --------------------------------------------------------------------------
    def __mapControlThicknessUpdate(self):

        ### 4. Project control thickness to obtain filtered thickness
        self.__ProjectThickness()

        ### 5. Filter thickness to obtain physical thickness
        if self.projection:
            self.variable_utils.CopyModelPartNodalVar(KSO.THICKNESS_CHANGE_CONTROL_PROJECTED, self.model_part_controller.GetOptimizationModelPart(),
                                                      self.initial_model_part_controller.GetOptimizationModelPart(), 0)
            self.initial_mapper.Map(KSO.THICKNESS_CHANGE_CONTROL_PROJECTED, KSO.THICKNESS_CHANGE)

        else:
            self.variable_utils.CopyModelPartNodalVar(KSO.THICKNESS_CHANGE_CONTROL, self.model_part_controller.GetOptimizationModelPart(),
                                                      self.initial_model_part_controller.GetOptimizationModelPart(), 0)
            self.initial_mapper.Map(KSO.THICKNESS_CHANGE_CONTROL, KSO.THICKNESS_CHANGE)

        # damping for thickness
        self.initial_model_part_controller.DampNodalThicknessVariableIfSpecified(KSO.THICKNESS_CHANGE)

        self.variable_utils.CopyModelPartNodalVar(KSO.THICKNESS_CHANGE, self.initial_model_part_controller.GetOptimizationModelPart(),
                                                  self.model_part_controller.GetOptimizationModelPart(), 0)


    # --------------------------------------------------------------------------
    def __mapControlShapeUpdate(self):

        self.shape_mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
        self.model_part_controller.DampNodalUpdateVariableIfSpecified(KSO.SHAPE_UPDATE)

    # --------------------------------------------------------------------------
    def __getActiveConstraints(self):
        active_constraint_values = []
        active_constraint_thickness_variables = []
        active_constraint_shape_variables = []

        for constraint in self.constraints:
            if self.__isConstraintActive(constraint):
                identifier = constraint["identifier"].GetString()
                constraint_value = self.communicator.getStandardizedValue(identifier)
                active_constraint_values.append(constraint_value)
                if self.projection:
                    active_constraint_thickness_variables.append(
                        self.constraint_gradient_variables[identifier]["projected_gradient"])
                else:
                    active_constraint_thickness_variables.append(
                        self.constraint_gradient_variables[identifier]["mapped_gradient"])

                if self.shape_opt:
                    active_constraint_shape_variables.append(
                        self.shape_constraint_gradient_variables[identifier]["mapped_gradient"])

        return active_constraint_values, active_constraint_thickness_variables, active_constraint_shape_variables

    # --------------------------------------------------------------------------
    def __isConstraintActive(self, constraint):
        identifier = constraint["identifier"].GetString()
        constraint_value = self.communicator.getStandardizedValue(identifier)
        if constraint["type"].GetString() == "=" or constraint_value >= 0:
            if self.projection:
                thickness_gradient_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(
                    self.design_surface, self.constraint_gradient_variables[identifier]["projected_gradient"]
                )
            else:
                thickness_gradient_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(
                    self.design_surface, self.constraint_gradient_variables[identifier]["mapped_gradient"]
                )
            shape_gradient_norm = 0.0
            if self.shape_opt:
                shape_gradient_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(
                    self.design_surface, self.shape_constraint_gradient_variables[identifier]["mapped_gradient"]
                )
            if math.isclose(thickness_gradient_norm, 0.0, abs_tol=1e-16) and math.isclose(shape_gradient_norm, 0.0, abs_tol=1e-16):
                KM.Logger.PrintWarning("Opt", f"Gradient for constraint {identifier} is 0.0 - will not be considered!")
                return False
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):

        for node in self.design_surface.Nodes:
            initial_thickness = node.GetSolutionStepValue(KSO.THICKNESS_INITIAL)
            thickness_change = node.GetSolutionStepValue(KSO.THICKNESS_CHANGE)
            old_thickness = node.GetSolutionStepValue(KSO.THICKNESS)

            new_thickness = initial_thickness + thickness_change
            thickness_update = new_thickness - old_thickness

            node.SetSolutionStepValue(KSO.THICKNESS_UPDATE, 0, thickness_update)

        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.step_size
        additional_values_to_log["inf_norm_s"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_SEARCH_DIRECTION)
        additional_values_to_log["inf_norm_c"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.THICKNESS_CORRECTION)

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
    def __updateBeta(self):

        if self.projection:
            self.beta *= self.q

    # --------------------------------------------------------------------------
    def __determineAbsoluteShapeChanges(self):
        self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.CONTROL_POINT_UPDATE, KSO.CONTROL_POINT_CHANGE)
        self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.SHAPE_UPDATE, KSO.SHAPE_CHANGE)

    # --------------------------------------------------------------------------
    def __mapElementGradientToNode(self, element_gradient_dict, gradient_variable):

        # reset variables
        for node in self.design_surface.Nodes:
            node.SetSolutionStepValue(gradient_variable, 0, 0.0)

        total_node_areas = dict()
        for condition in self.design_surface.Conditions:
            df_dt = element_gradient_dict[condition.Id]
            for node in condition.GetNodes():
                if node.Id in total_node_areas:
                    total_node_areas[node.Id] += condition.GetGeometry().Area()
                else:
                    total_node_areas[node.Id] = condition.GetGeometry().Area()
                df_dt_node = node.GetSolutionStepValue(gradient_variable)
                df_dt_node += df_dt * condition.GetGeometry().Area()
                node.SetSolutionStepValue(gradient_variable, 0, df_dt_node)

        for node in self.design_surface.Nodes:
            if node.Id in total_node_areas:
                total_node_area = total_node_areas[node.Id]
                df_dt_node = node.GetSolutionStepValue(gradient_variable)
                df_dt_node /= total_node_area
                node.SetSolutionStepValue(gradient_variable, 0, df_dt_node)

    # --------------------------------------------------------------------------
    def __InitializeThicknessField(self):

        element_thicknesses = dict()
        for condition in self.optimization_model_part.Conditions:
            element_thicknesses[condition.Id] = condition.Properties.GetValue(KM.THICKNESS)
            condition.SetValue(KSO.THICKNESS_ELEMENTAL_INITIAL, element_thicknesses[condition.Id])

        self.__mapElementGradientToNode(element_thicknesses, KSO.THICKNESS)
        self.__mapElementGradientToNode(element_thicknesses, KSO.THICKNESS_INITIAL)

    # --------------------------------------------------------------------------
    def __ProjectGradients(self, mapped_gradient_variable, projected_gradient_variable):

        if not self.projection:
            return

        for node in self.design_surface.Nodes:
            delta_t_control = node.GetSolutionStepValue(KSO.THICKNESS_CHANGE_CONTROL)
            mapped_gradient = node.GetSolutionStepValue(mapped_gradient_variable)

            delta_t_m = self.___GetInterval(node, delta_t_control)

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
                search_direction[i] = 0.0
                nabla_f[i] = 0.0
                for constraint_gradient in g_a.values():
                    constraint_gradient[i] = 0.0
            elif t >= delta_t_targets[-1] and search_direction[i] > 0.0:
                search_direction[i] = 0.0
                nabla_f[i] = 0.0
                for constraint_gradient in g_a.values():
                    constraint_gradient[i] = 0.0

        self.optimization_utilities.AssignVectorToVariable(self.design_surface, nabla_f, KSO.DF1DT_PROJECTED)
        for g_a_variable, constraint_gradient in g_a.items():
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, constraint_gradient, g_a_variable)

        return search_direction

    def __ProjectThickness(self):

        if not self.projection:
            return

        for node in self.design_surface.Nodes:
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

# ==============================================================================
