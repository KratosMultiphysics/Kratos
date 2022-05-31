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
from KratosMultiphysics.ShapeOptimizationApplication.algorithms.algorithm_base import OptimizationAlgorithm
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
from KratosMultiphysics.ShapeOptimizationApplication.loggers import data_logger_factory
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_variable_utilities import WriteDictionaryDataOnNodalVariable, WriteListToNodalVariable, ReadNodalVariableToList
from KratosMultiphysics.ShapeOptimizationApplication.utilities import custom_math as cm

# ==============================================================================
class AlgorithmPenalizedProjection(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = KM.Parameters("""
        {
            "name"                    : "penalized_projection",
            "correction_scaling"      : 1.0,
            "use_adaptive_correction" : true,
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
        self.mapper_settings = optimization_settings["design_variables"]["filter"]

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.design_surface = None
        self.mapper = None
        self.data_logger = None
        self.optimization_utilities = None
        self.constraint_value = 0.0
        self.correction_scaling = optimization_settings["optimization_algorithm"]["correction_scaling"].GetDouble()

        self.objectives = optimization_settings["objectives"]
        self.constraints = optimization_settings["constraints"]

        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.SEARCH_DIRECTION)

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Penalized projection algorithm only supports one objective function!")
        if self.constraints.size() == 0:
            raise RuntimeError("Penalized projection algorithm requires definition of a constraint!")
        if self.constraints.size() > 1:
            raise RuntimeError("Penalized projection algorithm only supports one constraint!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.model_part_controller.Initialize()

        self.analyzer.InitializeBeforeOptimizationLoop()

        self.design_surface = self.model_part_controller.GetDesignSurface()

        self.mapper = mapper_factory.CreateMapper(self.design_surface, self.design_surface, self.mapper_settings)
        self.mapper.Initialize()
        self.model_part_controller.InitializeDamping()

        self.data_logger = data_logger_factory.CreateDataLogger(self.model_part_controller, self.communicator, self.optimization_settings)
        self.data_logger.InitializeDataLogging()

        self.optimization_utilities = KSO.OptimizationUtilities

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

            if self.optimization_iteration > 1:
                self.__savePreviousGradientAndUpdate()

            self.__initializeNewShape()

            self.__analyzeShape()

            self.__computeShapeUpdate()

            if self.optimization_iteration > 1:
                self.__computeSensitivityHeatmap()

            self.__logCurrentOptimizationStep()

            KM.Logger.Print("")
            KM.Logger.PrintInfo("ShapeOpt", "Time needed for current optimization step = ", timer.GetLapTime(), "s")
            KM.Logger.PrintInfo("ShapeOpt", "Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

            if self.__isAlgorithmConverged():
                break
            else:
                self.__determineAbsoluteChanges()

    # --------------------------------------------------------------------------
    def __savePreviousGradientAndUpdate(self):
        # save previous search direction and objective gradient

        # previous control update
        self.d_prev_c = []
        # previous shape update
        self.d_prev_x = []
        # previous search direction
        self.prev_s = []
        # previous DF1DX in design space
        self.df_prev_x = []
        # previous DF1DX in control space
        self.df_prev_c = []
        # previous DC1DX in design space
        self.dc_prev_x = []
        # previous DC1DX in control space
        self.dc_prev_c = []

        for node in self.design_surface.Nodes:
            # The following variables are not yet updated and therefore contain the information from the previos step
            self.d_prev_x.append(node.GetSolutionStepValue(KSO.SHAPE_UPDATE))
            self.d_prev_c.append(node.GetSolutionStepValue(KSO.CONTROL_POINT_UPDATE))
            self.prev_s.append(node.GetSolutionStepValue(KSO.SEARCH_DIRECTION))
            self.df_prev_x.append(node.GetSolutionStepValue(KSO.DF1DX))
            self.df_prev_c.append(node.GetSolutionStepValue(KSO.DF1DX_MAPPED))
            self.dc_prev_x.append(node.GetSolutionStepValue(KSO.DC1DX))
            self.dc_prev_c.append(node.GetSolutionStepValue(KSO.DC1DX_MAPPED))

    # --------------------------------------------------------------------------
    def __computeSensitivityHeatmap(self):
        s = []
        df_dx = []
        df_dc = []
        dc_dx = []
        dc_dc = []
        for node in self.design_surface.Nodes:
            s.append(node.GetSolutionStepValue(KSO.SEARCH_DIRECTION))
            df_dx.append(node.GetSolutionStepValue(KSO.DF1DX))
            df_dc.append(node.GetSolutionStepValue(KSO.DF1DX_MAPPED))
            dc_dx.append(node.GetSolutionStepValue(KSO.DC1DX))
            dc_dc.append(node.GetSolutionStepValue(KSO.DC1DX_MAPPED))

        # Heatmap from search direction
        d_s = []
        hessian_diag_s = []
        inv_hessian_diag_s = []
        for i in range(len(self.design_surface.Nodes)):
            delta_s = cm.Minus(s[i], self.prev_s[i])
            if cm.Dot(delta_s, self.d_prev_c[i]) == 0.0:
            # if cm.Dot(delta_s, self.d_prev_c[i]) < 1e6:
                alpha_s = 1e8
            else:
                alpha_s = abs(cm.Dot(delta_s, delta_s) / cm.Dot(delta_s, self.d_prev_c[i]))
            d_s_temp = cm.ScalarVectorProduct(-1/alpha_s, s[i])
            d_s.append(d_s_temp[0])
            d_s.append(d_s_temp[1])
            d_s.append(d_s_temp[2])
            hessian_diag_s.append(alpha_s)
            inv_hessian_diag_s.append(1/alpha_s)

        WriteListToNodalVariable(hessian_diag_s, self.design_surface, KSO.HESSIAN_S, 1)
        WriteListToNodalVariable(inv_hessian_diag_s, self.design_surface, KSO.INV_HESSIAN_S, 1)
        WriteListToNodalVariable(d_s, self.design_surface, KSO.HEATMAP_S, 3)

        # Heatmap DF1DX
        heat_f_x = []
        hessian_diag_f_x = []
        inv_hessian_diag_f_x = []
        for i in range(len(self.design_surface.Nodes)):
            delta_f_x = cm.Minus(df_dx[i], self.df_prev_x[i])
            if cm.Dot(delta_f_x, self.d_prev_x[i]) == 0.0:
            # if cm.Dot(delta_f_x, self.d_prev_x[i]) < 1e6:
                alpha_f_x = 1e8
            else:
                alpha_f_x = abs(cm.Dot(delta_f_x, delta_f_x) / cm.Dot(delta_f_x, self.d_prev_x[i]))
            heat_f_x_i = cm.ScalarVectorProduct(-1/alpha_f_x, df_dx[i])
            heat_f_x.append(heat_f_x_i[0])
            heat_f_x.append(heat_f_x_i[1])
            heat_f_x.append(heat_f_x_i[2])
            hessian_diag_f_x.append(alpha_f_x)
            inv_hessian_diag_f_x.append(1/alpha_f_x)

        WriteListToNodalVariable(hessian_diag_f_x, self.design_surface, KSO.HESSIAN_DF1DX, 1)
        WriteListToNodalVariable(inv_hessian_diag_f_x, self.design_surface, KSO.INV_HESSIAN_DF1DX, 1)
        WriteListToNodalVariable(heat_f_x, self.design_surface, KSO.HEATMAP_DF1DX, 3)

        # Heatmap DF1DX_MAPPED
        heat_f_c = []
        hessian_diag_f_c = []
        inv_hessian_diag_f_c = []
        for i in range(len(self.design_surface.Nodes)):
            delta_f_c = cm.Minus(df_dc[i], self.df_prev_c[i])
            if cm.Dot(delta_f_c, self.d_prev_c[i]) == 0.0:
            # if cm.Dot(delta_f_c, self.d_prev_c[i]) < 1e6:
                alpha_f_c = 1e8
            else:
                alpha_f_c = abs(cm.Dot(delta_f_c, delta_f_c) / cm.Dot(delta_f_c, self.d_prev_c[i]))
            heat_f_c_i = cm.ScalarVectorProduct(-1/alpha_f_c, df_dc[i])
            heat_f_c.append(heat_f_c_i[0])
            heat_f_c.append(heat_f_c_i[1])
            heat_f_c.append(heat_f_c_i[2])
            hessian_diag_f_c.append(alpha_f_c)
            inv_hessian_diag_f_c.append(1/alpha_f_c)

        WriteListToNodalVariable(hessian_diag_f_c, self.design_surface, KSO.HESSIAN_DF1DX_MAPPED, 1)
        WriteListToNodalVariable(inv_hessian_diag_f_c, self.design_surface, KSO.INV_HESSIAN_DF1DX_MAPPED, 1)
        WriteListToNodalVariable(heat_f_c, self.design_surface, KSO.HEATMAP_DF1DX_MAPPED, 3)

        # Heatmap DC1DX
        heat_c_x = []
        hessian_diag_c_x = []
        inv_hessian_diag_c_x = []
        for i in range(len(self.design_surface.Nodes)):
            delta_c_x = cm.Minus(dc_dx[i], self.dc_prev_x[i])
            if cm.Dot(delta_c_x, self.d_prev_x[i]) == 0.0:
            # if cm.Dot(delta_c_x, self.d_prev_x[i]) < 1e6:
                alpha_c_x = 1e8
            else:
                alpha_c_x = abs(cm.Dot(delta_c_x, delta_c_x) / cm.Dot(delta_c_x, self.d_prev_x[i]))
            heat_c_x_i = cm.ScalarVectorProduct(-1/alpha_c_x, dc_dx[i])
            heat_c_x.append(heat_c_x_i[0])
            heat_c_x.append(heat_c_x_i[1])
            heat_c_x.append(heat_c_x_i[2])
            hessian_diag_c_x.append(alpha_c_x)
            inv_hessian_diag_c_x.append(1/alpha_c_x)

        WriteListToNodalVariable(hessian_diag_c_x, self.design_surface, KSO.HESSIAN_DC1DX, 1)
        WriteListToNodalVariable(inv_hessian_diag_c_x, self.design_surface, KSO.INV_HESSIAN_DC1DX, 1)
        WriteListToNodalVariable(heat_c_x, self.design_surface, KSO.HEATMAP_DC1DX, 3)

        # Heatmap DC1DX_MAPPED
        heat_c_c = []
        hessian_diag_c_c = []
        inv_hessian_diag_c_c = []
        for i in range(len(self.design_surface.Nodes)):
            delta_c_c = cm.Minus(dc_dc[i], self.dc_prev_c[i])
            if cm.Dot(delta_c_c, self.d_prev_c[i]) == 0.0:
            # if cm.Dot(delta_c_c, self.d_prev_c[i]) < 1e6:
                alpha_c_c = 1e8
            else:
                alpha_c_c = abs(cm.Dot(delta_c_c, delta_c_c) / cm.Dot(delta_c_c, self.d_prev_c[i]))
            heat_c_c_i = cm.ScalarVectorProduct(-1/alpha_c_c, dc_dc[i])
            heat_c_c.append(heat_c_c_i[0])
            heat_c_c.append(heat_c_c_i[1])
            heat_c_c.append(heat_c_c_i[2])
            hessian_diag_c_c.append(alpha_c_c)
            inv_hessian_diag_c_c.append(1/alpha_c_c)

        WriteListToNodalVariable(hessian_diag_c_c, self.design_surface, KSO.HESSIAN_DC1DX_MAPPED, 1)
        WriteListToNodalVariable(inv_hessian_diag_c_c, self.design_surface, KSO.INV_HESSIAN_DC1DX_MAPPED, 1)
        WriteListToNodalVariable(heat_c_c, self.design_surface, KSO.HEATMAP_DC1DX_MAPPED, 3)

        # Heatmap Max
        heat_max = []
        df_dx = ReadNodalVariableToList(self.design_surface, KSO.DF1DX)

        ### Start: Normalize objective gradient
        df_dx_norm = cm.NormInf3D(df_dx)
        if df_dx_norm != 0.0:
            df_dx_normalized = cm.ScalarVectorProduct(1/df_dx_norm, df_dx)
        else:
            df_dx_normalized = [0] * len(df_dx)
        ### End: Normalize objective gradient

        # ### Start: Normalize objective gradient by objective value
        # objective_value = self.communicator.getStandardizedValue(self.objectives[0]["identifier"].GetString())
        # df_dx_normalized = cm.ScalarVectorProduct(1/objective_value, df_dx)
        # ### End: Normalize objective gradient by objective value

        dc_dx = ReadNodalVariableToList(self.design_surface, KSO.DC1DX)
        ## Start: Normalize constraint gradient
        dc_dx_norm = cm.NormInf3D(dc_dx)
        if dc_dx_norm != 0.0:
            dc_dx_normalized = cm.ScalarVectorProduct(1/dc_dx_norm, dc_dx)
        else:
            dc_dx_normalized = [0] * len(dc_dx)
        ## End: Normalize constraint gradient

        # ### Start: Normalize constraint gradient by constraint value
        # constraint_value = self.communicator.getStandardizedValue(self.constraints[0]["identifier"].GetString())
        # if constraint_value != 0.0:
        #     dc_dx_normalized = cm.ScalarVectorProduct(1/constraint_value, dc_dx)
        # else:
        #     dc_dx_normalized = [0] * len(dc_dx)
        # ### End: Normalize constraint gradient by constraint value

        for i in range(len(self.design_surface.Nodes)):
            df_dx_i = df_dx_normalized[3*i:3*i+3]
            df_dx_i_norm = cm.Norm2(df_dx_i)
            dc_dx_i = dc_dx_normalized[3*i:3*i+3]
            dc_dx_i_norm = cm.Norm2(dc_dx_i)
            heat_max_i = max(df_dx_i_norm, dc_dx_i_norm)
            heat_max.append(heat_max_i)

        ### Start: Normalize heat sens field
        # heat_norm = cm.NormInf3D(heat)
        # heat = cm.ScalarVectorProduct(1/heat_norm, heat)
        ### End: Normalize heat sens field

        WriteListToNodalVariable(heat_max, self.design_surface, KSO.HEATMAP_MAX, 1)
        WriteListToNodalVariable(df_dx_normalized, self.design_surface, KSO.DF1DX_NORMALIZED, 3)
        WriteListToNodalVariable(dc_dx_normalized, self.design_surface, KSO.DC1DX_NORMALIZED, 3)

        # Heatmap Max Mapped
        heat_max_mapped = []
        df_dc = ReadNodalVariableToList(self.design_surface, KSO.DF1DX_MAPPED)

        ### Start: Normalize objective gradient
        df_dc_norm = cm.NormInf3D(df_dc)
        if df_dc_norm != 0.0:
            df_dc_normalized = cm.ScalarVectorProduct(1/df_dc_norm, df_dc)
        else:
            df_dc_normalized = [0] * len(df_dc)
        ### End: Normalize objective gradient

        # ### Start: Normalize objective gradient by objective value
        # objective_value = self.communicator.getStandardizedValue(self.objectives[0]["identifier"].GetString())
        # df_dc_normalized = cm.ScalarVectorProduct(1/objective_value, df_dc)
        # ### End: Normalize objective gradient by objective value

        dc_dc = ReadNodalVariableToList(self.design_surface, KSO.DC1DX_MAPPED)
        ## Start: Normalize constraint gradient
        dc_dc_norm = cm.NormInf3D(dc_dc)
        if dc_dc_norm != 0.0:
            dc_dc_normalized = cm.ScalarVectorProduct(1/dc_dc_norm, dc_dc)
        else:
            dc_dc_normalized = [0] * len(dc_dc)
        ## End: Normalize constraint gradient

        # ### Start: Normalize constraint gradient by constraint value
        # constraint_value = self.communicator.getStandardizedValue(self.constraints[0]["identifier"].GetString())
        # if constraint_value != 0.0:
        #     dc_dc_normalized = cm.ScalarVectorProduct(1/constraint_value, dc_dc)
        # else:
        #     dc_dc_normalized = [0] * len(dc_dc)
        # ### End: Normalize constraint gradient by constraint value

        for i in range(len(self.design_surface.Nodes)):
            df_dc_i = df_dc_normalized[3*i:3*i+3]
            df_dc_i_norm = cm.Norm2(df_dc_i)
            dc_dc_i = dc_dc_normalized[3*i:3*i+3]
            dc_dc_i_norm = cm.Norm2(dc_dc_i)
            heat_max_mapped_i = max(df_dc_i_norm, dc_dc_i_norm)
            heat_max_mapped.append(heat_max_mapped_i)

        ### Start: Normalize heat sens field
        # heat_norm = cm.NormInf3D(heat)
        # heat = cm.ScalarVectorProduct(1/heat_norm, heat)
        ### End: Normalize heat sens field

        WriteListToNodalVariable(heat_max_mapped, self.design_surface, KSO.HEATMAP_MAX_MAPPED, 1)
        WriteListToNodalVariable(df_dc_normalized, self.design_surface, KSO.DF1DX_MAPPED_NORMALIZED, 3)
        WriteListToNodalVariable(dc_dc_normalized, self.design_surface, KSO.DC1DX_MAPPED_NORMALIZED, 3)

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.data_logger.FinalizeDataLogging()
        self.analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeNewShape(self):
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)
        self.model_part_controller.UpdateMeshAccordingInputVariable(KSO.SHAPE_UPDATE)
        self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyzeShape(self):
        self.communicator.initializeCommunication()
        self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
        self.communicator.requestGradientOf(self.objectives[0]["identifier"].GetString())
        self.communicator.requestValueOf(self.constraints[0]["identifier"].GetString())
        self.communicator.requestGradientOf(self.constraints[0]["identifier"].GetString())

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.optimization_iteration, self.communicator)

        objGradientDict = self.communicator.getStandardizedGradient(self.objectives[0]["identifier"].GetString())
        conGradientDict = self.communicator.getStandardizedGradient(self.constraints[0]["identifier"].GetString())

        WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, KSO.DF1DX)
        WriteDictionaryDataOnNodalVariable(conGradientDict, self.optimization_model_part, KSO.DC1DX)

        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool() or self.constraints[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ComputeUnitSurfaceNormals()

        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(KSO.DF1DX)

        if self.constraints[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(KSO.DC1DX)

        self.model_part_controller.DampNodalSensitivityVariableIfSpecified(KSO.DF1DX)
        self.model_part_controller.DampNodalSensitivityVariableIfSpecified(KSO.DC1DX)

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):
        self.mapper.Update()
        self.mapper.InverseMap(KSO.DF1DX, KSO.DF1DX_MAPPED)
        self.mapper.InverseMap(KSO.DC1DX, KSO.DC1DX_MAPPED)
        is_adaptive = self.algorithm_settings["use_adaptive_correction"].GetBool()

        constraint_value = self.communicator.getStandardizedValue(self.constraints[0]["identifier"].GetString())
        if self.__isConstraintActive(constraint_value):
            self.optimization_utilities.ComputeProjectedSearchDirection(self.design_surface)
            self.correction_scaling = self.optimization_utilities.CorrectProjectedSearchDirection(self.design_surface, self.constraint_value, constraint_value, self.correction_scaling, is_adaptive)
            self.constraint_value = constraint_value
        else:
            self.optimization_utilities.ComputeSearchDirectionSteepestDescent(self.design_surface)

        normalize = self.algorithm_settings["line_search"]["normalize_search_direction"].GetBool()
        self.optimization_utilities.ComputeControlPointUpdate(self.design_surface, self.step_size, normalize)

        self.mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
        self.model_part_controller.DampNodalUpdateVariableIfSpecified(KSO.SHAPE_UPDATE)

    # --------------------------------------------------------------------------
    def __isConstraintActive(self, constraintValue):
        if self.constraints[0]["type"].GetString() == "=":
            return True
        elif constraintValue > 0:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        additional_values_to_log = {}
        additional_values_to_log["correction_scaling"] = self.correction_scaling
        additional_values_to_log["step_size"] = self.step_size
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
        self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.CONTROL_POINT_UPDATE, KSO.CONTROL_POINT_CHANGE)
        self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.SHAPE_UPDATE, KSO.SHAPE_CHANGE)

# ==============================================================================
