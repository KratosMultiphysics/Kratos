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
from KratosMultiphysics.LinearSolversApplication import dense_linear_solver_factory

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication.algorithms.algorithm_base import OptimizationAlgorithm
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
from KratosMultiphysics.ShapeOptimizationApplication.loggers import data_logger_factory
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_variable_utilities import WriteDictionaryDataOnNodalVariable, WriteListToNodalVariable, ReadNodalVariableToList
from KratosMultiphysics.ShapeOptimizationApplication.utilities import custom_math as cm

import math

# ==============================================================================
class AlgorithmGradientProjection(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = KM.Parameters("""
        {
            "name"                    : "penalized_projection",
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
        self.mapper_settings = optimization_settings["design_variables"]["filter"]

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
                    "gradient": KM.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX"),
                    "mapped_gradient": KM.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")
                }
            })
        self.max_correction_share = self.algorithm_settings["max_correction_share"].GetDouble()

        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.SEARCH_DIRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.CORRECTION)

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
        self.dc1_prev_x = []
        # previous DC1DX in control space
        self.dc1_prev_c = []
        # previous DC2DX in design space
        self.dc2_prev_x = []
        # previous DC2DX in control space
        self.dc2_prev_c = []
        # previous DC3DX in design space
        self.dc3_prev_x = []
        # previous DC3DX in control space
        self.dc3_prev_c = []


        for node in self.design_surface.Nodes:
            # The following variables are not yet updated and therefore contain the information from the previos step
            self.d_prev_x.append(node.GetSolutionStepValue(KSO.SHAPE_UPDATE))
            self.d_prev_c.append(node.GetSolutionStepValue(KSO.CONTROL_POINT_UPDATE))
            self.prev_s.append(node.GetSolutionStepValue(KSO.SEARCH_DIRECTION))
            self.df_prev_x.append(node.GetSolutionStepValue(KSO.DF1DX))
            self.df_prev_c.append(node.GetSolutionStepValue(KSO.DF1DX_MAPPED))
            self.dc1_prev_x.append(node.GetSolutionStepValue(KSO.DC1DX))
            self.dc1_prev_c.append(node.GetSolutionStepValue(KSO.DC1DX_MAPPED))
            self.dc2_prev_x.append(node.GetSolutionStepValue(KSO.DC2DX))
            self.dc2_prev_c.append(node.GetSolutionStepValue(KSO.DC2DX_MAPPED))
            self.dc3_prev_x.append(node.GetSolutionStepValue(KSO.DC3DX))
            self.dc3_prev_c.append(node.GetSolutionStepValue(KSO.DC3DX_MAPPED))

    # --------------------------------------------------------------------------
    def __computeSensitivityHeatmap(self):
        s = []
        # df_dx = []
        # df_dc = []
        # dc1_dx = []
        # dc1_dc = []
        # dc2_dx = []
        # dc2_dc = []
        # dc3_dx = []
        # dc3_dc = []
        for node in self.design_surface.Nodes:
            s.append(node.GetSolutionStepValue(KSO.SEARCH_DIRECTION))
            # df_dx.append(node.GetSolutionStepValue(KSO.DF1DX))
            # df_dc.append(node.GetSolutionStepValue(KSO.DF1DX_MAPPED))
            # dc1_dx.append(node.GetSolutionStepValue(KSO.DC1DX))
            # dc1_dc.append(node.GetSolutionStepValue(KSO.DC1DX_MAPPED))
            # dc2_dx.append(node.GetSolutionStepValue(KSO.DC2DX))
            # dc2_dc.append(node.GetSolutionStepValue(KSO.DC2DX_MAPPED))
            # dc3_dx.append(node.GetSolutionStepValue(KSO.DC3DX))
            # dc3_dc.append(node.GetSolutionStepValue(KSO.DC3DX_MAPPED))

        # Heatmap from search direction
        d_s = []
        hessian_diag_s = []
        max_step = 10000 * self.step_size
        min_step = 0.0001 * self.step_size
        inv_hessian_diag_s = []
        for i in range(len(self.design_surface.Nodes)):
            y_i = cm.Minus(self.prev_s[i], s[i])
            d_i = self.d_prev_c[i]
            if cm.Dot(y_i, y_i) < 1e-9:
                inv_hessian_i = max_step
            else:
                inv_hessian_i = abs(cm.Dot(d_i, y_i) / cm.Dot(y_i, y_i))

            if inv_hessian_i > max_step:
                inv_hessian_i = max_step
            if inv_hessian_i < min_step:
                inv_hessian_i = min_step
            d_s_temp = cm.ScalarVectorProduct(-inv_hessian_i, s[i])
            d_s.append(d_s_temp[0])
            d_s.append(d_s_temp[1])
            d_s.append(d_s_temp[2])
            inv_hessian_diag_s.append(inv_hessian_i)
            hessian_diag_s.append(1/inv_hessian_i)

        prev_inv_hessian_diag_s = KM.Vector()
        self.optimization_utilities.AssembleScalar(self.design_surface, prev_inv_hessian_diag_s, KSO.INV_HESSIAN_S)
        WriteListToNodalVariable(hessian_diag_s, self.design_surface, KSO.HESSIAN_S, 1)
        WriteListToNodalVariable(inv_hessian_diag_s, self.design_surface, KSO.INV_HESSIAN_S, 1)
        WriteListToNodalVariable(d_s, self.design_surface, KSO.HEATMAP_S, 3)

        # Heatmap Max
        heat_max = []
        df_dx = ReadNodalVariableToList(self.design_surface, KSO.DF1DX)

        df_dx_norm = cm.NormInf3D(df_dx)
        if df_dx_norm != 0.0:
            df_dx_normalized = cm.ScalarVectorProduct(1/df_dx_norm, df_dx)
        else:
            df_dx_normalized = [0] * len(df_dx)

        dc1_dx = ReadNodalVariableToList(self.design_surface, KSO.DC1DX)
        dc1_dx_norm = cm.NormInf3D(dc1_dx)
        if dc1_dx_norm != 0.0:
            dc1_dx_normalized = cm.ScalarVectorProduct(1/dc1_dx_norm, dc1_dx)
        else:
            dc1_dx_normalized = [0] * len(dc1_dx)

        dc2_dx = ReadNodalVariableToList(self.design_surface, KSO.DC2DX)
        dc2_dx_norm = cm.NormInf3D(dc2_dx)
        if dc2_dx_norm != 0.0:
            dc2_dx_normalized = cm.ScalarVectorProduct(1/dc2_dx_norm, dc2_dx)
        else:
            dc2_dx_normalized = [0] * len(dc2_dx)

        dc3_dx = ReadNodalVariableToList(self.design_surface, KSO.DC3DX)
        dc3_dx_norm = cm.NormInf3D(dc3_dx)
        if dc3_dx_norm != 0.0:
            dc3_dx_normalized = cm.ScalarVectorProduct(1/dc3_dx_norm, dc3_dx)
        else:
            dc3_dx_normalized = [0] * len(dc3_dx)

        for i in range(len(self.design_surface.Nodes)):
            df_dx_i = df_dx_normalized[3*i:3*i+3]
            df_dx_i_norm = cm.Norm2(df_dx_i)
            dc1_dx_i = dc1_dx_normalized[3*i:3*i+3]
            dc1_dx_i_norm = cm.Norm2(dc1_dx_i)
            dc2_dx_i = dc2_dx_normalized[3*i:3*i+3]
            dc2_dx_i_norm = cm.Norm2(dc2_dx_i)
            dc3_dx_i = dc3_dx_normalized[3*i:3*i+3]
            dc3_dx_i_norm = cm.Norm2(dc3_dx_i)

            heat_max_i = max(df_dx_i_norm, dc1_dx_i_norm, dc2_dx_i_norm, dc3_dx_i_norm)
            heat_max.append(heat_max_i)

        # Heatmap Max Relaxed
        prev_heat_max = KM.Vector()
        self.optimization_utilities.AssembleScalar(self.design_surface, prev_heat_max, KSO.HEATMAP_MAX)
        heat_max_relaxed = KM.Vector()
        self.optimization_utilities.AssembleScalar(self.design_surface, heat_max_relaxed, KSO.HEATMAP_MAX_RELAXED)
        relax_coeff = 0.5
        for i in range(len(self.design_surface.Nodes)):
            heat_max_relaxed[i] = relax_coeff * heat_max[i] + (1 - relax_coeff) * prev_heat_max[i]

        self.optimization_utilities.AssignScalarToVariable(self.design_surface, heat_max_relaxed, KSO.HEATMAP_MAX_RELAXED)
        WriteListToNodalVariable(heat_max, self.design_surface, KSO.HEATMAP_MAX, 1)

        # Heatmap Max Mapped
        heat_max_mapped = []
        df_dc = ReadNodalVariableToList(self.design_surface, KSO.DF1DX_MAPPED)

        df_dc_norm = cm.NormInf3D(df_dc)
        if df_dc_norm != 0.0:
            df_dc_normalized = cm.ScalarVectorProduct(1/df_dc_norm, df_dc)
        else:
            df_dc_normalized = [0] * len(df_dc)

        dc1_dc = ReadNodalVariableToList(self.design_surface, KSO.DC1DX_MAPPED)
        dc1_dc_norm = cm.NormInf3D(dc1_dc)
        if dc1_dc_norm != 0.0:
            dc1_dc_normalized = cm.ScalarVectorProduct(1/dc1_dc_norm, dc1_dc)
        else:
            dc1_dc_normalized = [0] * len(dc1_dc)

        dc2_dc = ReadNodalVariableToList(self.design_surface, KSO.DC2DX_MAPPED)
        dc2_dc_norm = cm.NormInf3D(dc2_dc)
        if dc2_dc_norm != 0.0:
            dc2_dc_normalized = cm.ScalarVectorProduct(1/dc2_dc_norm, dc2_dc)
        else:
            dc2_dc_normalized = [0] * len(dc2_dc)

        dc3_dc = ReadNodalVariableToList(self.design_surface, KSO.DC3DX_MAPPED)
        dc3_dc_norm = cm.NormInf3D(dc3_dc)
        if dc3_dc_norm != 0.0:
            dc3_dc_normalized = cm.ScalarVectorProduct(1/dc3_dc_norm, dc3_dc)
        else:
            dc3_dc_normalized = [0] * len(dc3_dc)

        for i in range(len(self.design_surface.Nodes)):
            df_dc_i = df_dc_normalized[3*i:3*i+3]
            df_dc_i_norm = cm.Norm2(df_dc_i)
            dc1_dc_i = dc1_dc_normalized[3*i:3*i+3]
            dc1_dc_i_norm = cm.Norm2(dc1_dc_i)
            dc2_dc_i = dc2_dc_normalized[3*i:3*i+3]
            dc2_dc_i_norm = cm.Norm2(dc2_dc_i)
            dc3_dc_i = dc3_dc_normalized[3*i:3*i+3]
            dc3_dc_i_norm = cm.Norm2(dc3_dc_i)

            heat_max_mapped_i = max(df_dc_i_norm, dc1_dc_i_norm, dc2_dc_i_norm, dc3_dc_i_norm)
            heat_max_mapped.append(heat_max_mapped_i)

        # Heatmap Max Mapped Relaxed
        prev_heat_max_mapped = KM.Vector()
        self.optimization_utilities.AssembleScalar(self.design_surface, prev_heat_max_mapped, KSO.HEATMAP_MAX_MAPPED)
        heat_max_mapped_relaxed = KM.Vector()
        self.optimization_utilities.AssembleScalar(self.design_surface, heat_max_mapped_relaxed, KSO.HEATMAP_MAX_MAPPED_RELAXED)
        relax_coeff = 0.5
        for i in range(len(self.design_surface.Nodes)):
            heat_max_mapped_relaxed[i] = relax_coeff * heat_max_mapped[i] + (1 - relax_coeff) * prev_heat_max_mapped[i]

        self.optimization_utilities.AssignScalarToVariable(self.design_surface, heat_max_mapped_relaxed, KSO.HEATMAP_MAX_MAPPED_RELAXED)
        WriteListToNodalVariable(heat_max_mapped, self.design_surface, KSO.HEATMAP_MAX_MAPPED, 1)

        # Heatmap L2 Mapped
        heat_L2_mapped = []

        df_dc = ReadNodalVariableToList(self.design_surface, KSO.DF1DX_MAPPED)
        df_dc_norm = cm.Norm2(df_dc)
        if df_dc_norm != 0.0:
            df_dc_normalized = cm.ScalarVectorProduct(1/df_dc_norm, df_dc)
        else:
            df_dc_normalized = [0] * len(df_dc)

        dc1_dc = ReadNodalVariableToList(self.design_surface, KSO.DC1DX_MAPPED)
        dc1_dc_norm = cm.Norm2(dc1_dc)
        if dc1_dc_norm != 0.0:
            dc1_dc_normalized = cm.ScalarVectorProduct(1/dc1_dc_norm, dc1_dc)
        else:
            dc1_dc_normalized = [0] * len(dc1_dc)

        dc2_dc = ReadNodalVariableToList(self.design_surface, KSO.DC2DX_MAPPED)
        dc2_dc_norm = cm.Norm2(dc2_dc)
        if dc2_dc_norm != 0.0:
            dc2_dc_normalized = cm.ScalarVectorProduct(1/dc2_dc_norm, dc2_dc)
        else:
            dc2_dc_normalized = [0] * len(dc2_dc)

        dc3_dc = ReadNodalVariableToList(self.design_surface, KSO.DC3DX_MAPPED)
        dc3_dc_norm = cm.Norm2(dc3_dc)
        if dc3_dc_norm != 0.0:
            dc3_dc_normalized = cm.ScalarVectorProduct(1/dc3_dc_norm, dc3_dc)
        else:
            dc3_dc_normalized = [0] * len(dc3_dc)

        for i in range(len(self.design_surface.Nodes)):
            df_dc_i = df_dc_normalized[3*i:3*i+3]
            df_dc_i_norm = cm.Norm2(df_dc_i)
            dc1_dc_i = dc1_dc_normalized[3*i:3*i+3]
            dc1_dc_i_norm = cm.Norm2(dc1_dc_i)
            dc2_dc_i = dc2_dc_normalized[3*i:3*i+3]
            dc2_dc_i_norm = cm.Norm2(dc2_dc_i)
            dc3_dc_i = dc3_dc_normalized[3*i:3*i+3]
            dc3_dc_i_norm = cm.Norm2(dc3_dc_i)

            heat_L2_mapped_i = max(df_dc_i_norm, dc1_dc_i_norm, dc2_dc_i_norm, dc3_dc_i_norm)
            heat_L2_mapped.append(heat_L2_mapped_i)

        # Heatmap L2 Mapped Relaxed
        prev_heat_L2_mapped = KM.Vector()
        self.optimization_utilities.AssembleScalar(self.design_surface, prev_heat_L2_mapped, KSO.HEATMAP_L2_MAPPED)
        heat_L2_mapped_relaxed = KM.Vector()
        self.optimization_utilities.AssembleScalar(self.design_surface, heat_L2_mapped_relaxed, KSO.HEATMAP_L2_MAPPED_RELAXED)
        relax_coeff = 0.5
        for i in range(len(self.design_surface.Nodes)):
            heat_L2_mapped_relaxed[i] = relax_coeff * heat_L2_mapped[i] + (1 - relax_coeff) * prev_heat_L2_mapped[i]

        self.optimization_utilities.AssignScalarToVariable(self.design_surface, heat_L2_mapped_relaxed, KSO.HEATMAP_L2_MAPPED_RELAXED)
        WriteListToNodalVariable(heat_L2_mapped, self.design_surface, KSO.HEATMAP_L2_MAPPED, 1)

        # Heatmap Value Norm Mapped
        f = self.communicator.getStandardizedValue(self.objectives[0]["identifier"].GetString())
        c1 = self.communicator.getStandardizedValue(self.constraints[0]["identifier"].GetString())
        c2 = self.communicator.getStandardizedValue(self.constraints[1]["identifier"].GetString())
        c3 = self.communicator.getStandardizedValue(self.constraints[2]["identifier"].GetString())

        heat_value_mapped = []

        df_dc = ReadNodalVariableToList(self.design_surface, KSO.DF1DX_MAPPED)
        if f != 0.0:
            df_dc_normalized = cm.ScalarVectorProduct(abs(1/f), df_dc)
        else:
            df_dc_normalized = [0] * len(df_dc)

        dc1_dc = ReadNodalVariableToList(self.design_surface, KSO.DC1DX_MAPPED)
        if c1 != 0.0:
            dc1_dc_normalized = cm.ScalarVectorProduct(abs(1/c1), dc1_dc)
        else:
            dc1_dc_normalized = [0] * len(dc1_dc)

        dc2_dc = ReadNodalVariableToList(self.design_surface, KSO.DC2DX_MAPPED)
        if c2 != 0.0:
            dc2_dc_normalized = cm.ScalarVectorProduct(abs(1/c2), dc2_dc)
        else:
            dc2_dc_normalized = [0] * len(dc2_dc)

        dc3_dc = ReadNodalVariableToList(self.design_surface, KSO.DC3DX_MAPPED)
        if c3 != 0.0:
            dc3_dc_normalized = cm.ScalarVectorProduct(abs(1/c3), dc3_dc)
        else:
            dc3_dc_normalized = [0] * len(dc3_dc)

        for i in range(len(self.design_surface.Nodes)):
            df_dc_i = df_dc_normalized[3*i:3*i+3]
            df_dc_i_norm = cm.Norm2(df_dc_i)
            dc1_dc_i = dc1_dc_normalized[3*i:3*i+3]
            dc1_dc_i_norm = cm.Norm2(dc1_dc_i)
            dc2_dc_i = dc2_dc_normalized[3*i:3*i+3]
            dc2_dc_i_norm = cm.Norm2(dc2_dc_i)
            dc3_dc_i = dc3_dc_normalized[3*i:3*i+3]
            dc3_dc_i_norm = cm.Norm2(dc3_dc_i)

            heat_value_mapped_i = max(df_dc_i_norm, dc1_dc_i_norm, dc2_dc_i_norm, dc3_dc_i_norm)
            heat_value_mapped.append(heat_value_mapped_i)

        # Heatmap Max Mapped Relaxed
        prev_heat_value_mapped = KM.Vector()
        self.optimization_utilities.AssembleScalar(self.design_surface, prev_heat_value_mapped, KSO.HEATMAP_VALUE_MAPPED)
        heat_value_mapped_relaxed = KM.Vector()
        self.optimization_utilities.AssembleScalar(self.design_surface, heat_value_mapped_relaxed, KSO.HEATMAP_VALUE_MAPPED_RELAXED)
        relax_coeff = 0.5
        for i in range(len(self.design_surface.Nodes)):
            heat_value_mapped_relaxed[i] = relax_coeff * heat_value_mapped[i] + (1 - relax_coeff) * prev_heat_value_mapped[i]

        self.optimization_utilities.AssignScalarToVariable(self.design_surface, heat_value_mapped_relaxed, KSO.HEATMAP_VALUE_MAPPED_RELAXED)
        WriteListToNodalVariable(heat_value_mapped, self.design_surface, KSO.HEATMAP_VALUE_MAPPED, 1)

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

        for constraint in self.constraints:
            con_id =  constraint["identifier"].GetString()
            self.communicator.requestValueOf(con_id)
            self.communicator.requestGradientOf(con_id)

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.optimization_iteration, self.communicator)

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

        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(KSO.DF1DX)

        self.model_part_controller.DampNodalSensitivityVariableIfSpecified(KSO.DF1DX)

        # project and damp constraint gradients
        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            conGradientDict = self.communicator.getStandardizedGradient(con_id)
            gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
            WriteDictionaryDataOnNodalVariable(conGradientDict, self.optimization_model_part, gradient_variable)

            if constraint["project_gradient_on_surface_normals"].GetBool():
                self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(gradient_variable)

            self.model_part_controller.DampNodalSensitivityVariableIfSpecified(gradient_variable)

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):
        self.mapper.Update()
        self.mapper.InverseMap(KSO.DF1DX, KSO.DF1DX_MAPPED)

        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
            mapped_gradient_variable = self.constraint_gradient_variables[con_id]["mapped_gradient"]
            self.mapper.InverseMap(gradient_variable, mapped_gradient_variable)

        self.__computeControlPointUpdate()

        self.mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
        self.model_part_controller.DampNodalUpdateVariableIfSpecified(KSO.SHAPE_UPDATE)

    # --------------------------------------------------------------------------
    def __computeControlPointUpdate(self):
        """adapted from https://msulaiman.org/onewebmedia/GradProj_2.pdf"""
        g_a, g_a_variables = self.__getActiveConstraints()

        KM.Logger.PrintInfo("ShapeOpt", "Assemble vector of objective gradient.")
        nabla_f = KM.Vector()
        s = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f, KSO.DF1DX_MAPPED)

        if len(g_a) == 0:
            KM.Logger.PrintInfo("ShapeOpt", "No constraints active, use negative objective gradient as search direction.")
            s = nabla_f * (-1.0)
            s *= self.step_size / s.norm_inf()
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
            if c.norm_inf() <= self.max_correction_share * self.step_size:
                delta = self.step_size - c.norm_inf()
                s *= delta/s.norm_inf()
            else:
                KM.Logger.PrintWarning("ShapeOpt", f"Correction is scaled down from {c.norm_inf()} to {self.max_correction_share * self.step_size}.")
                c *= self.max_correction_share * self.step_size / c.norm_inf()
                s *= (1.0 - self.max_correction_share) * self.step_size / s.norm_inf()
        else:
            s *= self.step_size / s.norm_inf()

        self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.SEARCH_DIRECTION)
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, c, KSO.CORRECTION)
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, s+c, KSO.CONTROL_POINT_UPDATE)

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
        additional_values_to_log["inf_norm_s"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.SEARCH_DIRECTION)
        additional_values_to_log["inf_norm_c"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.CORRECTION)
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
