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

            self.__initializeNewShape()

            self.__analyzeShape()

            self.__computeShapeUpdate()

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
    def __computeSensitivityHeatmap(self):

        def ___ComputeHeatmap(norm_type, sens_type):

            # relax_coeff = 0.8
            # reciprocal relaxation
            relax_coeff = 1 / self.optimization_iteration
            heat = []

            # read objective gradient
            # if sens_type == "raw":
            #     df_dx = ReadNodalVariableToList(self.design_surface, KSO.DF1DX)
            # elif sens_type == "mapped":
            df_dx = ReadNodalVariableToList(self.design_surface, KSO.DF1DX_MAPPED)

            # DF1DX individual heatmap
            heatmap_dfdx_name = "HEATMAP_DF1DX"
            if self.optimization_iteration == 1:
                heat_dfdx_relaxed = df_dx
            else:
                prev_heat_dfdx = KM.Vector()
                self.optimization_utilities.AssembleVector(self.design_surface, prev_heat_dfdx, KM.KratosGlobals.GetVariable(heatmap_dfdx_name))
                heat_dfdx_relaxed = []
                for i in range(len(self.design_surface.Nodes)):
                    for dim in range(3):
                        heat_dfdx_relaxed.append(relax_coeff * df_dx[3*i+dim] + (1 - relax_coeff) * prev_heat_dfdx[3*i+dim])

            WriteListToNodalVariable(heat_dfdx_relaxed, self.design_surface, KM.KratosGlobals.GetVariable(heatmap_dfdx_name))

            # normalize objective gradient
            if norm_type == "MAX":
                df_dx_norm = cm.NormInf3D(df_dx)
            elif norm_type == "L2":
                df_dx_norm = cm.Norm2(df_dx)
            elif norm_type == "VALUE":
                f = self.communicator.getStandardizedValue(self.objectives[0]["identifier"].GetString())
                df_dx_norm = abs(f)

            if df_dx_norm != 0.0:
                df_dx_normalized = cm.ScalarVectorProduct(1/df_dx_norm, df_dx)
            else:
                df_dx_normalized = [0] * len(df_dx)

            dc_dx_normalized = {}
            for itr, constraint in enumerate(self.constraints):
                # read constraint gradients
                con_id = constraint["identifier"].GetString()
                # if sens_type == "raw":
                #     gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
                # elif sens_type == "mapped":
                gradient_variable = self.constraint_gradient_variables[con_id]["mapped_gradient"]

                dci_dx = ReadNodalVariableToList(self.design_surface, gradient_variable)

                # normalize constraints
                if norm_type == "MAX":
                    dci_dx_norm = cm.NormInf3D(dci_dx)
                elif norm_type == "L2":
                    dci_dx_norm = cm.Norm2(dci_dx)
                elif norm_type == "VALUE":
                    c = self.communicator.getStandardizedValue(con_id)
                    dci_dx_norm = abs(c)

                if dci_dx_norm != 0.0:
                    dc_dx_normalized.update({con_id : cm.ScalarVectorProduct(1/dci_dx_norm, dci_dx)})
                else:
                    dc_dx_normalized.update({con_id : [0] * len(dci_dx)})

                # DCiDX individual heatmap
                heatmap_dcidx_name = "HEATMAP_" + "DC" + str(itr+1) + "DX"
                if self.optimization_iteration == 1:
                    heat_dcidx_relaxed = dci_dx
                else:
                    prev_heat_dcidx = KM.Vector()
                    self.optimization_utilities.AssembleVector(self.design_surface, prev_heat_dcidx, KM.KratosGlobals.GetVariable(heatmap_dcidx_name))
                    heat_dcidx_relaxed = []
                    for i in range(len(self.design_surface.Nodes)):
                        for dim in range(3):
                            heat_dcidx_relaxed.append(relax_coeff * dci_dx[3*i+dim] + (1 - relax_coeff) * prev_heat_dcidx[3*i+dim])

                WriteListToNodalVariable(heat_dcidx_relaxed, self.design_surface, KM.KratosGlobals.GetVariable(heatmap_dcidx_name))

            # fill heat map for each node
            for i in range(len(self.design_surface.Nodes)):
                df_dx_i = df_dx_normalized[3*i:3*i+3]
                df_dx_i_norm = cm.Norm2(df_dx_i)

                heat_i = df_dx_i_norm
                for dc_dx in dc_dx_normalized.values():
                    dc_dx_i = dc_dx[3*i:3*i+3]
                    dc_dx_i_norm = cm.Norm2(dc_dx_i)
                    heat_i = max(heat_i, dc_dx_i_norm)

                heat.append(heat_i)

            heat_map_name = "HEATMAP" + "_" + norm_type
            # if sens_type == "mapped":
            #     heat_map_name += "_" + "MAPPED"

            # Heatmap Relaxed
            if self.optimization_iteration == 1:
                heat_relaxed = heat
            else:
                prev_heat = KM.Vector()
                self.optimization_utilities.AssembleScalar(self.design_surface, prev_heat, KM.KratosGlobals.GetVariable(heat_map_name + "_RELAXED"))
                heat_relaxed = []
                # self.optimization_utilities.AssembleScalar(self.design_surface, heat_relaxed, KM.KratosGlobals.GetVariable(heat_map_name + "_RELAXED"))
                for i in range(len(self.design_surface.Nodes)):
                    heat_relaxed.append(relax_coeff * heat[i] + (1 - relax_coeff) * prev_heat[i])

            # self.optimization_utilities.AssignScalarToVariable(self.design_surface, heat_relaxed, KM.KratosGlobals.GetVariable(heat_map_name + "_RELAXED"))
            WriteListToNodalVariable(heat_relaxed, self.design_surface, KM.KratosGlobals.GetVariable(heat_map_name + "_RELAXED"), 1)
            WriteListToNodalVariable(heat, self.design_surface, KM.KratosGlobals.GetVariable(heat_map_name), 1)

        # ___ComputeHeatmap(norm_type="MAX", sens_type="raw")
        ___ComputeHeatmap(norm_type="MAX", sens_type="mapped")

        # ___ComputeHeatmap(norm_type="L2", sens_type="raw")
        ___ComputeHeatmap(norm_type="L2", sens_type="mapped")

        # ___ComputeHeatmap(norm_type="VALUE", sens_type="raw")
        ___ComputeHeatmap(norm_type="VALUE", sens_type="mapped")

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
