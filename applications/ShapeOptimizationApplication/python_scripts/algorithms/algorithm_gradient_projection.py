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
from KratosMultiphysics.ShapeOptimizationApplication.loggers import data_logger
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_variable_utilities import WriteDictionaryDataOnNodalVariable
from KratosMultiphysics.ShapeOptimizationApplication import OptimizationUtilities as optimization_utilities

# ==============================================================================
class AlgorithmGradientProjection(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        super().__init__(optimization_settings, analyzer, communicator, model_part_controller)

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

        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

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

        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.SEARCH_DIRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.CORRECTION)

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Gradient projection algorithm only supports one objective function!")
        if self.constraints.size() == 0:
            raise RuntimeError("Gradient projection algorithm requires definition of at least one constraint!")

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

            self._initializeNewShape()

            self.__analyzeShape()

            self.__computeShapeUpdate()

            self.__logCurrentOptimizationStep()

            KM.Logger.Print("")
            KM.Logger.PrintInfo("ShapeOpt", "Time needed for current optimization step = ", timer.GetLapTime(), "s")
            KM.Logger.PrintInfo("ShapeOpt", "Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

            if self.__isAlgorithmConverged():
                break
            else:
                self.__determineAbsoluteChanges()

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

        self.model_part_controller.DampNodalVariableIfSpecified(KSO.DF1DX)

        # project and damp constraint gradients
        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            conGradientDict = self.communicator.getStandardizedGradient(con_id)
            gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
            WriteDictionaryDataOnNodalVariable(conGradientDict, self.optimization_model_part, gradient_variable)

            if constraint["project_gradient_on_surface_normals"].GetBool():
                self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(gradient_variable)

            self.model_part_controller.DampNodalVariableIfSpecified(gradient_variable)

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):
        for mapper in self.mappers.values():
            mapper.Update()
            mapper.InverseMap(KSO.DF1DX, KSO.DF1DX_MAPPED)

        for constraint in self.constraints:
            con_id = constraint["identifier"].GetString()
            gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
            mapped_gradient_variable = self.constraint_gradient_variables[con_id]["mapped_gradient"]
            mapper.InverseMap(gradient_variable, mapped_gradient_variable)

        self.__computeControlPointUpdate()

        for mapper in self.mappers.values():
            mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
        
        self.model_part_controller.DampNodalVariableIfSpecified(KSO.SHAPE_UPDATE)

    # --------------------------------------------------------------------------
    def __computeControlPointUpdate(self):
        """adapted from https://msulaiman.org/onewebmedia/GradProj_2.pdf"""
        g_a, g_a_variables = self.__getActiveConstraints()
        if len(g_a) == 0:
            KM.Logger.PrintInfo("ShapeOpt", "No constraints active, use negative objective gradient as search direction.")
            for design_surface in self.design_surfaces.values():
                s = KM.Vector()
                s = nabla_f * (-1.0)
                s *= self.step_size / s.norm_inf()
                optimization_utilities.AssignVectorToVariable(design_surface, s, KSO.SEARCH_DIRECTION)
                optimization_utilities.AssignVectorToVariable(design_surface, [0.0]*len(s), KSO.CORRECTION)
                optimization_utilities.AssignVectorToVariable(design_surface, s, KSO.CONTROL_POINT_UPDATE)
            return

        nabla_f = KM.Vector()
        for design_surface in self.design_surfaces.values():
            KM.Logger.PrintInfo("ShapeOpt", "Assemble vector of objective gradient.")
            s = KM.Vector()
            optimization_utilities.AssembleVector(design_surface, nabla_f, KSO.DF1DX_MAPPED, True)

            KM.Logger.PrintInfo("ShapeOpt", "Assemble matrix of constraint gradient.")
            N = KM.Matrix()
            optimization_utilities.AssembleMatrix(design_surface, N, g_a_variables, True)  # TODO check if gradients are 0.0! - in cpp

        settings = KM.Parameters('{ "solver_type" : "LinearSolversApplication.dense_col_piv_householder_qr" }')
        solver = dense_linear_solver_factory.ConstructSolver(settings)

        KM.Logger.PrintInfo("ShapeOpt", "Calculate projected search direction and correction.")
        c = KM.Vector()
        optimization_utilities.CalculateProjectedSearchDirectionAndCorrection(
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

        previous_end = 0
        for design_surface in self.design_surfaces.values():
            initial_index = previous_end
            final_index = initial_index + design_surface.NumberOfNodes()*3
            previous_end = final_index

            optimization_utilities.AssignVectorToVariable(design_surface, s[initial_index, final_index], KSO.SEARCH_DIRECTION)
            optimization_utilities.AssignVectorToVariable(design_surface, c[initial_index, final_index], KSO.CORRECTION)
            optimization_utilities.AssignVectorToVariable(design_surface, (s+c)[initial_index, final_index], KSO.CONTROL_POINT_UPDATE)

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
        if constraint["type"].GetString() == "=":
            return True
        elif constraint_value >= 0:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        inf_norm_s = 0.0
        inf_norm_c = 0.0
        for name, design_surface in self.design_surfaces.items():
            inf_norm_s = inf_norm_s + optimization_utilities.ComputeMaxNormOfNodalVariable(design_surface, KSO.SEARCH_DIRECTION)
            inf_norm_c = inf_norm_c + optimization_utilities.ComputeMaxNormOfNodalVariable(design_surface, KSO.CORRECTION)

        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.step_size
        additional_values_to_log["inf_norm_s"] = inf_norm_s
        additional_values_to_log["inf_norm_c"] = inf_norm_c
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
        for name, design_surface in self.design_surfaces.items():
            optimization_utilities.AddFirstVariableToSecondVariable(design_surface, KSO.CONTROL_POINT_UPDATE, KSO.CONTROL_POINT_CHANGE)
            optimization_utilities.AddFirstVariableToSecondVariable(design_surface, KSO.SHAPE_UPDATE, KSO.SHAPE_CHANGE)

# ==============================================================================
