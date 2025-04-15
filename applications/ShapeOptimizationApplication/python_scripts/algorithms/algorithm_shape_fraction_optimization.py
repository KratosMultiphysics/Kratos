# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Schmölz David, https://github.com/dschmoelz
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
import KratosMultiphysics as Kratos
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.LinearSolversApplication import dense_linear_solver_factory

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication.utilities import custom_math as cm
from KratosMultiphysics.ShapeOptimizationApplication.algorithms.algorithm_base import OptimizationAlgorithm
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
from KratosMultiphysics.ShapeOptimizationApplication.loggers import data_logger_factory
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_variable_utilities import WriteDictionaryDataOnNodalVariable, ReadNodalVariableToList

import math
import time as timer

# ==============================================================================
class AlgorithmShapeFractionOptimization(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = Kratos.Parameters("""
        {
            "name"                    : "shape_fraction_optimization",
            "max_correction_share"    : 0.75,
            "max_iterations"          : 100,
            "relative_tolerance"      : 1e-3,
            "shape_fraction" : {
                "penalty_method"      : "exterior",
                "max_fraction"        : 0.5,
                "nodal_tolerance"     : 0.1,
                "inner_tolerance"     : 0.01,
                "max_inner_steps"     : 10,
                "initial_penalty_factor": 0.1,
                "penalty_scale_factor": 1.25
            },
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
        for itr, constraint in enumerate(self.constraints.values()):
            self.constraint_gradient_variables.update({
                constraint["identifier"].GetString() : {
                    "gradient": Kratos.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX"),
                    "mapped_gradient": Kratos.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")
                }
            })
        self.max_correction_share = self.algorithm_settings["max_correction_share"].GetDouble()

        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.SEARCH_DIRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.CORRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.DP1DX)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.DPF1DX)

        # shape fraction related settings
        self.penalty_method = self.algorithm_settings["shape_fraction"]["penalty_method"].GetString()
        self.max_shape_fraction = self.algorithm_settings["shape_fraction"]["max_fraction"].GetDouble()
        self.frac_tolerance = self.algorithm_settings["shape_fraction"]["nodal_tolerance"].GetDouble()
        self.inner_tolerance = self.algorithm_settings["shape_fraction"]["inner_tolerance"].GetDouble()
        self.max_inner_steps = self.algorithm_settings["shape_fraction"]["max_inner_steps"].GetDouble()
        self.initial_penalty_factor = self.algorithm_settings["shape_fraction"]["initial_penalty_factor"].GetDouble()
        self.gamma = self.algorithm_settings["shape_fraction"]["penalty_scale_factor"].GetDouble()
        if self.gamma <= 1.0:
            raise RuntimeError("Shape fraction algorithm: 'penalty_scale_factor' has to be larger than 1.0!")
        self.inner_step = 0
        self.penalty_factor = 0.0
        self.epsilon = -0.2
        self.previous_objective_value = None

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Shape fraction algorithm only supports one objective function!")

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
            Kratos.Logger.Print("")
            Kratos.Logger.Print("===============================================================================")
            Kratos.Logger.PrintInfo("ShapeOpt", timer.GetTimeStamp(), ": Starting optimization iteration ", self.optimization_iteration)
            Kratos.Logger.Print("===============================================================================\n")

            timer.StartNewLap()

            self.__initializeNewShape()

            self.__analyzeShape()

            self.__computeShapeUpdate()

            self.__logCurrentOptimizationStep()

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
    def __initializeNewShape(self):
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)
        self.model_part_controller.UpdateMeshAccordingInputVariable(KSO.SHAPE_UPDATE)
        self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyzeShape(self):
        self.communicator.initializeCommunication()
        self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
        self.communicator.requestGradientOf(self.objectives[0]["identifier"].GetString())

        for constraint in self.constraints.values():
            con_id =  constraint["identifier"].GetString()
            self.communicator.requestValueOf(con_id)
            self.communicator.requestGradientOf(con_id)

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.optimization_iteration, self.communicator)

        # compute normals only if required
        surface_normals_required = self.objectives[0]["project_gradient_on_surface_normals"].GetBool()
        for constraint in self.constraints.values():
            if constraint["project_gradient_on_surface_normals"].GetBool():
                surface_normals_required = True

        if surface_normals_required:
            self.model_part_controller.ComputeUnitSurfaceNormals()

        # objective value
        objective_value = self.communicator.getStandardizedValue(self.objectives[0]["identifier"].GetString())

        # objective gradient
        objGradientDict = self.communicator.getStandardizedGradient(self.objectives[0]["identifier"].GetString())
        WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, KSO.DF1DX)
        objective_gradient = Kratos.Vector()
        self.optimization_utilities.AssembleVector(self.optimization_model_part, objective_gradient, KSO.DF1DX)
        objective_gradient_norm = objective_gradient.norm_2()

        # response value g of penalty method
        response_value = self.__computeResponseValue()
        Kratos.Logger.PrintInfo("ShapeFractionOptimization", f"Shape Fraction value = {response_value}")
        # penalty value
        # exterior:             p = max(0, g)**2
        # extended interior:    p = -1/g                            if g =< epsilon (interior penalty)
        #                       p = - 2*epsilon - g / epsilon**2    if g > epsilon  (exterior penalty)
        self.penalty_value = self.__computePenaltyValue(response_value)

        # set up initial penalty factor
        if self.penalty_factor == 0.0:
            self.__SetUpInitialPenaltyFactor(objective_value, objective_gradient_norm)
        self.__IncrementInnerStep()

        # pseudo objective value
        pseudo_objective_value = self.__computePseudoObjectiveValue(objective_value)

        # check inner convergence
        inner_converged = self.__checkInnerConvergence(pseudo_objective_value)
        # update penalty factor if inner loop converged
        if inner_converged:
            Kratos.Logger.PrintInfo("ShapeFractionOptimization", "Updating penalty factor.")
            self.__updatePenaltyFactor()

        # update pseudo objective value
        self.pseudo_objective_value = self.__computePseudoObjectiveValue(objective_value)

        # compute penalty gradient
        # exterior:             p = 2 * g * dg/dx
        # extended interior:    p = -1/g                            if g =< epsilon (interior penalty)
        #                       p = - 2*epsilon - g / epsilon**2    if g > epsilon  (exterior penalty)
        penalty_gradient = self.__computePenaltyGradient(response_value)
        self.optimization_utilities.AssignVectorToVariable(self.optimization_model_part, penalty_gradient, KSO.DP1DX)

        # pseudo objective gradient
        pseudo_objective_gradient = objective_gradient + self.penalty_factor * penalty_gradient
        self.optimization_utilities.AssignVectorToVariable(self.optimization_model_part, pseudo_objective_gradient, KSO.DPF1DX)

        # project and damp objective gradients
        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(KSO.DPF1DX)

        self.model_part_controller.DampNodalSensitivityVariableIfSpecified(KSO.DPF1DX)

        # project and damp constraint gradients
        for constraint in self.constraints.values():
            con_id = constraint["identifier"].GetString()
            conGradientDict = self.communicator.getStandardizedGradient(con_id)
            gradient_variable = self.constraint_gradient_variables[con_id]["gradient"]
            WriteDictionaryDataOnNodalVariable(conGradientDict, self.optimization_model_part, gradient_variable)

            if constraint["project_gradient_on_surface_normals"].GetBool():
                self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(gradient_variable)

            self.model_part_controller.DampNodalSensitivityVariableIfSpecified(gradient_variable)

    # --------------------------------------------------------------------------
    def __computePenaltyValue(self, response_value):

        if self.penalty_method == "extended_interior":
            g = response_value
            if g > self.epsilon:
                penalty_value = - (2*self.epsilon - g) / self.epsilon**2
            else:
                penalty_value = - 1 / g

        elif self.penalty_method == "exterior":
            penalty_value = response_value

        return penalty_value

    # --------------------------------------------------------------------------
    def __computeResponseValue(self):

        def __calculateNodalValueExterior(self, node, quantile):
            shape_change = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
            norm = cm.Norm2(shape_change)

            if norm > quantile:
                return 0.0

            if norm > self.frac_tolerance:
                return norm
            else:
                return 0.0

        def __calculateNodalValueExtendedInterior(self, node, quantile):
            shape_change = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
            norm = cm.Norm2(shape_change)

            if norm > quantile:
                return 0.0
            else:
                return norm

        def __calculateIntegralToleranceOfNode(self, node, quantile):
            shape_change = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
            norm = cm.Norm2(shape_change)

            if norm > quantile:
                return 0.0
            else:
                return norm

        Kratos.Logger.PrintInfo("ShapeFractionOptimization", "Starting calculation of penalty response value:")

        startTime = timer.time()
        total_integral = 0.0
        integral_tol = 0.0
        quantile = self.__getQuantile()
        for node in self.design_surface.Nodes:

            if self.penalty_method == "exterior":
                g_i = __calculateNodalValueExterior(self, node, quantile)
            elif self.penalty_method == "extended_interior":
                g_i = __calculateNodalValueExtendedInterior(self, node, quantile)
                g_tol_i = __calculateIntegralToleranceOfNode(self, node, quantile)
                integral_tol += g_tol_i

            total_integral += g_i

        if self.penalty_method == "extended_interior":
            g = total_integral - integral_tol

        elif self.penalty_method == "exterior":
            g = total_integral

        Kratos.Logger.PrintInfo("ShapeFractionOptimization", "Time needed for calculating the penalty response value = ",round(timer.time() - startTime,2),"s")

        return g


    # --------------------------------------------------------------------------
    def __computePenaltyGradient(self, penalty_value):

        def __calculateNodalGradient(self, node, quantile):
            shape_change = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
            norm = cm.Norm2(shape_change)

            if norm > self.frac_tolerance:
                if norm > quantile:
                    return [0.0, 0.0, 0.0]
                else:
                    return [
                        shape_change[0] / norm,
                        shape_change[1] / norm,
                        shape_change[2] / norm
                    ]
            else:
                return [0.0, 0.0, 0.0]

        Kratos.Logger.PrintInfo("ShapeFractionOptimization", "Starting calculation of penalty gradient:")

        startTime = timer.time()
        penalty_gradient = Kratos.Vector(3*len(self.design_surface.Nodes))
        quantile = self.__getQuantile()

        for i, node in enumerate(self.design_surface.Nodes):
            if self.penalty_method == "exterior":
                if penalty_value > 0.0:
                    gradient_i = __calculateNodalGradient(self, node, quantile)
                else:
                    gradient_i = [0.0, 0.0, 0.0]

            elif self.penalty_method == "extended_interior":
                gradient_i = __calculateNodalGradient(self, node, quantile)
                if penalty_value > self.epsilon:
                    gradient_i[0] /= self.epsilon**2
                    gradient_i[1] /= self.epsilon**2
                    gradient_i[2] /= self.epsilon**2
                else:
                    gradient_i[0] /= penalty_value**2
                    gradient_i[1] /= penalty_value**2
                    gradient_i[2] /= penalty_value**2

            penalty_gradient[3*i:3*i+3] = gradient_i

        Kratos.Logger.PrintInfo("ShapeFractionOptimization", "Time needed for calculating the penalty gradient = ",round(timer.time() - startTime,2),"s")

        return penalty_gradient

    def __getQuantile(self):

        nodal_variable = Kratos.KratosGlobals.GetVariable("SHAPE_CHANGE")
        shape_change = ReadNodalVariableToList(self.design_surface, nodal_variable)

        shape_change_norm = []

        for i in range(int(len(shape_change)/3)):
            shape_change_norm.append(cm.Norm2(shape_change[3*i:3*i+3]))

        shape_change_norm.sort()

        index = round(len(shape_change_norm)*(1-self.max_shape_fraction))

        quantile = shape_change_norm[index]

        Kratos.Logger.PrintInfo(f"ShapeFractionOptimization: Shape change quantile {(1-self.max_shape_fraction)} = {quantile}")

        return quantile

    def __computeShapeFraction(self):

        n_active = 0
        for node in self.design_surface.Nodes:
            shape_change = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
            norm = cm.Norm2(shape_change)

            if norm > self.frac_tolerance:
                n_active += 1

        total_nodes = len(self.design_surface.Nodes)

        shape_fraction = n_active / total_nodes

        return shape_fraction

    def __computePseudoObjectiveValue(self, objective_value):
        if self.penalty_method == "exterior":
            pseudo_objective_value = objective_value + self.penalty_factor * max(0, self.penalty_value)**2
        elif self.penalty_method == "extended_interior":
            pseudo_objective_value = objective_value + self.penalty_factor * self.penalty_value

        return pseudo_objective_value

    # --------------------------------------------------------------------------
    def __IncrementInnerStep(self):
        if self.penalty_method == "exterior" and self.penalty_value > 0:
            self.inner_step += 1
        elif self.penalty_method == "extended_interior":
            self.inner_step += 1

    # --------------------------------------------------------------------------
    def __checkInnerConvergence(self, objective_value):

        tolerance = self.inner_tolerance

        is_converged = False
        if self.previous_objective_value is None:
            self.reference_objective_value = objective_value
            self.previous_objective_value = objective_value
            self.relative_change = 0
        elif self.inner_step >= self.max_inner_steps:
            self.previous_objective_value = objective_value
            is_converged = True
        else:
            self.relative_change = (objective_value - self.previous_objective_value) / self.reference_objective_value
            if abs(self.relative_change) < tolerance:
                is_converged = True

            self.previous_objective_value = objective_value

        if is_converged:
            self.inner_step = 0
            self.reference_objective_value = objective_value

        return is_converged

    # --------------------------------------------------------------------------
    def __SetUpInitialPenaltyFactor(self, objective_value, objective_gradient_norm):
        if self.penalty_method == "exterior":
            # compute penalty gradients for start penalty factor
            penalty_gradient = self.__computePenaltyGradient(self.penalty_value)
            penalty_gradient_norm = penalty_gradient.norm_2()
            if penalty_gradient_norm > 0:
                self.penalty_factor = self.initial_penalty_factor * objective_gradient_norm / penalty_gradient_norm
        elif self.penalty_method == "extended_interior":
            self.penalty_factor = self.initial_penalty_factor * abs(objective_value / self.penalty_value)
            a = 0.5
            self.C = - self.epsilon / (self.penalty_factor**a)

    # --------------------------------------------------------------------------
    def __updatePenaltyFactor(self):

        if self.penalty_method == "exterior" and self.penalty_value > 0.0:
            self.penalty_factor *= self.gamma # increase penalty factor by 25%

        elif self.penalty_method == "extended_interior":
            self.penalty_factor *= self.gamma # increase penalty factor by 25%
            a = 0.5
            self.epsilon = - self.C * (self.penalty_factor**a)

        Kratos.Logger.PrintInfo("ShapeOpt", f"New penalty factor = {self.penalty_factor}")

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):
        self.mapper.Update()
        self.mapper.InverseMap(KSO.DPF1DX, KSO.DF1DX_MAPPED)

        for constraint in self.constraints.values():
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

        Kratos.Logger.PrintInfo("ShapeOpt", "Assemble vector of objective gradient.")
        nabla_f = Kratos.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f, KSO.DF1DX_MAPPED)

        s = Kratos.Vector()

        if len(g_a) == 0:
            Kratos.Logger.PrintInfo("ShapeOpt", "No constraints active, use negative objective gradient as search direction.")
            s = nabla_f * (-1.0)
            s *= self.step_size / s.norm_inf()
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.SEARCH_DIRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, [0.0]*len(s), KSO.CORRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.CONTROL_POINT_UPDATE)
            return


        Kratos.Logger.PrintInfo("ShapeOpt", "Assemble matrix of constraint gradient.")
        N = Kratos.Matrix()
        self.optimization_utilities.AssembleMatrix(self.design_surface, N, g_a_variables)

        settings = Kratos.Parameters('{ "solver_type" : "LinearSolversApplication.dense_col_piv_householder_qr" }')
        solver = dense_linear_solver_factory.ConstructSolver(settings)

        Kratos.Logger.PrintInfo("ShapeOpt", "Calculate projected search direction and correction.")
        c = Kratos.Vector()
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
                Kratos.Logger.PrintWarning("ShapeOpt", f"Correction is scaled down from {c.norm_inf()} to {self.max_correction_share * self.step_size}.")
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

        for constraint in self.constraints.values():
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
                Kratos.Logger.PrintWarning("ShapeOpt", f"Gradient for constraint {identifier} is 0.0 - will not be considered!")
                return False
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        additional_values_to_log = {}
        additional_values_to_log["penalty_value"] = self.penalty_value
        additional_values_to_log["penalty_factor"] = self.penalty_factor
        additional_values_to_log["f_p"] = self.pseudo_objective_value
        additional_values_to_log["df_rel_p"] = self.relative_change
        additional_values_to_log["shape_fraction"] = self.__computeShapeFraction()
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
                Kratos.Logger.Print("")
                Kratos.Logger.PrintInfo("ShapeOpt", "Maximal iterations of optimization problem reached!")
                return True

            # Check for relative tolerance
            relative_change_of_objective_value = self.data_logger.GetValues("rel_change_objective")[self.optimization_iteration]
            if abs(relative_change_of_objective_value) < self.relative_tolerance:
                Kratos.Logger.Print("")
                Kratos.Logger.PrintInfo("ShapeOpt", "Optimization problem converged within a relative objective tolerance of ",self.relative_tolerance,"%.")
                return True

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges(self):
        self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.CONTROL_POINT_UPDATE, KSO.CONTROL_POINT_CHANGE)
        self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.SHAPE_UPDATE, KSO.SHAPE_CHANGE)

# ==============================================================================
