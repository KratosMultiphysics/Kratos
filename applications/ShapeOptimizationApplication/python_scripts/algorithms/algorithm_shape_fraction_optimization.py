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

from sympy import false

# Kratos Core and Apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics.LinearSolversApplication import dense_linear_solver_factory
import KratosMultiphysics.eigen_solver_factory as eigen_solver_factory

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
        default_algorithm_settings = KM.Parameters("""
        {
            "name"                    : "shape_fraction_optimization",
            "max_correction_share"    : 0.75,
            "max_iterations"          : 100,
            "relative_tolerance"      : 1e-3,
            "shape_fraction" : {
                "max_fraction"  : 0.5,
                "tolerance"     : 0.1,
                "inner_tolerance": 0.01,
                "max_inner_steps": 10
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
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.DP1DX)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.DPF1DX)

        self.max_shape_fraction = self.algorithm_settings["shape_fraction"]["max_fraction"].GetDouble()
        self.frac_tolerance = self.algorithm_settings["shape_fraction"]["tolerance"].GetDouble()
        self.inner_tolerance = self.algorithm_settings["shape_fraction"]["inner_tolerance"].GetDouble()
        self.max_inner_steps = self.algorithm_settings["shape_fraction"]["max_inner_steps"].GetDouble()
        self.inner_step = 0
        self.penalty_factor = 0.0
        self.previous_objective_value = None

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Gradient projection algorithm only supports one objective function!")
        # if self.constraints.size() == 0:
        #     raise RuntimeError("Gradient projection algorithm requires definition of at least one constraint!")

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

        # Compute pseudo objective
        objective_value = self.communicator.getStandardizedValue(self.objectives[0]["identifier"].GetString())

        objGradientDict = self.communicator.getStandardizedGradient(self.objectives[0]["identifier"].GetString())
        WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, KSO.DF1DX)
        objective_gradient = KM.Vector()
        self.optimization_utilities.AssembleVector(self.optimization_model_part, objective_gradient, KSO.DF1DX)
        objective_gradient_norm = objective_gradient.norm_2()
        KM.Logger.Print("objective_gradient_norm: {}".format(objective_gradient_norm))

        self.penalty_value, penaltyGradientDict = self.__computePenalty()

        WriteDictionaryDataOnNodalVariable(penaltyGradientDict, self.optimization_model_part, KSO.DP1DX)
        penalty_gradient = KM.Vector()
        self.optimization_utilities.AssembleVector(self.optimization_model_part, penalty_gradient, KSO.DP1DX)
        penalty_gradient_norm = penalty_gradient.norm_2()
        KM.Logger.Print("penalty_gradient_norm: {}".format(penalty_gradient_norm))

        KM.Logger.PrintInfo("objective_value: {}".format(objective_value))
        KM.Logger.PrintInfo("self.penalty_value: {}".format(self.penalty_value))

        # use penalty factor from start
        if self.penalty_factor == 0.0 and penalty_gradient_norm > 0:
            self.penalty_factor = 0.1 * objective_gradient_norm / penalty_gradient_norm

        if penalty_gradient_norm > 0:
            self.inner_step += 1
        # use pseudo objective for inner convergence check
        pseudo_objective_value = objective_value + self.penalty_factor * max(0, self.penalty_value)**2
        inner_converged = self.__checkInnerConvergence(pseudo_objective_value)

        # use only objective for inner convergence check
        # inner_converged = self.__checkInnerConvergence(objective_value)

        if inner_converged:
            KM.Logger.PrintInfo("ShapeFractionOptimization", "Updating penalty factor:")
            self.__updatePenaltyFactor(objective_gradient_norm, penalty_gradient_norm)
            KM.Logger.PrintInfo("self.penalty_factor: {}".format(self.penalty_factor))

        pseudo_objective_value = objective_value + self.penalty_factor * max(0, self.penalty_value)**2
        KM.Logger.PrintInfo("pseudo_objective_value: {}".format(pseudo_objective_value))

        pseudo_objective_gradient = objective_gradient + self.penalty_factor * penalty_gradient
        self.optimization_utilities.AssignVectorToVariable(self.optimization_model_part, pseudo_objective_gradient, KSO.DPF1DX)

        pseudo_objective_gradient_norm = pseudo_objective_gradient.norm_2()
        KM.Logger.Print("pseudo_objective_gradient_norm: {}".format(pseudo_objective_gradient_norm))

        # project and damp objective gradients
        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(KSO.DPF1DX)

        self.model_part_controller.DampNodalSensitivityVariableIfSpecified(KSO.DPF1DX)

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
    def __computePenalty(self):

        def __calculateNodalValue(self, node):
            shape_change = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
            norm = cm.Norm2(shape_change)

            if norm > self.frac_tolerance:
                return 1.0
            else:
                return 0.0

        def __calculateNodalGradient(self, node, quantile, neglect_feasible_nodes):
            shape_change = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
            norm = cm.Norm2(shape_change)

            if norm > self.frac_tolerance:
                if neglect_feasible_nodes and norm > quantile:
                    return [0.0, 0.0, 0.0]
                else:
                    return [
                        (shape_change[0] / norm)*(quantile / norm)**1,
                        (shape_change[1] / norm)*(quantile / norm)**1,
                        (shape_change[2] / norm)*(quantile / norm)**1
                    ]
            else:
                return [0.0, 0.0, 0.0]

        def __getQuantile(self):

            nodal_variable = KM.KratosGlobals.GetVariable("SHAPE_CHANGE")
            shape_change = ReadNodalVariableToList(self.design_surface, nodal_variable)

            shape_change_norm = []

            for i in range(int(len(shape_change)/3)):
                shape_change_norm.append(cm.Norm2(shape_change[3*i:3*i+2]))

            shape_change_norm.sort()

            index = round(len(shape_change_norm)*(1-self.max_shape_fraction))

            quantile = shape_change_norm[index]

            KM.Logger.PrintInfo("ShapeChange Quantile {}: {}".format((1-self.max_shape_fraction), quantile))

            return quantile

        KM.Logger.PrintInfo("ShapeFractionOptimization", "Starting calculation of penalty value and gradient:")

        startTime = timer.time()
        penalty_value = 0.0
        penalty_gradient = {}
        model_part = self.design_surface
        total_nodes = len(model_part.Nodes)
        g = 0.0
        quantile = 1.0
        # quantile = __getQuantile(self)
        for node in model_part.Nodes:
            g_i = __calculateNodalValue(self, node)
            # gradient_i = __calculateNodalGradient(self, node, quantile)
            # penalty_gradient[node.Id] = gradient_i
            g += g_i

        penalty_value = g / total_nodes - self.max_shape_fraction
        KM.Logger.PrintInfo("ShapeFractionOptimization," "Current shape fraction = {}".format(g / total_nodes))
        # penalty_value = max(0, penalty_value)

        neglect_feasible_nodes = False
        for node in model_part.Nodes:
            if penalty_value > 0.0:
                gradient_i = __calculateNodalGradient(self, node, quantile, neglect_feasible_nodes)
            else:
                gradient_i = [0.0, 0.0, 0.0]
            penalty_gradient[node.Id] = gradient_i

        KM.Logger.PrintInfo("ShapeFractionOptimization", "Time needed for calculating the penalty value and gradient = ",round(timer.time() - startTime,2),"s")

        return penalty_value, penalty_gradient

    # --------------------------------------------------------------------------
    def __checkInnerConvergence(self, objective_value):

        # tolerance = 1e-3
        tolerance = self.inner_tolerance

        is_converged = False
        if self.previous_objective_value == None:
            self.reference_objective_value = objective_value
            self.previous_objective_value = objective_value
            self.relative_change = 0
        elif self.inner_step >= self.max_inner_steps:
            self.previous_objective_value = objective_value
            is_converged = True
        else:
            # if self.previous_objective_value != 0:
            self.relative_change = (objective_value - self.previous_objective_value) / self.reference_objective_value
            KM.Logger.PrintInfo("ShapeFractionOptimization", "relative_change: {}".format(self.relative_change))
            if abs(self.relative_change) < tolerance:
                is_converged = True
            # else:
            #     if objective_value < tolerance:
            #         is_converged = True

            self.previous_objective_value = objective_value

        if is_converged:
            self.max_inner_steps = 0

        return is_converged

    # --------------------------------------------------------------------------
    def __updatePenaltyFactor(self, objective_gradient_norm, penalty_gradient_norm):

        if penalty_gradient_norm > 0.0:
            if self.penalty_factor == 0.0:
                # initialize penalty factor such that the penalty gradient norm is 10% of the objective gradient norm
                self.penalty_factor = 0.1 * objective_gradient_norm / penalty_gradient_norm
            else:
                self.penalty_factor *= 1.25 # increase penalty factor by 25%
        # not necessary to start again
        # else:
        #     self.penalty_factor = 0.0

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):
        self.mapper.Update()
        self.mapper.InverseMap(KSO.DPF1DX, KSO.DF1DX_MAPPED)

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
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f, KSO.DF1DX_MAPPED)

        "active subspace strategy"
        model_part = self.design_surface
        total_nodes = len(model_part.Nodes)
        number_of_active_nodes = int(total_nodes * self.max_shape_fraction)

        # TODO: fix compilierung mit DUSE_EIGEN_FEAST=ON
        # eigen_settings = KM.Parameters(
        #     '''{
        #             "solver_type" : "feast",
        #             "symmetric": true,
        #             "search_highest_eigenvalues": true,
        #             "number_of_eigenvalues": 0
        #         }''')
        eigen_settings = KM.Parameters(
            '''{
                    "solver_type" : "spectra_sym_g_eigs_shift",
                    "number_of_eigenvalues": 1,
                    "shift": 0.0
                }''')
        eigen_settings["number_of_eigenvalues"].SetDouble(total_nodes - number_of_active_nodes)
        KM.Logger.PrintInfo("ShapeOpt", "Construct solver.")
        eigen_solver = eigen_solver_factory.ConstructSolver(eigen_settings)

        active_nodes = KM.Vector()
        startTime = timer.time()
        self.optimization_utilities.ConstructActiveSubspace(nabla_f, active_nodes, eigen_solver)
        KM.Logger.PrintInfo("ShapeFractionOptimization", "Time needed for calculating the active subspace = ",round(timer.time() - startTime,2),"s")

        s = KM.Vector()

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
        additional_values_to_log["penalty_value"] = self.penalty_value
        additional_values_to_log["penalty_factor"] = self.penalty_factor
        additional_values_to_log["df_rel_p"] = self.relative_change
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
