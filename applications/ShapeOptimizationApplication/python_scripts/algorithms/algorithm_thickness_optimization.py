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
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_variable_utilities import WriteDictionaryDataOnNodalVariable

import math
import numpy as np
from collections import OrderedDict

# ==============================================================================
class AlgorithmThicknessOptimization(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = KM.Parameters("""
        {
            "name"                    : "thickness_optimization",
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

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.design_surface = None
        self.data_logger = None
        self.optimization_utilities = None

        self.objectives = optimization_settings["objectives"]
        self.objective_gradient_dict = None
        self.constraints = optimization_settings["constraints"]
        self.constraint_gradients_dict = {}
        for constraint in self.constraints:
            self.constraint_gradients_dict.update({
                constraint["identifier"].GetString() : 0.0
            })
        self.max_correction_share = self.algorithm_settings["max_correction_share"].GetDouble()

        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_SEARCH_DIRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CORRECTION)

        self.number_of_design_variables = None

        self.thickness_update_dict = {}

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

            if self.optimization_iteration == 1:
                self.__initializeFirstIteration()
            else:
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
    def __initializeFirstIteration(self):
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)
        self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __initializeNewThickness(self):
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)
        print(f"Properties: {self.optimization_model_part.Properties}")
        for property in self.optimization_model_part.Properties:
            if property.Id in self.thickness_update_dict:
                current_thickness = property.GetValue(KM.THICKNESS)
                new_thickness = current_thickness + self.thickness_update_dict[property.Id]
                property.SetValue(KM.THICKNESS, new_thickness)

        i = 0
        for condition in self.optimization_model_part.Conditions:
            i += 1
            if i < 11:
                print(f"Condition: {condition.Id} Thickness: {condition.Properties.GetValue(KM.THICKNESS)}")
            else:
                break

        i = 0
        for property in self.optimization_model_part.Properties:
            i += 1
            if i < 11:
                print(f"Property: {property.Id} Thickness: {property.GetValue(KM.THICKNESS)}")
            else:
                break

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

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.optimization_iteration, self.communicator)

        # project and damp objective gradients
        objGradientDict = self.communicator.getStandardizedThicknessGradient(self.objectives[0]["identifier"].GetString())
        self.objective_gradient_dict = OrderedDict(sorted(objGradientDict.items()))
        self.property_ids = list(self.objective_gradient_dict.keys())
        # for visualization on mesh
        objElementGradientDict = dict()
        self.__mapPropertyGradientToElement(self.objective_gradient_dict, objElementGradientDict)
        self.__mapElementGradientToNode(objElementGradientDict, KSO.DF1DT)

        self.model_part_controller.DampNodalSensitivityVariableIfSpecified(KSO.DF1DT)

        if not self.number_of_design_variables:
            self.number_of_design_variables = len(self.objective_gradient_dict.values())

        # project and damp constraint gradients
        for itr, constraint in enumerate(self.constraints):
            con_id = constraint["identifier"].GetString()
            conGradientDict = self.communicator.getStandardizedThicknessGradient(con_id)
            self.constraint_gradients_dict[con_id] = OrderedDict(sorted(conGradientDict.items()))
            # for visualization on mesh
            gradient_variable = KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT")
            conElementGradientDict = dict()
            self.__mapPropertyGradientToElement(conGradientDict, conElementGradientDict)
            self.__mapElementGradientToNode(conElementGradientDict, gradient_variable)

            self.model_part_controller.DampNodalSensitivityVariableIfSpecified(gradient_variable)

    # --------------------------------------------------------------------------
    def __computeThicknessUpdate(self):

        thickness_update = self.__computeControlThicknessUpdate()
        for itr, key in enumerate(self.objective_gradient_dict):
            self.thickness_update_dict[key] = thickness_update[itr]

        self.model_part_controller.DampNodalUpdateVariableIfSpecified(KSO.THICKNESS_UPDATE)

    # --------------------------------------------------------------------------
    def __computeControlThicknessUpdate(self):
        """adapted from https://msulaiman.org/onewebmedia/GradProj_2.pdf"""
        g_a, g_a_gradients = self.__getActiveConstraints()

        KM.Logger.PrintInfo("ShapeOpt", "Assemble vector of objective gradient.")
        nabla_f = KM.Vector(list(self.objective_gradient_dict.values()))
        s = KM.Vector(self.number_of_design_variables)
        print(f"GradientProjection:: nabla_f: {nabla_f}")
        for itr, gradient in enumerate(g_a_gradients):
            print(f"GradientProjection:: dg{itr}_dt: {gradient}")

        if len(g_a) == 0:
            KM.Logger.PrintInfo("ShapeOpt", "No constraints active, use negative objective gradient as search direction.")
            s = nabla_f * (-1.0)
            s *= self.step_size / s.norm_inf()

            # for visualization
            self.__mapDesignVariableVectorToNodalVariable(s, KSO.THICKNESS_SEARCH_DIRECTION)
            self.__mapDesignVariableVectorToNodalVariable([0.0]*len(s), KSO.THICKNESS_CORRECTION)
            self.__mapDesignVariableVectorToNodalVariable(s, KSO.THICKNESS_UPDATE)
            print(f"SteepestDescent:: search direction: {s}")
            return s


        KM.Logger.PrintInfo("ShapeOpt", "Assemble matrix of constraint gradient.")
        N = KM.Matrix()
        self.optimization_utilities.AssembleMatrixFromGradientVectors(self.design_surface, N, g_a_gradients)

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
        print(f"GradientProjection:: Vor Skalierung")
        print(f"GradientProjection:: search direction: {s}")
        print(f"GradientProjection:: correction: {c}")

        if c.norm_inf() != 0.0:
            if c.norm_inf() <= self.max_correction_share * self.step_size:
                delta = self.step_size - c.norm_inf()
                if s.norm_inf() > 0:
                    s *= delta/s.norm_inf()
            else:
                KM.Logger.PrintWarning("ShapeOpt", f"Correction is scaled down from {c.norm_inf()} to {self.max_correction_share * self.step_size}.")
                c *= self.max_correction_share * self.step_size / c.norm_inf()
                s *= (1.0 - self.max_correction_share) * self.step_size / s.norm_inf()
        else:
            if s.norm_inf() > 0:
                s *= self.step_size / s.norm_inf()

        print(f"GradientProjection:: search direction: {s}")
        print(f"GradientProjection:: correction: {c}")
        print(f"GradientProjection:: control update: {s+c}")
        # for visualization
        self.__mapDesignVariableVectorToNodalVariable(s, KSO.THICKNESS_SEARCH_DIRECTION)
        self.__mapDesignVariableVectorToNodalVariable(c, KSO.THICKNESS_CORRECTION)
        self.__mapDesignVariableVectorToNodalVariable(s+c, KSO.THICKNESS_UPDATE)

        return s+c

    # --------------------------------------------------------------------------
    def __getActiveConstraints(self):
        active_constraint_values = []
        active_constraint_gradients = []

        for constraint in self.constraints:
            if self.__isConstraintActive(constraint):
                identifier = constraint["identifier"].GetString()
                constraint_value = self.communicator.getStandardizedValue(identifier)
                active_constraint_values.append(constraint_value)
                active_constraint_gradients.append(
                    KM.Vector(list(self.constraint_gradients_dict[identifier].values())))

        return active_constraint_values, active_constraint_gradients

    # --------------------------------------------------------------------------
    def __isConstraintActive(self, constraint):
        identifier = constraint["identifier"].GetString()
        constraint_value = self.communicator.getStandardizedValue(identifier)
        if constraint["type"].GetString() == "=" or constraint_value >= 0:
            gradient = KM.Vector(list(self.constraint_gradients_dict[identifier].values()))
            gradient_norm = gradient.norm_inf()
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

        # TODO: move to c++
        for node in self.optimization_model_part.Nodes:
            # absolute_control_update = node.GetSolutionStepValue(KSO.THICKNESS_CONTROL_CHANGE)
            # absolute_control_update += node.GetSolutionStepValue(KSO.THICKNESS_CONTROL_UPDATE)
            # node.SetSolutionStepValue(KSO.THICKNESS_CONTROL_CHANGE, 0, absolute_control_update)

            absolute_update = node.GetSolutionStepValue(KSO.THICKNESS_CHANGE)
            absolute_update += node.GetSolutionStepValue(KSO.THICKNESS_UPDATE)
            node.SetSolutionStepValue(KSO.THICKNESS_CHANGE, 0, absolute_update)

            thickness = node.GetSolutionStepValue(KSO.THICKNESS)
            thickness += node.GetSolutionStepValue(KSO.THICKNESS_UPDATE)
            node.SetSolutionStepValue(KSO.THICKNESS, 0, thickness)


        # self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.THICKNESS_CONTROL_UPDATE, KSO.THICKNESS_CONTROL_CHANGE)
        # self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.THICKNESS_UPDATE, KSO.THICKNESS_CHANGE)
        return

    # --------------------------------------------------------------------------
    def __mapDesignVariableVectorToNodalVariable(self, design_variable_vector, gradient_variable):

        property_dict = dict()
        element_gradient_dict = dict()
        self.__mapDesignVariableVectorToPropertyDict(design_variable_vector, property_dict)
        self.__mapPropertyGradientToElement(property_dict, element_gradient_dict)
        self.__mapElementGradientToNode(element_gradient_dict, gradient_variable)


    # --------------------------------------------------------------------------
    def __mapDesignVariableVectorToPropertyDict(self, design_variable_vector, property_dict):

        for i in range(len(self.property_ids)):
            key = self.property_ids[i]
            property_dict[key] = design_variable_vector[i]

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

        self.__mapElementGradientToNode(element_thicknesses, KSO.THICKNESS)

# ==============================================================================
