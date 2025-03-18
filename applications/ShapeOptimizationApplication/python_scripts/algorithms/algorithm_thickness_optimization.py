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
from KratosMultiphysics.ShapeOptimizationApplication.loggers import data_logger_factory
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer

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
        if optimization_settings["design_variables"][0].Has("projection"):
            self.projection = optimization_settings["design_variables"][0]["projection"]
        else:
            self.projection = False

        if self.projection:
            self.thickness_targets = optimization_settings["design_variables"][0]["projection_settings"]["available_thicknesses"].GetVector()
            self.beta = optimization_settings["design_variables"][0]["projection_settings"]["initial_beta"].GetDouble()
            self.beta_max = optimization_settings["design_variables"][0]["projection_settings"]["max_beta"].GetDouble()
            self.q = optimization_settings["design_variables"][0]["projection_settings"]["increase_beta_factor"].GetDouble()

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.design_surface = None
        self.data_logger = None
        self.optimization_utilities = None

        self.objectives = optimization_settings["objectives"]
        self.objective_gradient_dict = {
            "gradient": 0.0,
            "projected_gradient": 0.0
        }
        self.constraints = optimization_settings["constraints"]
        self.constraint_gradients_dict = {}
        for constraint in self.constraints:
            self.constraint_gradients_dict.update({
                constraint["identifier"].GetString() : {
                    "gradient": 0.0,
                    "projected_gradient": 0.0
                }
            })
        self.max_correction_share = self.algorithm_settings["max_correction_share"].GetDouble()

        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_SEARCH_DIRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CORRECTION)

        self.number_of_design_variables = None
        self.property_ids = None

        self.thickness_dict = OrderedDict()
        self.thickness_update_dict = OrderedDict()
        self.control_thickness_dict = OrderedDict()
        self.control_update_dict = OrderedDict()

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Gradient projection algorithm only supports one objective function!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.model_part_controller.Initialize()

        self.analyzer.InitializeBeforeOptimizationLoop()

        self.design_surface = self.model_part_controller.GetDesignSurface()

        self.data_logger = data_logger_factory.CreateDataLogger(self.model_part_controller, self.communicator, self.optimization_settings)
        self.data_logger.InitializeDataLogging()

        self.optimization_utilities = KSO.OptimizationUtilities

        thickness_dict = dict()
        for condition in self.design_surface.Conditions:
            property = condition.Properties
            if property.Has(KM.THICKNESS):
                if property.Id not in thickness_dict:
                    thickness_dict[property.Id] = property.GetValue(KM.THICKNESS)

        self.number_of_design_variables = len(thickness_dict.values())

        self.thickness_dict = OrderedDict(sorted(thickness_dict.items()))
        self.control_thickness_dict = self.thickness_dict.copy()
        self.property_ids = list(self.thickness_dict.keys())

        self.__mapPropertyToNode(self.thickness_dict, KSO.THICKNESS)
        self.__mapPropertyToNode(self.control_thickness_dict, KSO.THICKNESS_CONTROL)

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
                self.__updateBeta()
                self.__determineAbsoluteChanges()

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.data_logger.FinalizeDataLogging()
        self.analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeNewThickness(self):
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)

        for property in self.optimization_model_part.Properties:
            if property.Id in self.thickness_update_dict:
                property.SetValue(KM.THICKNESS, self.thickness_dict[property.Id])

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
        objElementGradientDict = self.communicator.getStandardizedThicknessGradient(self.objectives[0]["identifier"].GetString())
        objGradientDict = dict()
        self.__mapElementToProperty(objElementGradientDict, objGradientDict)

        self.objective_gradient_dict["gradient"] = OrderedDict(sorted(objGradientDict.items()))

        # for visualization on mesh
        self.__mapPropertyToNode(self.objective_gradient_dict["gradient"], KSO.DF1DT)

        # project and damp constraint gradients
        for itr, constraint in enumerate(self.constraints):
            con_id = constraint["identifier"].GetString()
            conElementGradientDict = self.communicator.getStandardizedThicknessGradient(con_id)
            conGradientDict = dict()
            self.__mapElementToProperty(conElementGradientDict, conGradientDict)
            self.constraint_gradients_dict[con_id]["gradient"] = OrderedDict(sorted(conGradientDict.items()))

            # for visualization on mesh
            gradient_variable = KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT")
            self.__mapPropertyToNode(self.constraint_gradients_dict[con_id]["gradient"], gradient_variable)

    # --------------------------------------------------------------------------
    def __computeThicknessUpdate(self):

        ### 1. Project gradients to control thickness
        self.objective_gradient_dict["projected_gradient"] = self.__ProjectGradients(self.objective_gradient_dict["gradient"])

        for itr, constraint in enumerate(self.constraints):
            con_id = constraint["identifier"].GetString()
            self.constraint_gradients_dict[con_id]["projected_gradient"] = self.__ProjectGradients(
                self.constraint_gradients_dict[con_id]["gradient"])

        ### 2. Compute control thickness update via gradient projection
        control_thickness_update = self.__computeControlThicknessUpdate()
        for itr, key in enumerate(self.control_thickness_dict):
            self.control_update_dict[key] = control_thickness_update[itr]
            self.control_thickness_dict[key] += self.control_update_dict[key]

        ### 3. Project control thickness to obtain physical thickness
        thickness = self.__ProjectThickness(self.control_thickness_dict)
        for key in self.thickness_dict:
            self.thickness_update_dict[key] = thickness[key] - self.thickness_dict[key]
            self.thickness_dict[key] = thickness[key]

        self.__mapPropertyToNode(self.thickness_update_dict, KSO.THICKNESS_UPDATE)

    # --------------------------------------------------------------------------
    def __computeControlThicknessUpdate(self):
        """adapted from https://msulaiman.org/onewebmedia/GradProj_2.pdf"""

        prev_s = KM.Vector(self.number_of_design_variables)

        max_inner_iter = 10
        for inner_iter in range(max_inner_iter):
            KM.Logger.PrintInfo("ThicknessOpt", f"Inner Iteration: {inner_iter+1}")

            g_a, g_a_gradients = self.__getActiveConstraints()

            KM.Logger.PrintInfo("ShapeOpt", "Assemble vector of objective gradient.")
            nabla_f = KM.Vector(list(self.objective_gradient_dict["projected_gradient"].values()))

            s = KM.Vector(self.number_of_design_variables)

            if len(g_a) == 0:
                KM.Logger.PrintInfo("ShapeOpt", "No constraints active, use negative objective gradient as search direction.")
                s = nabla_f * (-1.0)

                if self.projection:
                    s = self.__ProjectSearchDirectionAndGradients(s)

                if not self.projection or (s - prev_s).norm_inf() == 0.0 or inner_iter == max_inner_iter - 1:
                    if s.norm_inf() > 0:
                        s *= self.step_size / s.norm_inf()

                    # for visualization
                    self.__mapDesignVariableVectorToNodalVariable(s, KSO.THICKNESS_SEARCH_DIRECTION)
                    self.__mapDesignVariableVectorToNodalVariable([0.0]*len(s), KSO.THICKNESS_CORRECTION)
                    self.__mapDesignVariableVectorToNodalVariable(s, KSO.THICKNESS_CONTROL_UPDATE)

                    # for visualization on mesh
                    objGradientDict = self.objective_gradient_dict["projected_gradient"]
                    self.__mapPropertyToNode(objGradientDict, KSO.DF1DT_PROJECTED)
                    for itr, constraint in enumerate(self.constraints):
                        con_id = constraint["identifier"].GetString()
                        gradient_variable = KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT_PROJECTED")
                        conGradientDict = self.constraint_gradients_dict[con_id]["projected_gradient"]
                        self.__mapPropertyToNode(conGradientDict, gradient_variable)

                    return s
                else:
                    prev_s = s
                    continue

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

            if self.projection:
                s = self.__ProjectSearchDirectionAndGradients(s)

            if not self.projection or (s - prev_s).norm_inf() == 0.0 or inner_iter == max_inner_iter - 1:

                if c.norm_inf() != 0.0:
                    if c.norm_inf() <= self.max_correction_share * self.step_size:
                        delta = self.step_size - c.norm_inf()
                        if s.norm_inf() > 0:
                            s *= delta/s.norm_inf()
                    else:
                        KM.Logger.PrintWarning("ShapeOpt", f"Correction is scaled down from {c.norm_inf()} to {self.max_correction_share * self.step_size}.")
                        c *= self.max_correction_share * self.step_size / c.norm_inf()
                        if s.norm_inf() > 0:
                            s *= (1.0 - self.max_correction_share) * self.step_size / s.norm_inf()
                else:
                    if s.norm_inf() > 0:
                        s *= self.step_size / s.norm_inf()

                # for visualization
                self.__mapDesignVariableVectorToNodalVariable(s, KSO.THICKNESS_SEARCH_DIRECTION)
                self.__mapDesignVariableVectorToNodalVariable(c, KSO.THICKNESS_CORRECTION)
                self.__mapDesignVariableVectorToNodalVariable(s+c, KSO.THICKNESS_CONTROL_UPDATE)

                # for visualization on mesh
                objGradientDict = self.objective_gradient_dict["projected_gradient"]
                self.__mapPropertyToNode(objGradientDict, KSO.DF1DT_PROJECTED)
                for itr, constraint in enumerate(self.constraints):
                    con_id = constraint["identifier"].GetString()
                    gradient_variable = KM.KratosGlobals.GetVariable(f"DC{(itr+1)}DT_PROJECTED")
                    conGradientDict = self.constraint_gradients_dict[con_id]["projected_gradient"]
                    self.__mapPropertyToNode(conGradientDict, gradient_variable)

                return s+c

            else:
                prev_s = s
                continue

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
                    KM.Vector(list(self.constraint_gradients_dict[identifier]["projected_gradient"].values())))

        return active_constraint_values, active_constraint_gradients

    # --------------------------------------------------------------------------
    def __isConstraintActive(self, constraint):
        identifier = constraint["identifier"].GetString()
        constraint_value = self.communicator.getStandardizedValue(identifier)
        if constraint["type"].GetString() == "=" or constraint_value >= 0:
            gradient = KM.Vector(list(self.constraint_gradients_dict[identifier]["projected_gradient"].values()))
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
    def __updateBeta(self):

        if self.projection:
            self.beta *= self.q
            if self.beta > self.beta_max:
                self.beta = self.beta_max

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges(self):

        for node in self.design_surface.Nodes:
            absolute_control_update = node.GetSolutionStepValue(KSO.THICKNESS_CONTROL_CHANGE)
            absolute_control_update += node.GetSolutionStepValue(KSO.THICKNESS_CONTROL_UPDATE)
            node.SetSolutionStepValue(KSO.THICKNESS_CONTROL_CHANGE, 0, absolute_control_update)

            control_thickness = node.GetSolutionStepValue(KSO.THICKNESS_CONTROL)
            control_thickness += node.GetSolutionStepValue(KSO.THICKNESS_CONTROL_UPDATE)
            node.SetSolutionStepValue(KSO.THICKNESS_CONTROL, 0, control_thickness)

            absolute_update = node.GetSolutionStepValue(KSO.THICKNESS_CHANGE)
            absolute_update += node.GetSolutionStepValue(KSO.THICKNESS_UPDATE)
            node.SetSolutionStepValue(KSO.THICKNESS_CHANGE, 0, absolute_update)

            thickness = node.GetSolutionStepValue(KSO.THICKNESS)
            thickness += node.GetSolutionStepValue(KSO.THICKNESS_UPDATE)
            node.SetSolutionStepValue(KSO.THICKNESS, 0, thickness)

        return

    # --------------------------------------------------------------------------
    def __mapDesignVariableVectorToNodalVariable(self, design_variable_vector, gradient_variable):

        property_dict = dict()
        element_gradient_dict = dict()
        self.__mapDesignVariableVectorToPropertyDict(design_variable_vector, property_dict)
        self.__mapPropertyToElement(property_dict, element_gradient_dict)
        self.__mapElementDataToNode(element_gradient_dict, gradient_variable)

    # --------------------------------------------------------------------------
    def __mapDesignVariableVectorToPropertyDict(self, design_variable_vector, property_dict):

        for i in range(len(self.property_ids)):
            key = self.property_ids[i]
            property_dict[key] = design_variable_vector[i]

    # --------------------------------------------------------------------------
    def __mapPropertyToNode(self, property_dict, nodal_variable):
        element_dict = dict()
        self.__mapPropertyToElement(property_dict, element_dict)
        self.__mapElementDataToNode(element_dict, nodal_variable)

    # --------------------------------------------------------------------------
    def __mapElementToProperty(self, element_dict, property_dict):

        for condition in self.design_surface.Conditions:
            if condition.Properties.Id not in property_dict:
                property_dict[condition.Properties.Id] = 0.0
            property_dict[condition.Properties.Id] += element_dict[condition.Id]

    # --------------------------------------------------------------------------
    def __mapPropertyToElement(self, property_dict, element_dict):

        for condition in self.design_surface.Conditions:
            element_dict[condition.Id] = property_dict[condition.Properties.Id]

    # --------------------------------------------------------------------------
    def __mapElementDataToNode(self, element_data_dict, nodal_variable):

        # reset variables
        for node in self.design_surface.Nodes:
            node.SetSolutionStepValue(nodal_variable, 0, 0.0)

        for condition in self.design_surface.Conditions:
            condition.SetValue(nodal_variable, element_data_dict[condition.Id])

        total_node_areas = dict()
        for condition in self.design_surface.Conditions:
            element_data = condition.GetValue(nodal_variable)
            for node in condition.GetNodes():
                if node.Id in total_node_areas:
                    total_node_areas[node.Id] += condition.GetGeometry().Area()
                else:
                    total_node_areas[node.Id] = condition.GetGeometry().Area()
                node_data = node.GetSolutionStepValue(nodal_variable)
                node_data += element_data * condition.GetGeometry().Area()
                node.SetSolutionStepValue(nodal_variable, 0, node_data)

        for node in self.design_surface.Nodes:
            total_node_area = total_node_areas[node.Id]
            node_data = node.GetSolutionStepValue(nodal_variable)
            node_data /= total_node_area
            node.SetSolutionStepValue(nodal_variable, 0, node_data)

    # --------------------------------------------------------------------------
    def __ProjectGradients(self, gradient_dict):

        if not self.projection:
            return gradient_dict

        projected_gradient_dict = OrderedDict()

        for id, t in self.control_thickness_dict.items():
            t_m = self.___GetInterval(t)

            w = (t - t_m[0]) / (t_m[1] - t_m[0])
            tanh_numerator = 2 * np.tanh(self.beta * 0.5)
            dw_dt_p = 1 / (t_m[1] - t_m[0])
            dt_p_dw = ((t_m[1] - t_m[0]) / tanh_numerator) \
                * (1 - np.tanh(self.beta * (w - 0.5))) * self.beta
            projected_gradient_dict[id] = gradient_dict[id] * dw_dt_p * dt_p_dw

        return projected_gradient_dict

    def __ProjectSearchDirectionAndGradients(self, search_direction):

        if not self.projection:
            return search_direction

        for itr, key in enumerate(self.control_thickness_dict.keys()):
            t = self.control_thickness_dict[key]

            if t <= self.thickness_targets[0] and search_direction[itr] < 0.0:
                search_direction[itr] = 0.0
                self.objective_gradient_dict["projected_gradient"][key] = 0.0
                for constraint_gradient in self.constraint_gradients_dict.values():
                    constraint_gradient["projected_gradient"][key] = 0.0

            elif t >= self.thickness_targets[len(self.thickness_targets)-1] and search_direction[itr] > 0.0:
                search_direction[itr] = 0.0
                self.objective_gradient_dict["projected_gradient"][key] = 0.0
                for constraint_gradient in self.constraint_gradients_dict.values():
                    constraint_gradient["projected_gradient"][key] = 0.0

        return search_direction

    def __ProjectThickness(self, control_thickness_dict):

        if not self.projection:
            return control_thickness_dict

        thickness_dict = OrderedDict()

        for id, t in control_thickness_dict.items():

            if t <= self.thickness_targets[0]:
                thickness_dict[id] = self.thickness_targets[0]

            elif t >= self.thickness_targets[len(self.thickness_targets)-1]:
                thickness_dict[id] = self.thickness_targets[len(self.thickness_targets)-1]

            else:
                t_m = self.___GetInterval(t)
                w = (t - t_m[0]) / (t_m[1] - t_m[0])
                tanh_denominator = np.tanh(self.beta * 0.5) + np.tanh(self.beta * (w - 0.5))
                tanh_numerator = 2 * np.tanh(self.beta * 0.5)
                tanh_term = tanh_denominator / tanh_numerator
                thickness_dict[id] = t_m[0] + (t_m[1] - t_m[0]) * tanh_term

        return thickness_dict

    def ___GetInterval(self, t):

        if t <= self.thickness_targets[0]:
            return np.array([self.thickness_targets[0], self.thickness_targets[1]])
        if t >= self.thickness_targets[len(self.thickness_targets)-1]:
            return np.array([self.thickness_targets[len(self.thickness_targets)-2], self.thickness_targets[len(self.thickness_targets)-1]])

        for i in range(len(self.thickness_targets)-1):
            if t >= self.thickness_targets[i] and t <= self.thickness_targets[i+1]:
                return np.array([self.thickness_targets[i], self.thickness_targets[i+1]])
# ==============================================================================
