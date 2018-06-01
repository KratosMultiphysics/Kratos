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
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from algorithm_base import OptimizationAlgorithm
from custom_math import NormInf3D, DotProduct, ScalarVectorProduct
from custom_variable_utilities import WriteDictionaryDataOnNodalVariable, ReadNodalVariableToList, WriteNodeCoordinatesToList, WriteListToNodalVariable
from custom_timer import Timer
import mapper_factory

# ==============================================================================
class AlgorithmTrustRegion(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = Parameters("""
        {
            "name"                        : "trust_region",
            "max_step_length"             : 1.0,
            "step_length_tolerance"       : 1e-3,
            "step_length_decrease_factor" : 2.0,
            "step_length_increase_factor" : 1.2,
            "feasibility_frequency"       : 5,
            "max_iterations"              : 10
        }""")
        self.algorithm_settings =  optimization_settings["optimization_algorithm"]
        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.objectives, self.equality_constraints, self.inequality_constraints = communicator.GetInfoAboutResponses()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surface = model_part_controller.GetDesignSurface()

        self.mapper = mapper_factory.CreateMapper(self.design_surface, optimization_settings["design_variables"]["filter"])

        self.GeometryUtilities = GeometryUtilities(self.design_surface)

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if len(self.objectives) > 1:
            raise RuntimeError("Trust-region algorithm only supports one objective function!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.max_step_length = self.algorithm_settings["max_step_length"].GetDouble()
        self.step_length_tolerance = self.algorithm_settings["step_length_tolerance"].GetDouble()
        self.step_length_decrease_factor = self.algorithm_settings["step_length_decrease_factor"].GetDouble()
        self.step_length_increase_factor = self.algorithm_settings["step_length_increase_factor"].GetDouble()
        self.feasibility_frequency = self.algorithm_settings["feasibility_frequency"].GetDouble()
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt()

        self.algo_is_in_termination_phase = False
        self.algo_is_in_feasibility_phase = False
        self.next_iteration_to_enforce_feasibility = 9999
        self.last_itr_with_step_length_update = 0

        self.number_of_design_variables = 3*self.design_surface.NumberOfNodes()
        self.x_init = WriteNodeCoordinatesToList(self.design_surface)

        self.value_history = {}
        self.value_history["val_obj"] = []
        self.value_history["val_eqs"] = []
        self.value_history["val_ineqs"] = []
        self.value_history["len_eqs"] = []
        self.value_history["len_ineqs"] = []
        self.value_history["step_length"] = []

        self.model_part_controller.ImportOptimizationModelPart()
        self.model_part_controller.InitializeMeshController()
        self.mapper.InitializeMapping()
        self.analyzer.InitializeBeforeOptimizationLoop()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        for self.opt_iteration in range(1,self.algorithm_settings["max_iterations"].GetInt()):
            print("\n>===================================================================")
            print("> ",timer.GetTimeStamp(),": Starting optimization iteration ",self.opt_iteration)
            print(">===================================================================\n")

            timer.StartNewLap()

            self.__InitializeNewShape()

            self.__AnalyzeShape()

            val_obj, val_eqs, val_ineqs = self.__ProcessAnalysisResults()

            if self.__CheckConvergence(val_eqs, val_ineqs):
                break

            dir_objs, len_eqs, dir_eqs, len_ineqs, dir_ineqs =  self.__ConvertAnalysisResultsToLengthDirectionFormat()

            step_length = self.__ApplyStepLengthRule(val_obj, len_eqs, len_ineqs)

            len_bar_eqs = ScalarVectorProduct(1/step_length, len_eqs)
            len_bar_ineqs = ScalarVectorProduct(1/step_length, len_ineqs)

            # # convert (stepLength,lInequality,lEquality) into delta_x
            # deltaXBar,lambdaBarObjective,lambdaBarInequality,lambdaBarEquality,sInit,sDx = stepDirectionRule(eObjective,lBarInequality,eInequality,lBarEquality,eEquality,self)
            # deltaX = [stepLength*deltaXBar[i] for i in range(self.n)]
            dx = [0.0 for i in range(self.number_of_design_variables)]

            self.__UpdateMesh(dx)

            self.value_history["val_obj"].append(val_obj)
            self.value_history["val_eqs"].append(val_ineqs)
            self.value_history["val_ineqs"].append(val_ineqs)
            self.value_history["len_eqs"].append(len_eqs)
            self.value_history["len_ineqs"].append(len_ineqs)
            self.value_history["step_length"].append(step_length)

            print("--------------------------------------------")
            print("val_obj = ", val_obj)
            print("val_eqs = ", val_eqs)
            print("val_ineqs = ", val_ineqs)
            print("len_eqs = ", len_eqs)
            print("len_ineqs = ", len_ineqs)
            print("step_length = ", step_length)
            print("--------------------------------------------")

            print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
            print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __InitializeNewShape(self):
            self.model_part_controller.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
            self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __AnalyzeShape(self):
            self.communicator.initializeCommunication()

            for obj in self.objectives:
                obj_id = obj["settings"]["identifier"].GetString()
                self.communicator.requestValueOf(obj_id)
                self.communicator.requestGradientOf(obj_id)

            for eq in self.equality_constraints:
                eq_id = eq["settings"]["identifier"].GetString()
                self.communicator.requestValueOf(eq_id)
                self.communicator.requestGradientOf(eq_id)

            for ineq in self.inequality_constraints:
                ineq_id = ineq["settings"]["identifier"].GetString()
                self.communicator.requestValueOf(ineq_id)
                self.communicator.requestGradientOf(ineq_id)

            self.analyzer.AnalyzeDesignAndReportToCommunicator(self.design_surface, self.opt_iteration, self.communicator)

            # Reset possible shape modifications during analysis
            self.model_part_controller.SetMeshToReferenceMesh()
            self.model_part_controller.SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    def __ProcessAnalysisResults(self):
        # Process objective results
        obj_value = None

        obj = self.objectives[0]
        obj_number = obj["number"]
        obj_id = obj["settings"]["identifier"].GetString()

        obj_value = self.communicator.getStandardizedValue(obj_id)
        obj_gradients_dict = self.communicator.getStandardizedGradient(obj_id)

        nodal_variable = KratosGlobals.GetVariable("DF"+str(obj_number+1)+"DX")
        WriteDictionaryDataOnNodalVariable(obj_gradients_dict, self.optimization_model_part, nodal_variable)

        if obj["settings"]["project_gradient_on_surface_normals"].GetBool():
            self.__PerformProjectionOnNormals(nodal_variable)

        nodal_variable_mapped = KratosGlobals.GetVariable("DF"+str(obj_number+1)+"DX_MAPPED")
        self.__PerformMapping(nodal_variable, nodal_variable_mapped)

        # Process equality constraint results
        eq_values = []
        for eq in self.equality_constraints:
            eq_number = eq["number"]
            eq_id = eq["settings"]["identifier"].GetString()

            eq_values.append(self.communicator.getStandardizedValue(eq_id))
            eq_gradients_dict = self.communicator.getStandardizedGradient(eq_id)

            nodal_variable = KratosGlobals.GetVariable("DC"+str(eq_number+1)+"DX")
            WriteDictionaryDataOnNodalVariable(eq_gradients_dict, self.optimization_model_part, nodal_variable)

            if eq["settings"]["project_gradient_on_surface_normals"].GetBool():
                self.__PerformProjectionOnNormals(nodal_variable)

            nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(eq_number+1)+"DX_MAPPED")
            self.__PerformMapping(nodal_variable, nodal_variable_mapped)

        # Process inequality constraint results
        ineq_values = []
        for ineq in self.inequality_constraints:
            ineq_number = ineq["number"]
            ineq_id = ineq["settings"]["identifier"].GetString()

            ineq_values.append(self.communicator.getStandardizedValue(ineq_id))
            ineq_gradients_dict = self.communicator.getStandardizedGradient(ineq_id)

            nodal_variable = KratosGlobals.GetVariable("DC"+str(ineq_number+1)+"DX")
            WriteDictionaryDataOnNodalVariable(ineq_gradients_dict, self.optimization_model_part, nodal_variable)

            if self.ineq["settings"]["project_gradient_on_surface_normals"].GetBool():
                self.__PerformProjectionOnNormals(nodal_variable)

            nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(ineq_number+1)+"DX_MAPPED")
            self.__PerformMapping(nodal_variable, nodal_variable_mapped)

        return obj_value, eq_values, ineq_values

    # --------------------------------------------------------------------------
    def __PerformMapping(self, nodal_variable, nodal_variable_mapped):
        self.mapper.MapToDesignSpace(nodal_variable, nodal_variable_mapped)
        self.mapper.MapToGeometrySpace(nodal_variable_mapped, nodal_variable_mapped)

    # --------------------------------------------------------------------------
    def __PerformProjectionOnNormals(self, nodal_variable):
        self.GeometryUtilities.ComputeUnitSurfaceNormals()
        self.GeometryUtilities.ProjectNodalVariableOnUnitSurfaceNormals(nodal_variable)

    # --------------------------------------------------------------------------
    def __CheckConvergence(self, val_eqs, val_ineqs):

        # Check termination conditions
        if self.opt_iteration == self.max_iterations:
            self.algo_is_in_termination_phase = True
            print("> Algorithm: Begin enforcing feasibility for termination (max number of optimization iterations reached)!")

        if self.opt_iteration > 1 and self.value_history["step_length"][-1] < self.step_length_tolerance:
            self.algo_is_in_termination_phase = True
            print("> Algorithm: Begin enforcing feasibility for termination (step length < tolerance)!")

        # Check feasibility
        if self.algo_is_in_termination_phase:
            self.algo_is_in_feasibility_phase = True

        if self.opt_iteration == self.next_iteration_to_enforce_feasibility:
            self.algo_is_in_feasibility_phase = True
            self.reduceOscillationFactor = 1
            print("> Algorithm: Begin enforcing feasibility according specified frequency!")

        feasibility_is_enforced = self.__IsFeasibilityEnforced(val_eqs, val_ineqs)

        if feasibility_is_enforced and self.algo_is_in_feasibility_phase:
            print("> Algorithm: Feasibility was enforced")
            self.algo_is_in_feasibility_phase = False
            self.next_iteration_to_enforce_feasibility = iterCounter + self.feasibility_frequency

        if feasibility_is_enforced and not self.algo_is_in_feasibility_phase:
            print("> Algorithm: Feasibility is spontaneously enforced.")

        if self.algo_is_in_termination_phase and feasibility_is_enforced:
            print("> Algorithm: Converged!")
            return True
        else:
            False

    # --------------------------------------------------------------------------
    def __IsFeasibilityEnforced(self, val_eqs, val_ineqs):
        # Return false if any equality constraint is violated
        for itr, eq in enumerate(self.equality_constraints):
            eq_value = val_eqs[itr]
            eq_id = eq["settings"]["identifier"].GetString()
            eq_tolerance = eq["settings"]["tolerance"].GetDouble()
            eq_reference = self.communicator.getReferenceValue(eq_id)
            if abs(eq_value) > eq_tolerance*eq_reference:
                return False

        # Return false if any inequality constraint is violated
        for val_ineq in val_ineqs:
            if val_ineq > 0:
                return False

        return True

    # --------------------------------------------------------------------------
    def __ConvertAnalysisResultsToLengthDirectionFormat(self):
        # Convert objective results
        dir_objs = []
        for obj in self.objectives:
            obj_number = obj["number"]
            obj_id = obj["settings"]["identifier"].GetString()

            nodal_variable = KratosGlobals.GetVariable("DF"+str(obj_number+1)+"DX")
            nodal_variable_mapped = KratosGlobals.GetVariable("DF"+str(obj_number+1)+"DX_MAPPED")

            obj_gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
            obj_gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

            dir_obj = self.__ConvertToLengthDirectionFormat(obj_gradient, obj_gradient_mapped)
            dir_objs.append(dir_obj)

        # Convert equality constraint results
        len_eqs = []
        dir_eqs = []
        for eq in self.equality_constraints:
            eq_number = eq["number"]
            eq_id = eq["settings"]["identifier"].GetString()

            nodal_variable = KratosGlobals.GetVariable("DC"+str(eq_number+1)+"DX")
            nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(eq_number+1)+"DX_MAPPED")

            eq_value = self.communicator.getStandardizedValue(eq_id)
            eq_gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
            eq_gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

            dir_eq, len_eq = self.__ConvertToLengthDirectionFormat(eq_gradient, eq_gradient_mapped, eq_value)

            dir_eqs.append(dir_eq)
            len_eqs.append(len_eq)

        # Convert inequality constraint results
        len_ineqs = []
        dir_ineqs = []
        for ineq in self.inequality_constraints:
            ineq_number = ineq["number"]
            ineq_id = ineq["settings"]["identifier"].GetString()

            nodal_variable = KratosGlobals.GetVariable("DC"+str(ineq_number+1)+"DX")
            nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(ineq_number+1)+"DX_MAPPED")

            ineq_value = self.communicator.getStandardizedValue(id)
            ineq_gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
            ineq_gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

            dir_ineq, len_ineq = self.__ConvertToLengthDirectionFormat(ineq_gradient, ineq_gradient_mapped, ineq_value)

            dir_ineqs.append(dir_ineq)
            len_ineqs.append(len_ineq)

        return dir_objs, len_eqs, dir_eqs, len_ineqs, dir_ineqs

    # --------------------------------------------------------------------------
    def __ConvertToLengthDirectionFormat(self, gradient, modified_gradient, value=None):
        norm_inf = NormInf3D(modified_gradient)
        direction = [-modified_gradient[itr]/norm_inf for itr in range(len(modified_gradient))]

        if value == None:
            return direction
        else:
            grad_dot_dir = DotProduct(gradient, direction)
            length = -value/grad_dot_dir
            return direction, length

    # --------------------------------------------------------------------------
    def __ApplyStepLengthRule(self, val_obj, len_eqs, len_ineqs):

        step_length = self.max_step_length

        # dont change step lenght with final feasibility is enforced
        if self.algo_is_in_feasibility_phase:
            self.last_itr_with_step_length_update = self.opt_iteration
            return step_length

        # Wait with an updat of the step length for min three iterations
        if self.opt_iteration<3 or self.opt_iteration<self.last_itr_with_step_length_update+2:
            return step_length

        # Decrease step length if necessary conditions are met
        if self.__IsStepLengthToBeDecreased(val_obj, len_eqs, len_ineqs):
            self.last_itr_with_step_length_update = self.opt_iteration
            print("> Algorithm: Decreasing step length!")
            return step_length / self.step_length_decrease_factor

        # Increase slowly step length if necessary conditions are met
        if self.__IsStepLengthToBeIncreased(val_obj, len_eqs, len_ineqs):
            self.last_itr_with_step_length_update = self.opt_iteration
            print("> Algorithm: Increasing step length")
            return min(self.max_step_length, step_length * self.step_length_increase_factor)

        #do nothing
        return step_length

    # --------------------------------------------------------------------------
    def __IsStepLengthToBeDecreased(self, val_obj, len_eqs, len_ineqs):

        # Check if enough values are already computed
        if self.opt_iteration < 3:
            return False

        v2, v1, v0 = self.value_history["val_obj"][-2:] + [val_obj]

        # Check if objective oscillates
        if (v1 < v2 and v1 < v0) and v0 > (v2-0.1*(v2-v1)):
            objective_oscillates = True
        else:
            objective_oscillates = False

        # Objective does not improve
        if v0>v2:
            no_progress_in_obj = True
        else:
            no_progress_in_obj = False

        # No progress in the distance to the feasible domain
        s_old = self.value_history["step_length"][-2]

        old_len_eqs = self.value_history["len_eqs"][-2]
        old_len_ineqs = self.value_history["len_ineqs"][-2]

        dist_feasible_domoain = self.__SimpleDistanceToFeasibleDomain(len_eqs,len_ineqs)
        dist_feasible_domoain_old = self.__SimpleDistanceToFeasibleDomain(old_len_eqs , old_len_ineqs)

        if dist_feasible_domoain > (dist_feasible_domoain_old - 0.2*s_old):
            no_progress_towards_feasible_domain = True
        else:
            no_progress_towards_feasible_domain = False

        # Evaluate individual conditions
        if (no_progress_in_obj or objective_oscillates) and no_progress_towards_feasible_domain:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __SimpleDistanceToFeasibleDomain(self, len_equality, len_inequality):
        return max( [max(l,0) for l in len_inequality] + [abs(l) for l in len_equality] + [0] ) # add [0] for case with no constraint

    # --------------------------------------------------------------------------
    def __IsStepLengthToBeIncreased(self, val_obj, len_eqs, len_ineqs):
        # Progress of objective function is monotonic over a few iterations
        number_of_smooth_steps = 5
        if self.opt_iteration < number_of_smooth_steps:
            return False
        else:
            val_hist = self.value_history["val_obj"][-number_of_smooth_steps+1:] + [val_obj]

            # Check if values are changing monotonically
            is_oscillating = []
            for i in range(number_of_smooth_steps-2):
                is_oscillating.append( (val_hist[i+1]-val_hist[i])*(val_hist[i+2]-val_hist[i+1])<0 )

            if sum(is_oscillating)==0: # no oscillation at all
                return True
            else:
                return False

    # --------------------------------------------------------------------------
    def __UpdateMesh(self, dx):
        WriteListToNodalVariable(dx, self.design_surface, SHAPE_UPDATE)
        self.model_part_controller.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
        self.model_part_controller.SetReferenceMeshToMesh()

        x = WriteNodeCoordinatesToList(self.design_surface)
        dxAbsolute = [x[i]+dx[i] - self.x_init[i] for i in range(self.number_of_design_variables)]
        WriteListToNodalVariable(dxAbsolute, self.design_surface, SHAPE_CHANGE)

# ==============================================================================
