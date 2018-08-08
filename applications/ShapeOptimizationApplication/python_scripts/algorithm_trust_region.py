# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from algorithm_base import OptimizationAlgorithm
from custom_math import NormInf3D, Dot, ScalarVectorProduct, Norm2, RowSize, HorzCat, Minus, Plus, Trans, Prod, ZeroVector, ZeroMatrix, IsEmpty, SafeConvertVectorToMatrix
from custom_math import TranslateToNewBasis, TranslateToOriginalBasis, QuadProg, PerformBisectioning, SolveLinearSystem
from custom_variable_utilities import WriteDictionaryDataOnNodalVariable, ReadNodalVariableToList, WriteNodeCoordinatesToList, WriteListToNodalVariable
from custom_timer import Timer
import mapper_factory
import data_logger_factory
import copy

# ==============================================================================
class AlgorithmTrustRegion(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = Parameters("""
        {
            "name"                          : "trust_region",
            "max_step_length"               : 1.0,
            "step_length_tolerance"         : 1e-3,
            "step_length_reduction_factor"  : 0.5,
            "min_share_objective"           : 0.1,
            "max_iterations"                : 10,
            "far_away_length"               : 2.0,
            "subopt_max_itr"                : 50,
            "subopt_tolerance"              : 1e-10,
            "bisectioning_max_itr"          : 30,
            "bisectioning_tolerance"        : 1e-2,
            "obj_share_during_correction"   : 1
        }""")
        self.algorithm_settings =  optimization_settings["optimization_algorithm"]
        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

        self.specified_objectives = optimization_settings["objectives"]
        self.specified_constraints = optimization_settings["constraints"]

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surface = model_part_controller.GetDesignSurface()

        self.mapper = mapper_factory.CreateMapper(self.design_surface, optimization_settings["design_variables"]["filter"])
        self.data_logger = data_logger_factory.CreateDataLogger(model_part_controller, communicator, optimization_settings)

        self.geometry_utilities = GeometryUtilities(self.design_surface)
        self.optimization_utilities = OptimizationUtilities(self.design_surface, optimization_settings)

        self.is_damping_specified = optimization_settings["design_variables"]["damping"]["perform_damping"].GetBool()
        if self.is_damping_specified:
            damping_regions = self.model_part_controller.GetDampingRegions()
            self.damping_utilities = DampingUtilities(self.design_surface, damping_regions, optimization_settings)

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.specified_objectives.size() > 1:
            raise RuntimeError("Trust-region algorithm only supports one objective function!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.number_of_design_variables = 3*self.design_surface.NumberOfNodes()

        self.x_init = WriteNodeCoordinatesToList(self.design_surface)

        self.model_part_controller.ImportOptimizationModelPart()
        self.model_part_controller.InitializeMeshController()
        self.mapper.InitializeMapping()
        self.analyzer.InitializeBeforeOptimizationLoop()
        self.data_logger.InitializeDataLogging()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        for self.opt_iteration in range(1,self.algorithm_settings["max_iterations"].GetInt()+1):
            print("\n>===================================================================")
            print("> ",timer.GetTimeStamp(),": Starting optimization iteration ",self.opt_iteration)
            print(">===================================================================\n")

            timer.StartNewLap()

            self.__InitializeNewShape()

            self.__AnalyzeShape()

            self.__ResetPossibleShapeModificationsDuringAnalysis()

            self.__PostProcessGradientsObtainedFromAnalysis()

            len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs =  self.__ConvertAnalysisResultsToLengthDirectionFormat()

            step_length = self.__DetermineMaxStepLength()

            len_bar_obj, len_bar_eqs, len_bar_ineqs = self.__ExpressInStepLengthUnit(len_obj, len_eqs, len_ineqs, step_length)

            dx_bar, test_norm_dx_bar, bi_itrs, bi_err, adj_len_bar_obj, adj_len_bar_eqs, adj_len_bar_ineqs = self.__DetermineStep(len_bar_obj, dir_obj, len_bar_eqs, dir_eqs, len_bar_ineqs, dir_ineqs)

            norm_dx = self.__ComputeShapeUpdate(dx_bar, step_length)

            self.__LogCurrentOptimizationStep(step_length, len_bar_obj, len_bar_eqs, len_bar_ineqs, test_norm_dx_bar, bi_itrs, bi_err, adj_len_bar_obj, adj_len_bar_eqs, adj_len_bar_ineqs, norm_dx)

            print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
            print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.analyzer.FinalizeAfterOptimizationLoop()
        self.data_logger.FinalizeDataLogging()

    # --------------------------------------------------------------------------
    def __InitializeNewShape(self):
            self.model_part_controller.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
            self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __AnalyzeShape(self):
            self.communicator.initializeCommunication()

            obj_id = self.specified_objectives[0]["identifier"].GetString()
            self.communicator.requestValueOf(obj_id)
            self.communicator.requestGradientOf(obj_id)

            for itr in range(self.specified_constraints.size()):
                con = self.specified_constraints[itr]

                if con["type"].GetString()=="=":
                    eq_id = con["identifier"].GetString()
                    self.communicator.requestValueOf(eq_id)
                    self.communicator.requestGradientOf(eq_id)
                else:
                    ineq_id = con["identifier"].GetString()
                    self.communicator.requestValueOf(ineq_id)
                    self.communicator.requestGradientOf(ineq_id)

            self.analyzer.AnalyzeDesignAndReportToCommunicator(self.design_surface, self.opt_iteration, self.communicator)

    # --------------------------------------------------------------------------
    def __ResetPossibleShapeModificationsDuringAnalysis(self):
        self.model_part_controller.SetMeshToReferenceMesh()
        self.model_part_controller.SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    def __PostProcessGradientsObtainedFromAnalysis(self):
        self.geometry_utilities.ComputeUnitSurfaceNormals()

        # Process objective gradients
        obj = self.specified_objectives[0]
        obj_id = obj["identifier"].GetString()

        obj_gradients_dict = self.communicator.getStandardizedGradient(obj_id)

        nodal_variable = KratosGlobals.GetVariable("DF1DX")
        WriteDictionaryDataOnNodalVariable(obj_gradients_dict, self.optimization_model_part, nodal_variable)

        # Projection on surface normals
        if obj["project_gradient_on_surface_normals"].GetBool():
            self.geometry_utilities.ProjectNodalVariableOnUnitSurfaceNormals(nodal_variable)

        # Damping
        if self.is_damping_specified:
            self.damping_utilities.DampNodalVariable(nodal_variable)

        # Mapping
        nodal_variable_mapped = KratosGlobals.GetVariable("DF1DX_MAPPED")
        self.mapper.MapToDesignSpace(nodal_variable, nodal_variable_mapped)
        self.mapper.MapToGeometrySpace(nodal_variable_mapped, nodal_variable_mapped)

        # Damping
        if self.is_damping_specified:
            self.damping_utilities.DampNodalVariable(nodal_variable_mapped)

        # Process constraint gradients
        for itr in range(self.specified_constraints.size()):
            con = self.specified_constraints[itr]

            # Process equality constraints
            if con["type"].GetString()=="=":
                eq_id = con["identifier"].GetString()

                eq_gradients_dict = self.communicator.getStandardizedGradient(eq_id)

                nodal_variable = KratosGlobals.GetVariable("DC"+str(itr+1)+"DX")
                WriteDictionaryDataOnNodalVariable(eq_gradients_dict, self.optimization_model_part, nodal_variable)

                # Projection on surface normals
                if con["project_gradient_on_surface_normals"].GetBool():
                    self.geometry_utilities.ProjectNodalVariableOnUnitSurfaceNormals(nodal_variable)

                # Damping
                if self.is_damping_specified:
                    self.damping_utilities.DampNodalVariable(nodal_variable)

                # Mapping
                nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")
                self.mapper.MapToDesignSpace(nodal_variable, nodal_variable_mapped)
                self.mapper.MapToGeometrySpace(nodal_variable_mapped, nodal_variable_mapped)

                # Damping
                if self.is_damping_specified:
                    self.damping_utilities.DampNodalVariable(nodal_variable_mapped)

            # Process inequality constraints
            else:
                ineq_id = con["identifier"].GetString()

                ineq_gradients_dict = self.communicator.getStandardizedGradient(ineq_id)

                nodal_variable = KratosGlobals.GetVariable("DC"+str(itr+1)+"DX")
                WriteDictionaryDataOnNodalVariable(ineq_gradients_dict, self.optimization_model_part, nodal_variable)

                # Projection on surface normals
                if con["project_gradient_on_surface_normals"].GetBool():
                    self.geometry_utilities.ProjectNodalVariableOnUnitSurfaceNormals(nodal_variable)

                # Damping
                if self.is_damping_specified:
                    self.damping_utilities.DampNodalVariable(nodal_variable)

                # Mapping
                nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")
                self.mapper.MapToDesignSpace(nodal_variable, nodal_variable_mapped)
                self.mapper.MapToGeometrySpace(nodal_variable_mapped, nodal_variable_mapped)

                # Damping
                if self.is_damping_specified:
                    self.damping_utilities.DampNodalVariable(nodal_variable_mapped)

    # --------------------------------------------------------------------------
    def __ConvertAnalysisResultsToLengthDirectionFormat(self):
        # Convert objective results
        obj = self.specified_objectives[0]
        obj_id = obj["identifier"].GetString()

        nodal_variable = KratosGlobals.GetVariable("DF1DX")
        nodal_variable_mapped = KratosGlobals.GetVariable("DF1DX_MAPPED")

        obj_value = self.communicator.getStandardizedValue(obj_id)
        obj_gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
        obj_gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

        dir_obj, len_obj = self.__ConvertToLengthDirectionFormat(obj_value, obj_gradient, obj_gradient_mapped)
        dir_obj = dir_obj

        # Convert constraints
        len_eqs = []
        dir_eqs = []
        len_ineqs = []
        dir_ineqs = []

        for itr in range(self.specified_constraints.size()):
            con = self.specified_constraints[itr]

            # Convert equality constraints
            if con["type"].GetString()=="=":
                eq_id = con["identifier"].GetString()

                nodal_variable = KratosGlobals.GetVariable("DC"+str(itr+1)+"DX")
                nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")

                eq_value = self.communicator.getStandardizedValue(eq_id)
                eq_gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
                eq_gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

                dir_eq, len_eq = self.__ConvertToLengthDirectionFormat(eq_value, eq_gradient, eq_gradient_mapped)

                dir_eqs.append(dir_eq)
                len_eqs.append(len_eq)

            # Convert inequality constraints
            else:
                ineq_id = con["identifier"].GetString()

                nodal_variable = KratosGlobals.GetVariable("DC"+str(itr+1)+"DX")
                nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")

                ineq_value = self.communicator.getStandardizedValue(ineq_id)
                ineq_gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
                ineq_gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

                dir_ineq, len_ineq = self.__ConvertToLengthDirectionFormat(ineq_value, ineq_gradient, ineq_gradient_mapped)

                dir_ineqs.append(dir_ineq)
                len_ineqs.append(len_ineq)

        return len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs

    # --------------------------------------------------------------------------
    def __ConvertToLengthDirectionFormat(self, value, gradient, modified_gradient):
        norm_inf = NormInf3D(modified_gradient)
        if norm_inf > 1e-12:
            direction = ScalarVectorProduct(-1/norm_inf,modified_gradient)
            length = -value/Dot(gradient, direction)
        else:
            print("\nWarning! Vanishing norm-infinity for gradient detected!")
            direction = modified_gradient
            length = 0.0

        return direction, length

    # --------------------------------------------------------------------------
    def __DetermineMaxStepLength(self):
        if self.opt_iteration < 4:
            return self.algorithm_settings["max_step_length"].GetDouble()
        else:
            obj_id = self.specified_objectives[0]["identifier"].GetString()
            current_obj_val = self.communicator.getStandardizedValue(obj_id)
            obj_history = self.data_logger.GetValueHistory()[obj_id]
            step_history = self.data_logger.GetValueHistory()["step_length"]

            objective_is_oscillating = False
            is_decrease_1 = (current_obj_val - obj_history[self.opt_iteration-1])< 0
            is_decrease_2 = (obj_history[self.opt_iteration-1] - obj_history[self.opt_iteration-2])<0
            is_decrease_3 = (current_obj_val - obj_history[self.opt_iteration-3])< 0
            if (is_decrease_1 and is_decrease_2== False and is_decrease_3) or (is_decrease_1== False and is_decrease_2 and is_decrease_3==False):
                objective_is_oscillating = True

            if objective_is_oscillating:
                return step_history[self.opt_iteration-1]*self.algorithm_settings["step_length_reduction_factor"].GetDouble()
            else:
                return step_history[self.opt_iteration-1]

    # --------------------------------------------------------------------------
    def __ExpressInStepLengthUnit(self, len_obj, len_eqs, len_ineqs, step_length):
        len_bar_obj = 1/step_length * len_obj
        len_bar_eqs = ScalarVectorProduct(1/step_length, len_eqs)
        len_bar_ineqs = ScalarVectorProduct(1/step_length, len_ineqs)
        return len_bar_obj, len_bar_eqs, len_bar_ineqs

    # --------------------------------------------------------------------------
    def __DetermineStep(self, len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs):
        print("\n> Starting determination of step...")

        timer = Timer()
        timer.StartTimer()

        # Create projector object wich can do the projection in the orthogonalized subspace
        projector = Projector(len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs, self.algorithm_settings)

        # 1. Test projection if there is room for objective improvement
        # I.e., the actual step length to become feasible for an inactive threshold is smaller than 1 and hence a part of the step can be dedicated to objective improvement
        nargout = 2
        len_obj_test = 0.01
        inactive_threshold = 100
        test_norm_dX, is_projection_sucessfull = projector.RunProjection(len_obj_test, inactive_threshold, nargout)

        print("> Time needed for one projection step = ", timer.GetTotalTime(), "s")

        # 2. Determine step following two different modes depending on the previos found step length to the feasible domain
        if is_projection_sucessfull:
            if test_norm_dX < 1: # Minimizing mode
                print ("\n> Computing projection case 1...")

                func = lambda len_obj: projector.RunProjection(len_obj, inactive_threshold, nargout)

                len_obj_min = len_obj_test
                len_obj_max = 1.3
                bi_target = 1
                bi_tolerance = self.algorithm_settings["bisectioning_tolerance"].GetDouble()
                bi_max_itr = self.algorithm_settings["bisectioning_max_itr"].GetInt()
                len_obj_result, bi_itrs, bi_err = PerformBisectioning(func, len_obj_min, len_obj_max, bi_target, bi_tolerance, bi_max_itr)

                nargout = 6
                norm_dX, dX, is_projection_sucessfull, adj_len_obj, adj_len_eqs, adj_len_ineqs = projector.RunProjection(len_obj_result, inactive_threshold, nargout)

            else: # Correction mode
                print ("\n> Computing projection case 2...")

                len_obj = self.algorithm_settings["obj_share_during_correction"].GetDouble()
                func = lambda threshold: projector.RunProjection(len_obj, threshold, nargout)

                threshold_min = 0
                threshold_max = 1.3
                bi_target = 1
                bi_tolerance = self.algorithm_settings["bisectioning_tolerance"].GetDouble()
                bi_max_itr = self.algorithm_settings["bisectioning_max_itr"].GetInt()
                l_threshold_result, bi_itrs, bi_err = PerformBisectioning(func, threshold_min, threshold_max, bi_target, bi_tolerance, bi_max_itr)

                nargout = 6
                norm_dX, dX, is_projection_sucessfull, adj_len_obj, adj_len_eqs, adj_len_ineqs = projector.RunProjection(len_obj, l_threshold_result, nargout)
        else:
            raise RuntimeError("Case of not converged test projection not yet implemented yet!")

        print("\n> Time needed for determining step = ", timer.GetTotalTime(), "s")

        return dX, test_norm_dX, bi_itrs, bi_err, adj_len_obj, adj_len_eqs, adj_len_ineqs

    # --------------------------------------------------------------------------
    def __ComputeShapeUpdate(self, dx_bar, step_length):
        # Compute update in regular units
        dx = ScalarVectorProduct(step_length,dx_bar)

        WriteListToNodalVariable(dx, self.design_surface, SHAPE_UPDATE)
        self.optimization_utilities.AddFirstVariableToSecondVariable(SHAPE_UPDATE, SHAPE_CHANGE)

        return NormInf3D(dx)

    # --------------------------------------------------------------------------
    def __LogCurrentOptimizationStep(self, step_length, len_bar_obj, len_bar_eqs, len_bar_ineqs, test_norm_dx_bar, bi_itrs, bi_err, adj_len_bar_obj, adj_len_bar_eqs, adj_len_bar_ineqs, norm_dx):
        additional_values_to_log = {}
        additional_values_to_log["len_bar_obj"] = len_bar_obj
        additional_values_to_log["adj_len_bar_obj"] = adj_len_bar_obj

        len_bar_cons = self.__CombineConstraintDataToOrderedList(len_bar_eqs, len_bar_ineqs)
        adj_len_bar_cons = self.__CombineConstraintDataToOrderedList(adj_len_bar_eqs, adj_len_bar_ineqs)
        additional_values_to_log["len_bar_cons"] = len_bar_cons
        additional_values_to_log["adj_len_bar_cons"] = adj_len_bar_cons

        additional_values_to_log["test_norm_dx_bar"] = test_norm_dx_bar
        additional_values_to_log["bi_itrs"] = bi_itrs
        additional_values_to_log["bi_err"] = bi_err
        additional_values_to_log["norm_dx"] = norm_dx
        additional_values_to_log["step_length"] = step_length

        self.data_logger.LogCurrentValues(self.opt_iteration, additional_values_to_log)
        self.data_logger.LogCurrentDesign(self.opt_iteration)

    # --------------------------------------------------------------------------
    def __CombineConstraintDataToOrderedList(self, eqs_data_list, ineqs_data_list):
        num_eqs = 0
        num_ineqs = 0
        combined_list = []

        # Order is given by appearance of constraints in optimization settings
        for itr in range(self.specified_constraints.size()):
            if self.specified_constraints[itr]["type"].GetString()=="=":
                combined_list.append(eqs_data_list[num_eqs])
                num_eqs == num_eqs+1
            else:
                combined_list.append(ineqs_data_list[num_ineqs])
                num_ineqs == num_ineqs+1

        return combined_list

# ==============================================================================
class Projector():
    # --------------------------------------------------------------------------
    def __init__(self, len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs, settings):

        self.far_away_length = settings["far_away_length"].GetDouble()
        self.subopt_max_itr = settings["subopt_max_itr"].GetInt()
        self.subopt_tolerance = settings["subopt_tolerance"].GetDouble()

        # Reduce input data to relevant info
        len_ineqs, dir_ineqs =  self.__RemoveInactiveConstraintsFarAway(len_ineqs, dir_ineqs)
        len_eqs, dir_eqs, len_ineqs, dir_ineqs = self.__RemoveConstraintsWithNoGradientInfo(len_eqs, dir_eqs, len_ineqs, dir_ineqs)

        # Store some working variables
        self.len_eqs = len_eqs
        self.len_ineqs = len_ineqs
        self.num_eqs = len(len_eqs)
        self.num_ineqs = len(len_ineqs)
        self.num_unknowns = 1 + self.num_eqs + self.num_ineqs

        # Create orthogonal basis
        self.ortho_basis = self.__PerformGramSchmidtOrthogonalization(dir_obj, dir_eqs, dir_ineqs)

        # Transform directions to orthogonal space since they don't change with different projections
        self.dir_obj_o = TranslateToNewBasis(dir_obj, self.ortho_basis)
        self.dir_eqs_o = TranslateToNewBasis(dir_eqs, self.ortho_basis)
        self.dir_ineqs_o = TranslateToNewBasis(dir_ineqs, self.ortho_basis)

        # Make sure directions of constraints are stored as matrix
        self.dir_eqs_o = SafeConvertVectorToMatrix(self.dir_eqs_o)
        self.dir_ineqs_o = SafeConvertVectorToMatrix(self.dir_ineqs_o)

    # --------------------------------------------------------------------------
    def RunProjection(self, len_obj, threshold, nargout):

        # Adjust halfspaces according input
        len_obj, len_eqs, len_ineqs = self.__AdjustHalfSpacesAndHyperplanes(len_obj, threshold)

        # Determine position of border of halfspaces and hyperplanes
        pos_obj_o, pos_eqs_o, pos_ineqs_o = self.__DetermineConstraintBorders(len_obj, len_eqs, len_ineqs)

        # Project current position onto intersection of
        current_position = ZeroVector(self.num_unknowns)
        dlambda_hp = self.__ProjectToHyperplanes(current_position, self.dir_eqs_o, pos_eqs_o)

        # Project position and direction of halfspaces onto intersection of hyperplanes
        zero_position_eqs_o = ZeroMatrix(self.num_unknowns,self.num_eqs)

        pos_obj_hp = self.__ProjectToHyperplanes(pos_obj_o, HorzCat(self.dir_eqs_o, self.dir_obj_o), HorzCat(pos_eqs_o, pos_obj_o))
        dir_obj_hp = self.__ProjectToHyperplanes(self.dir_obj_o, self.dir_eqs_o, zero_position_eqs_o)

        pos_ineqs_hp = []
        dir_ineqs_hp = []
        for itr in range(self.num_ineqs):
            pos_ineqs_hp_i = self.__ProjectToHyperplanes(pos_ineqs_o[itr], HorzCat(self.dir_eqs_o, self.dir_ineqs_o[itr]), HorzCat(pos_eqs_o, pos_ineqs_o[itr]))
            dir_ineqs_hp_i = self.__ProjectToHyperplanes(self.dir_ineqs_o[itr], self.dir_eqs_o, zero_position_eqs_o)

            pos_ineqs_hp.append(pos_ineqs_hp_i)
            dir_ineqs_hp.append(dir_ineqs_hp_i)

        # Project onto adjusted halfspaces along the intersection of hyperplanes
        dX_o, subopt_itr, error, exit_code = self.__ProjectToHalfSpaces(dlambda_hp, HorzCat(pos_ineqs_hp, pos_obj_hp), HorzCat(dir_ineqs_hp, dir_obj_hp))

        # Determine return values
        if exit_code == 0:
            is_projection_sucessfull = True

            # Backtransformation and multiplication with -1 because the direction vectors are chosen opposite to the gradients such that the lengths are positive if violated
            dX = ScalarVectorProduct(-1, TranslateToOriginalBasis(dX_o, self.ortho_basis))
            norm_dX = NormInf3D(dX)
        else:
            is_projection_sucessfull = False

            dX = []
            norm_dX = 1e10

        # Allow to specify output arguments to allow a more general definition of algorithms utilizing this function
        if nargout == 2:
            return norm_dX, is_projection_sucessfull
        elif nargout == 3:
            return norm_dX, dX, is_projection_sucessfull
        else:
            return norm_dX, dX, is_projection_sucessfull, len_obj, len_eqs, len_ineqs

    # --------------------------------------------------------------------------
    def __RemoveInactiveConstraintsFarAway(self, len_ineqs, dir_ineqs):
        len_ineqs_relevant = []
        dir_ineqs_relevant = []

        for itr in range(len(len_ineqs)):
            len_i = len_ineqs[itr]
            dir_i = dir_ineqs[itr]

            if len_i > -self.far_away_length:
                len_ineqs_relevant.append(len_i)
                dir_ineqs_relevant.append(dir_i)

        return len_ineqs_relevant, dir_ineqs_relevant

    # --------------------------------------------------------------------------
    def __RemoveConstraintsWithNoGradientInfo(self, len_eqs, dir_eqs, len_ineqs, dir_ineqs):
        len_eqs_relevant = []
        dir_eqs_relevant = []

        for itr in range(len(dir_eqs)):
            len_i = len_eqs[itr]
            dir_i = dir_eqs[itr]

            if NormInf3D(dir_i) > 1e-13:
                len_eqs_relevant.append(len_i)
                dir_eqs_relevant.append(dir_i)

        len_ineqs_relevant = []
        dir_ineqs_relevant = []

        for itr in range(len(dir_ineqs)):
            len_i = len_ineqs[itr]
            dir_i = dir_ineqs[itr]

            if NormInf3D(dir_i) > 1e-13:
                len_ineqs_relevant.append(len_i)
                dir_ineqs_relevant.append(dir_i)

        return len_eqs_relevant, dir_eqs_relevant, len_ineqs_relevant, dir_ineqs_relevant

    # --------------------------------------------------------------------------
    def __PerformGramSchmidtOrthogonalization(self, dir_obj, dir_eqs, dir_ineqs):
        # Determine vector space
        V = [dir_obj]
        for itr in range(len(dir_eqs)):
            V = HorzCat(V, dir_eqs[itr])
        for itr in range(len(dir_ineqs)):
            V = HorzCat(V, dir_ineqs[itr])

        B = []

        # Orthogonalization
        norm2_V0 = Norm2(V[0])
        B.append( ScalarVectorProduct(1/norm2_V0,V[0]) )
        for v in V[1:]:
            for b in B:
                norm2_b = Norm2(b)
                v = Minus( v , ScalarVectorProduct( Dot(v,b)/norm2_b**2 , b ) )

            # Add only if vector is independent
            norm2_v = Norm2(v)
            if norm2_v>1e-10:
                B.append( ScalarVectorProduct(1/norm2_v,v) )
            else:
                print("Zero basis vector after Gram-Schmidt orthogonalization!")
                B.append(v)

        return B

    # --------------------------------------------------------------------------
    def __AdjustHalfSpacesAndHyperplanes(self, len_obj, threshold):
        if threshold<len_obj:
            len_obj = threshold

        len_eqs = copy.deepcopy(self.len_eqs)
        len_ineqs = copy.deepcopy(self.len_ineqs)

        for itr in range(self.num_eqs):
            len_eqs[itr] = min(max(self.len_eqs[itr],-threshold),threshold)

        for itr in range(self.num_ineqs):
            len_ineqs[itr] = min(self.len_ineqs[itr],threshold)

        return len_obj, len_eqs, len_ineqs

    # --------------------------------------------------------------------------
    def __DetermineConstraintBorders(self, len_obj, len_eqs, len_ineqs):
        pos_obj = ScalarVectorProduct(-len_obj,self.dir_obj_o)

        pos_eqs = []
        pos_ineqs = []
        for i in range(self.num_eqs):
            pos_eqs.append(ScalarVectorProduct(-len_eqs[i],self.dir_eqs_o[i]))

        for i in range(self.num_ineqs):
            pos_ineqs.append(ScalarVectorProduct(-len_ineqs[i],self.dir_ineqs_o[i]))

        return pos_obj, pos_eqs, pos_ineqs

    # --------------------------------------------------------------------------
    def __ProjectToHyperplanes(self, vector, dir_hps, pos_hps):
            if IsEmpty(dir_hps):
                return vector

            num_hps = len(dir_hps)

            tmp_mat = Prod(Trans(dir_hps),dir_hps)
            tmp_vec = [ Dot(dir_hps[j],Minus(pos_hps[j],vector)) for j in range(num_hps) ]

            tmp_solution = SolveLinearSystem(tmp_mat,tmp_vec)

            return Plus(Prod(dir_hps,tmp_solution),vector)

    # --------------------------------------------------------------------------
    def __ProjectToHalfSpaces(self, dx0, pos_hss, dir_hss):
        A = Trans(dir_hss)
        b = [ Dot(pos_hss[i],dir_hss[i]) for i in range(RowSize(A)) ]

        dX_o, subopt_itr, error, exit_code = QuadProg(A, b, self.subopt_max_itr, self.subopt_tolerance)

        # Consider initial delta
        dX_o = Plus(dX_o,dx0)

        return dX_o, subopt_itr, error, exit_code

# ==============================================================================