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
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
from .algorithm_base import OptimizationAlgorithm
from . import mapper_factory
from . import data_logger_factory
from . import custom_math as cm
from .custom_timer import Timer
from .custom_variable_utilities import WriteDictionaryDataOnNodalVariable, ReadNodalVariableToList, WriteListToNodalVariable
import copy

# ==============================================================================
class AlgorithmTrustRegion(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = KM.Parameters("""
        {
            "name"                          : "trust_region",
            "max_step_length"               : 1.0,
            "step_length_tolerance"         : 1e-3,
            "step_length_reduction_factor"  : 0.5,
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

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Trust-region algorithm only supports one objective function!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.model_part_controller.Initialize()

        self.analyzer.InitializeBeforeOptimizationLoop()

        self.design_surface = self.model_part_controller.GetDesignSurface()

        self.mapper = mapper_factory.CreateMapper(self.design_surface, self.design_surface, self.mapper_settings)
        self.mapper.Initialize()

        self.data_logger = data_logger_factory.CreateDataLogger(self.model_part_controller, self.communicator, self.optimization_settings)
        self.data_logger.InitializeDataLogging()

        self.optimization_utilities = KSO.OptimizationUtilities(self.design_surface, self.optimization_settings)

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        for self.opt_iteration in range(1,self.algorithm_settings["max_iterations"].GetInt()+1):
            KM.Logger.Print("")
            KM.Logger.Print("===============================================================================")
            KM.Logger.PrintInfo("ShapeOpt", timer.GetTimeStamp(), ": Starting optimization iteration ",self.opt_iteration)
            KM.Logger.Print("===============================================================================\n")

            timer.StartNewLap()

            self.__InitializeNewShape()

            self.__AnalyzeShape()

            self.__PostProcessGradientsObtainedFromAnalysis()

            len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs =  self.__ConvertAnalysisResultsToLengthDirectionFormat()

            step_length = self.__DetermineMaxStepLength()

            len_bar_obj, len_bar_eqs, len_bar_ineqs = self.__ExpressInStepLengthUnit(len_obj, len_eqs, len_ineqs, step_length)

            dX_bar, process_details = self.__DetermineStep(len_bar_obj, dir_obj, len_bar_eqs, dir_eqs, len_bar_ineqs, dir_ineqs)

            dX = self.__ComputeShapeUpdate(dX_bar, step_length)

            values_to_be_logged = {}
            values_to_be_logged["len_bar_obj"] = len_bar_obj
            values_to_be_logged["len_bar_cons"] = self.__CombineConstraintDataToOrderedList(len_bar_eqs, len_bar_ineqs)
            values_to_be_logged["step_length"] = step_length
            values_to_be_logged["test_norm_dX_bar"] = process_details["test_norm_dX"]
            values_to_be_logged["bi_itrs"] = process_details["bi_itrs"]
            values_to_be_logged["bi_err"] = process_details["bi_err"]
            values_to_be_logged["adj_len_bar_obj"] = process_details["adj_len_obj"]
            values_to_be_logged["adj_len_bar_cons"] = self.__CombineConstraintDataToOrderedList(process_details["adj_len_eqs"], process_details["adj_len_ineqs"])
            values_to_be_logged["norm_dX"] = cm.NormInf3D(dX)

            self.__LogCurrentOptimizationStep(values_to_be_logged)

            KM.Logger.Print("")
            KM.Logger.PrintInfo("ShapeOpt", "Time needed for current optimization step = ", timer.GetLapTime(), "s")
            KM.Logger.PrintInfo("ShapeOpt", "Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.analyzer.FinalizeAfterOptimizationLoop()
        self.data_logger.FinalizeDataLogging()

    # --------------------------------------------------------------------------
    def __InitializeNewShape(self):
        self.model_part_controller.UpdateTimeStep(self.opt_iteration)
        self.model_part_controller.UpdateMeshAccordingInputVariable(KSO.SHAPE_UPDATE)
        self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __AnalyzeShape(self):
            self.communicator.initializeCommunication()

            obj_id = self.objectives[0]["identifier"].GetString()
            self.communicator.requestValueOf(obj_id)
            self.communicator.requestGradientOf(obj_id)

            for itr in range(self.constraints.size()):
                con_id =  self.constraints[itr]["identifier"].GetString()
                self.communicator.requestValueOf(con_id)
                self.communicator.requestGradientOf(con_id)

            self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.opt_iteration, self.communicator)

    # --------------------------------------------------------------------------
    def __PostProcessGradientsObtainedFromAnalysis(self):
        # Compute surface normals if required
        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ComputeUnitSurfaceNormals()
        else:
            for itr in range(self.constraints.size()):
                if self.constraints[itr]["project_gradient_on_surface_normals"].GetBool():
                    self.model_part_controller.ComputeUnitSurfaceNormals()

        # Process objective gradients
        obj = self.objectives[0]
        obj_id = obj["identifier"].GetString()

        obj_gradients_dict = self.communicator.getStandardizedGradient(obj_id)

        nodal_variable = KM.KratosGlobals.GetVariable("DF1DX")
        WriteDictionaryDataOnNodalVariable(obj_gradients_dict, self.optimization_model_part, nodal_variable)

        # Projection on surface normals
        if obj["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(nodal_variable)

        # Damping
        self.model_part_controller.DampNodalVariableIfSpecified(nodal_variable)

        # Mapping
        nodal_variable_mapped = KM.KratosGlobals.GetVariable("DF1DX_MAPPED")
        self.mapper.Update()
        self.mapper.InverseMap(nodal_variable, nodal_variable_mapped)
        self.mapper.Map(nodal_variable_mapped, nodal_variable_mapped)

        # Damping
        self.model_part_controller.DampNodalVariableIfSpecified(nodal_variable_mapped)

        # Process constraint gradients
        for itr in range(self.constraints.size()):
            con = self.constraints[itr]
            con_id = con["identifier"].GetString()

            eq_gradients_dict = self.communicator.getStandardizedGradient(con_id)

            nodal_variable = KM.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX")
            WriteDictionaryDataOnNodalVariable(eq_gradients_dict, self.optimization_model_part, nodal_variable)

            # Projection on surface normals
            if con["project_gradient_on_surface_normals"].GetBool():
                self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(nodal_variable)

            # Damping
            self.model_part_controller.DampNodalVariableIfSpecified(nodal_variable)

            # Mapping
            nodal_variable_mapped = KM.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")
            self.mapper.InverseMap(nodal_variable, nodal_variable_mapped)
            self.mapper.Map(nodal_variable_mapped, nodal_variable_mapped)

            # Damping
            self.model_part_controller.DampNodalVariableIfSpecified(nodal_variable_mapped)

    # --------------------------------------------------------------------------
    def __ConvertAnalysisResultsToLengthDirectionFormat(self):
        # Convert objective results
        obj = self.objectives[0]
        obj_id = obj["identifier"].GetString()

        nodal_variable = KM.KratosGlobals.GetVariable("DF1DX")
        nodal_variable_mapped = KM.KratosGlobals.GetVariable("DF1DX_MAPPED")

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

        for itr in range(self.constraints.size()):
            con = self.constraints[itr]
            con_id = con["identifier"].GetString()

            nodal_variable = KM.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX")
            nodal_variable_mapped = KM.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")

            value = self.communicator.getStandardizedValue(con_id)
            gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
            gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

            direction, length = self.__ConvertToLengthDirectionFormat(value, gradient, gradient_mapped)

            if con["type"].GetString()=="=":
                dir_eqs.append(direction)
                len_eqs.append(length)
            else:
                dir_ineqs.append(direction)
                len_ineqs.append(length)

        return len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs

    # --------------------------------------------------------------------------
    @staticmethod
    def __ConvertToLengthDirectionFormat(value, gradient, modified_gradient):
        norm_inf = cm.NormInf3D(modified_gradient)
        if norm_inf > 1e-12:
            direction = cm.ScalarVectorProduct(-1/norm_inf,modified_gradient)
            length = -value/cm.Dot(gradient, direction)
        else:
            KM.Logger.PrintWarning("ShapeOpt::AlgorithmTrustRegion", "Vanishing norm-infinity for gradient detected!")
            direction = modified_gradient
            length = 0.0

        return direction, length

    # --------------------------------------------------------------------------
    def __DetermineMaxStepLength(self):
        if self.opt_iteration < 4:
            return self.algorithm_settings["max_step_length"].GetDouble()
        else:
            obj_id = self.objectives[0]["identifier"].GetString()
            current_obj_val = self.communicator.getStandardizedValue(obj_id)
            obj_history = self.data_logger.GetValues("response_value")[obj_id]
            step_history = self.data_logger.GetValues("step_length")

            # Check for osciallation
            objective_is_oscillating = False
            is_decrease_1 = (current_obj_val - obj_history[self.opt_iteration-1])< 0
            is_decrease_2 = (obj_history[self.opt_iteration-1] - obj_history[self.opt_iteration-2])<0
            is_decrease_3 = (current_obj_val - obj_history[self.opt_iteration-3])< 0
            if (is_decrease_1 and is_decrease_2== False and is_decrease_3) or (is_decrease_1== False and is_decrease_2 and is_decrease_3==False):
                objective_is_oscillating = True

            # Reduce step length if certain conditions are fullfilled
            if objective_is_oscillating:
                return step_history[self.opt_iteration-1]*self.algorithm_settings["step_length_reduction_factor"].GetDouble()
            else:
                return step_history[self.opt_iteration-1]

    # --------------------------------------------------------------------------
    @staticmethod
    def __ExpressInStepLengthUnit(len_obj, len_eqs, len_ineqs, step_length):
        len_bar_obj = 1/step_length * len_obj
        len_bar_eqs = cm.ScalarVectorProduct(1/step_length, len_eqs)
        len_bar_ineqs = cm.ScalarVectorProduct(1/step_length, len_ineqs)
        return len_bar_obj, len_bar_eqs, len_bar_ineqs

    # --------------------------------------------------------------------------
    def __DetermineStep(self, len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs):
        KM.Logger.Print("")
        KM.Logger.PrintInfo("ShapeOpt", "Starting determination of step...")

        timer = Timer()
        timer.StartTimer()

        # Create projector object wich can do the projection in the orthogonalized subspace
        projector = Projector(len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs, self.algorithm_settings)

        # 1. Test projection if there is room for objective improvement
        # I.e., the actual step length to become feasible for an inactive threshold is smaller than 1 and hence a part of the step can be dedicated to objective improvement
        len_obj_test = 0.01
        inactive_threshold = 100
        test_norm_dX, is_projection_sucessfull = projector.RunProjection(len_obj_test, inactive_threshold)

        KM.Logger.PrintInfo("ShapeOpt", "Time needed for one projection step = ", timer.GetTotalTime(), "s")

        # 2. Determine step following two different modes depending on the previos found step length to the feasible domain
        if is_projection_sucessfull:
            if test_norm_dX < 1: # Minimizing mode
                print ("\n> Computing projection case 1...")

                func = lambda len_obj: projector.RunProjection(len_obj, inactive_threshold)

                len_obj_min = len_obj_test
                len_obj_max = 1.3
                bi_target = 1
                bi_tolerance = self.algorithm_settings["bisectioning_tolerance"].GetDouble()
                bi_max_itr = self.algorithm_settings["bisectioning_max_itr"].GetInt()
                len_obj_result, bi_itrs, bi_err = cm.PerformBisectioning(func, len_obj_min, len_obj_max, bi_target, bi_tolerance, bi_max_itr)

                projection_results = projector.GetDetailedResultsOfLatestProjection()

            else: # Correction mode
                print ("\n> Computing projection case 2...")

                len_obj = self.algorithm_settings["obj_share_during_correction"].GetDouble()
                func = lambda threshold: projector.RunProjection(len_obj, threshold)

                threshold_min = 0
                threshold_max = 1.3
                bi_target = 1
                bi_tolerance = self.algorithm_settings["bisectioning_tolerance"].GetDouble()
                bi_max_itr = self.algorithm_settings["bisectioning_max_itr"].GetInt()
                l_threshold_result, bi_itrs, bi_err = cm.PerformBisectioning(func, threshold_min, threshold_max, bi_target, bi_tolerance, bi_max_itr)

                projection_results = projector.GetDetailedResultsOfLatestProjection()
        else:
            raise RuntimeError("Case of not converged test projection not yet implemented yet!")

        KM.Logger.Print("")
        KM.Logger.PrintInfo("ShapeOpt", "Time needed for determining step = ", timer.GetTotalTime(), "s")

        process_details = { "test_norm_dX": test_norm_dX,
                            "bi_itrs":bi_itrs,
                            "bi_err":bi_err,
                            "adj_len_obj": projection_results["adj_len_obj"],
                            "adj_len_eqs": projection_results["adj_len_eqs"],
                            "adj_len_ineqs": projection_results["adj_len_ineqs"] }

        return projection_results["dX"], process_details

    # --------------------------------------------------------------------------
    def __ComputeShapeUpdate(self, dX_bar, step_length):
        # Compute update in regular units
        dX = cm.ScalarVectorProduct(step_length,dX_bar)

        WriteListToNodalVariable(dX, self.design_surface, KSO.SHAPE_UPDATE)
        self.optimization_utilities.AddFirstVariableToSecondVariable(KSO.SHAPE_UPDATE, KSO.SHAPE_CHANGE)

        return dX

    # --------------------------------------------------------------------------
    def __LogCurrentOptimizationStep(self, additional_values_to_log):
        self.data_logger.LogCurrentValues(self.opt_iteration, additional_values_to_log)
        self.data_logger.LogCurrentDesign(self.opt_iteration)

    # --------------------------------------------------------------------------
    def __CombineConstraintDataToOrderedList(self, eqs_data_list, ineqs_data_list):
        num_eqs = 0
        num_ineqs = 0
        combined_list = []

        # Order is given by appearance of constraints in optimization settings
        for itr in range(self.constraints.size()):
            if self.constraints[itr]["type"].GetString()=="=":
                combined_list.append(eqs_data_list[num_eqs])
                num_eqs = num_eqs+1
            else:
                combined_list.append(ineqs_data_list[num_ineqs])
                num_ineqs = num_ineqs+1

        return combined_list

# ==============================================================================
class Projector():
    # --------------------------------------------------------------------------
    def __init__(self, len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs, settings):

        # Store settings
        self.far_away_length = settings["far_away_length"].GetDouble()
        self.subopt_max_itr = settings["subopt_max_itr"].GetInt()
        self.subopt_tolerance = settings["subopt_tolerance"].GetDouble()

        # Initialize projection results
        self.are_projection_restuls_stored = False
        self.projection_results = {}

        # Reduce input data to relevant info
        self.input_len_obj = len_obj
        self.input_len_eqs = len_eqs
        self.input_len_ineqs = len_ineqs

        len_eqs, dir_eqs, remaining_eqs_entries = self.__ReduceToRelevantEqualityConstraints(len_eqs, dir_eqs)
        len_ineqs, dir_ineqs, remaining_ineqs_entries = self.__ReduceToRelevantInequalityConstraints(len_ineqs, dir_ineqs)

        # Store some working variables depening on data reduction
        self.remaining_eqs_entries = remaining_eqs_entries
        self.remaining_ineqs_entries = remaining_ineqs_entries
        self.len_eqs = len_eqs
        self.len_ineqs = len_ineqs
        self.num_eqs = len(len_eqs)
        self.num_ineqs = len(len_ineqs)
        self.num_unknowns = 1 + self.num_eqs + self.num_ineqs

        # Create orthogonal basis
        vector_space = [dir_obj]
        for itr in range(len(dir_eqs)):
            vector_space = cm.HorzCat(vector_space, dir_eqs[itr])
        for itr in range(len(dir_ineqs)):
            vector_space = cm.HorzCat(vector_space, dir_ineqs[itr])
        self.ortho_basis = cm.PerformGramSchmidtOrthogonalization(vector_space)

        # Transform directions to orthogonal space since they don't change with different projections
        self.dir_obj_o = cm.TranslateToNewBasis(dir_obj, self.ortho_basis)
        self.dir_eqs_o = cm.TranslateToNewBasis(dir_eqs, self.ortho_basis)
        self.dir_ineqs_o = cm.TranslateToNewBasis(dir_ineqs, self.ortho_basis)

        # Make sure directions of constraints are stored as matrix
        self.dir_eqs_o = cm.SafeConvertVectorToMatrix(self.dir_eqs_o)
        self.dir_ineqs_o = cm.SafeConvertVectorToMatrix(self.dir_ineqs_o)

    # --------------------------------------------------------------------------
    def RunProjection(self, len_obj, threshold):

        # Adjust halfspaces according input
        adj_len_obj, adj_len_eqs, adj_len_ineqs = self.__AdjustHalfSpacesAndHyperplanes(len_obj, threshold)

        # Determine position of border of halfspaces and hyperplanes
        pos_obj_o, pos_eqs_o, pos_ineqs_o = self.__DetermineConstraintBorders(adj_len_obj, adj_len_eqs, adj_len_ineqs)

        # Project current position onto intersection of
        current_position = cm.ZeroVector(self.num_unknowns)
        dlambda_hp = self.__ProjectToHyperplanes(current_position, self.dir_eqs_o, pos_eqs_o)

        # Project position and direction of halfspaces onto intersection of hyperplanes
        zero_position_eqs_o = cm.ZeroMatrix(self.num_unknowns,self.num_eqs)

        pos_obj_hp = self.__ProjectToHyperplanes(pos_obj_o, cm.HorzCat(self.dir_eqs_o, self.dir_obj_o), cm.HorzCat(pos_eqs_o, pos_obj_o))
        dir_obj_hp = self.__ProjectToHyperplanes(self.dir_obj_o, self.dir_eqs_o, zero_position_eqs_o)

        pos_ineqs_hp = []
        dir_ineqs_hp = []
        for itr in range(self.num_ineqs):
            pos_ineqs_hp_i = self.__ProjectToHyperplanes(pos_ineqs_o[itr], cm.HorzCat(self.dir_eqs_o, self.dir_ineqs_o[itr]), cm.HorzCat(pos_eqs_o, pos_ineqs_o[itr]))
            dir_ineqs_hp_i = self.__ProjectToHyperplanes(self.dir_ineqs_o[itr], self.dir_eqs_o, zero_position_eqs_o)

            pos_ineqs_hp.append(pos_ineqs_hp_i)
            dir_ineqs_hp.append(dir_ineqs_hp_i)

        # Project onto adjusted halfspaces along the intersection of hyperplanes
        dX_o, _, _, exit_code = self.__ProjectToHalfSpaces(dlambda_hp, cm.HorzCat(pos_ineqs_hp, pos_obj_hp), cm.HorzCat(dir_ineqs_hp, dir_obj_hp))

        # Determine return values
        if exit_code == 0:
            is_projection_sucessfull = True

            # Backtransformation and multiplication with -1 because the direction vectors are chosen opposite to the gradients such that the lengths are positive if violated
            dX = cm.ScalarVectorProduct(-1, cm.TranslateToOriginalBasis(dX_o, self.ortho_basis))
            norm_dX = cm.NormInf3D(dX)
        else:
            is_projection_sucessfull = False

            dX = []
            norm_dX = 1e10

        self.__StoreProjectionResults(norm_dX, dX, is_projection_sucessfull, adj_len_obj, adj_len_eqs, adj_len_ineqs)

        return norm_dX, is_projection_sucessfull


    # --------------------------------------------------------------------------
    def GetDetailedResultsOfLatestProjection(self):
        if self.are_projection_restuls_stored == False:
            raise RuntimeError("Projector::__StoreProjectionResults: No projection results stored yet!")

        return self.projection_results

    # --------------------------------------------------------------------------
    @staticmethod
    def __ReduceToRelevantEqualityConstraints(len_eqs, dir_eqs):
        len_eqs_relevant = []
        dir_eqs_relevant = []
        remaining_entries = []

        for itr in range(len(dir_eqs)):
            len_i = len_eqs[itr]
            dir_i = dir_eqs[itr]

            is_no_gradient_info_available = cm.NormInf3D(dir_i) < 1e-13

            if is_no_gradient_info_available:
                pass
            else:
                remaining_entries.append(itr)
                len_eqs_relevant.append(len_i)
                dir_eqs_relevant.append(dir_i)

        return len_eqs_relevant, dir_eqs_relevant, remaining_entries

    # --------------------------------------------------------------------------
    def __ReduceToRelevantInequalityConstraints(self, len_ineqs, dir_ineqs):
        len_ineqs_relevant = []
        dir_ineqs_relevant = []
        remaining_entries = []

        for itr in range(len(dir_ineqs)):
            len_i = len_ineqs[itr]
            dir_i = dir_ineqs[itr]

            is_no_gradient_info_available = cm.NormInf3D(dir_i) < 1e-13
            is_constraint_inactive_and_far_away = len_i < -self.far_away_length

            if is_no_gradient_info_available or is_constraint_inactive_and_far_away:
                pass
            else:
                remaining_entries.append(itr)
                len_ineqs_relevant.append(len_i)
                dir_ineqs_relevant.append(dir_i)

        return len_ineqs_relevant, dir_ineqs_relevant, remaining_entries

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
        pos_obj = cm.ScalarVectorProduct(-len_obj,self.dir_obj_o)

        pos_eqs = []
        pos_ineqs = []
        for i in range(self.num_eqs):
            pos_eqs.append(cm.ScalarVectorProduct(-len_eqs[i],self.dir_eqs_o[i]))

        for i in range(self.num_ineqs):
            pos_ineqs.append(cm.ScalarVectorProduct(-len_ineqs[i],self.dir_ineqs_o[i]))

        return pos_obj, pos_eqs, pos_ineqs

    # --------------------------------------------------------------------------
    @staticmethod
    def __ProjectToHyperplanes(vector, dir_hps, pos_hps):
        if cm.IsEmpty(dir_hps):
            return vector

        num_hps = len(dir_hps)

        tmp_mat = cm.Prod(cm.Trans(dir_hps),dir_hps)
        tmp_vec = [ cm.Dot(dir_hps[j],cm.Minus(pos_hps[j],vector)) for j in range(num_hps) ]

        tmp_solution = cm.SolveLinearSystem(tmp_mat,tmp_vec)

        return cm.Plus(cm.Prod(dir_hps,tmp_solution),vector)

    # --------------------------------------------------------------------------
    def __ProjectToHalfSpaces(self, dX0, pos_hss, dir_hss):
        A = cm.Trans(dir_hss)
        b = [ cm.Dot(pos_hss[i],dir_hss[i]) for i in range(cm.RowSize(A)) ]

        dX_o, subopt_itr, error, exit_code = cm.QuadProg(A, b, self.subopt_max_itr, self.subopt_tolerance)

        # Consider initial delta
        dX_o = cm.Plus(dX_o,dX0)

        return dX_o, subopt_itr, error, exit_code

    # --------------------------------------------------------------------------
    def __StoreProjectionResults(self, norm_dX, dX, is_projection_sucessfull, adj_len_obj, adj_len_eqs, adj_len_ineqs):
        self.are_projection_restuls_stored = True
        self.projection_results["norm_dX"] = norm_dX
        self.projection_results["dX"] = dX
        self.projection_results["is_projection_sucessfull"] = is_projection_sucessfull
        self.projection_results["adj_len_obj"] = adj_len_obj
        self.projection_results["adj_len_eqs"], self.projection_results["adj_len_ineqs"] = self.__CompleteConstraintLengthsWithRemovedEntries(adj_len_eqs, adj_len_ineqs)

    # --------------------------------------------------------------------------
    def __CompleteConstraintLengthsWithRemovedEntries(self, len_eqs, len_ineqs):
        # Complete list of eqs
        complete_list_eqs = copy.deepcopy(self.input_len_eqs)
        for itr in range(len(len_eqs)):
            original_eq_number = self.remaining_eqs_entries[itr]
            complete_list_eqs[original_eq_number] = len_eqs[itr]

        # Complete list of ineqs
        complete_list_ineqs = copy.deepcopy(self.input_len_ineqs)
        for itr in range(len(len_ineqs)):
            original_eq_number = self.remaining_ineqs_entries[itr]
            complete_list_ineqs[original_eq_number] = len_ineqs[itr]

        return complete_list_eqs, complete_list_ineqs

# ==============================================================================