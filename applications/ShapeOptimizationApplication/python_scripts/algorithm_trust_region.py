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
from custom_math import NormInf3D, Dot, ScalarVectorProduct, Norm2, RowSize, CollSize, HorzCat, Minus, TranslateToNewBasis, TranslateToOriginalBasis, Trans, QuadProg, Prod, PerformBisectioning
from custom_variable_utilities import WriteDictionaryDataOnNodalVariable, ReadNodalVariableToList, WriteNodeCoordinatesToList, WriteListToNodalVariable
from custom_timer import Timer
import mapper_factory

# ==============================================================================
class AlgorithmTrustRegion(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = Parameters("""
        {
            "name"                          : "trust_region",
            "max_step_length"               : 1.0,
            "step_length_tolerance"         : 1e-3,
            "step_length_decrease_factor"   : 2.0,
            "step_length_increase_factor"   : 1.2,
            "min_share_objective"           : 0.1,
            "max_iterations"                : 10,
            "far_away_length"               : 2.0,
            "subopt_max_itr"                : 50,
            "subopt_tolerance"              : 1e-10,
            "bisectioning_max_itr"          : 30,
            "bisectioning_tolerance"        : 1e-4
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

        self.GeometryUtilities = GeometryUtilities(self.design_surface)

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.specified_objectives.size() > 1:
            raise RuntimeError("Trust-region algorithm only supports one objective function!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
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

            len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs =  self.__ConvertAnalysisResultsToLengthDirectionFormat()

            len_ineqs, dir_ineqs =  self.__ReduceToRelevantInequalityConstraints(len_ineqs, dir_ineqs)

            # Apply step length control
            step_length = self.algorithm_settings["max_step_length"].GetDouble()

            # Express lengths in step length unit (indicated by the term "bar")
            len_bar_obj = 1/step_length * len_obj
            len_bar_eqs = ScalarVectorProduct(1/step_length, len_eqs)
            len_bar_ineqs = ScalarVectorProduct(1/step_length, len_ineqs)

            # Compute step direction in step length units
            x_bar = self.__DetermineStep(len_obj, dir_obj, len_bar_eqs, dir_eqs, len_bar_ineqs, dir_ineqs)

            # Compute actual step
            dx = [step_length*x_bar[i] for i in range(self.n)]

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

            # Reset possible shape modifications during analysis
            self.model_part_controller.SetMeshToReferenceMesh()
            self.model_part_controller.SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    def __ProcessAnalysisResults(self):
        # Process objective results
        obj = self.specified_objectives[0]
        obj_id = obj["identifier"].GetString()

        obj_value = self.communicator.getStandardizedValue(obj_id)
        obj_gradients_dict = self.communicator.getStandardizedGradient(obj_id)

        nodal_variable = KratosGlobals.GetVariable("DF1DX")
        WriteDictionaryDataOnNodalVariable(obj_gradients_dict, self.optimization_model_part, nodal_variable)

        if obj["project_gradient_on_surface_normals"].GetBool():
            self.__PerformProjectionOnNormals(nodal_variable)

        nodal_variable_mapped = KratosGlobals.GetVariable("DF1DX_MAPPED")
        self.__PerformMapping(nodal_variable, nodal_variable_mapped)

        # Process constraint results
        eq_values = []
        ineq_values = []

        for itr in range(self.specified_constraints.size()):
            con = self.specified_constraints[itr]

            # Process equality constraints
            if con["type"].GetString()=="=":
                eq_id = con["identifier"].GetString()

                eq_values.append(self.communicator.getStandardizedValue(eq_id))
                eq_gradients_dict = self.communicator.getStandardizedGradient(eq_id)

                nodal_variable = KratosGlobals.GetVariable("DC"+str(itr+1)+"DX")
                WriteDictionaryDataOnNodalVariable(eq_gradients_dict, self.optimization_model_part, nodal_variable)

                if con["project_gradient_on_surface_normals"].GetBool():
                    self.__PerformProjectionOnNormals(nodal_variable)

                nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")
                self.__PerformMapping(nodal_variable, nodal_variable_mapped)

            # Process inequality constraints
            else:
                ineq_id = con["identifier"].GetString()

                ineq_values.append(self.communicator.getStandardizedValue(ineq_id))
                ineq_gradients_dict = self.communicator.getStandardizedGradient(ineq_id)

                nodal_variable = KratosGlobals.GetVariable("DC"+str(itr+1)+"DX")
                WriteDictionaryDataOnNodalVariable(ineq_gradients_dict, self.optimization_model_part, nodal_variable)

                if con["project_gradient_on_surface_normals"].GetBool():
                    self.__PerformProjectionOnNormals(nodal_variable)

                nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")
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
        if norm_inf < 1e-12:
            raise RuntimeError("Vanishing norm-infinity detected!")

        direction = [-modified_gradient[itr]/norm_inf for itr in range(len(modified_gradient))]

        if value == None:
            return direction
        else:
            grad_dot_dir = Dot(gradient, direction)
            length = -value/grad_dot_dir
            return direction, length

    # --------------------------------------------------------------------------
    def __ReduceToRelevantInequalityConstraints(self, len_ineqs, dir_ineqs):
        len_ineqs_relevant = []
        dir_ineqs_relevant = []

        for itr in range(len(len_ineqs)):
            len_i = len_ineqs[itr]
            dir_i = dir_ineqs[itr]

            if len_i < self.algorithm_settings["far_away_length"].GetDouble():
                len_ineqs_relevant.append(len_i)
                dir_ineqs_relevant.append(dir_i)

        return len_ineqs_relevant, dir_ineqs_relevant

    # --------------------------------------------------------------------------
    def __DetermineStep(self, len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs):

        timer = Timer()
        timer.StartTimer()


        # Create projector object wich can do the projection in the orthogonalized subspace
        projector = Projector(len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs)

        # 1. Test projection if there is room for objective improvement
        # I.e., the actual step length to become feasible for an inactive threshold is smaller than 1 and hence a part of the step can be dedicated to objective improvement
        len_obj_test = 0.01
        inactive_threshold = 100
        nargout = 2

        norm_inf_dX, is_projection_sucessfull = projector.RunProjection(len_obj, inactive_threshold, nargout)

        print("\n> Time needed for one projection step = ", timer.GetTotalTime(), "s")

        # 2. Determine step following two different modes depending on the previos found step length to the feasible domain

        bi_tolerance = self.algorithm_settings["bisectioning_tolerance"].GetDouble()
        bi_max_itr = self.algorithm_settings["bisectioning_max_itr"].GetInt()
        bi_target = 1

        if is_projection_sucessfull:
            if norm_inf_dX < 1: # Minimizing mode
                print ("> Computing projection case 1...")

                len_obj_min = len_obj_test
                len_obj_max = 1.1

                func = lambda len_obj: projector.RunProjection(len_obj, inactive_threshold, nargout)

                len_obj_result, bi_itrs, bi_err = PerformBisectioning(func, len_obj_min, len_obj_max, bi_target, bi_tolerance, bi_max_itr)

                print(bi_itrs)
                print(bi_err)

        #     else: # Correction mode


        # else:









        # # Corrective behavior --> choose only a minimal objective share and adjust the constraint shares within the limits of the trust region (3d infinity norm of dx = 1)
        # if norm_inf_dx>1:
        #     func = *NormInf3D(*self.__ApplyProjection(-10, dir_obj, len_share_eqs, dir_eqs, len_share_ineqs, threshold))

        # # Minimizing behavior --> increase the objective share as much as possible within the limits of the trust region (3d infinity norm of dx = 1)
        # elif norm_inf_dx<1:
        #     pss

        # # No room for adjustments
        # else:
        #     pass

        print("\n> Time needed for determining step = ", timer.GetTotalTime(), "s")

        return dX

    # --------------------------------------------------------------------------
    def __UpdateMesh(self, dx):
        WriteListToNodalVariable(dx, self.design_surface, SHAPE_UPDATE)
        self.model_part_controller.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
        self.model_part_controller.SetReferenceMeshToMesh()

        x = WriteNodeCoordinatesToList(self.design_surface)
        dxAbsolute = [x[i]+dx[i] - self.x_init[i] for i in range(self.number_of_design_variables)]
        WriteListToNodalVariable(dxAbsolute, self.design_surface, SHAPE_CHANGE)

# ==============================================================================
class Projector():
    # --------------------------------------------------------------------------
    def __init__(self, len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs):

        # Store some working variables
        self.num_eqs = len(len_eqs)
        self.num_ineqs = len(len_ineqs)

        # Store lengths with which initial position of constraint border is computed as reference
        self.reference_len_obj = len_obj
        self.reference_len_eqs = len_eqs
        self.reference_len_ineqs = len_ineqs

        # Create orthogonal basis
        self.ortho_basis = self.__PerformGramSchmidtOrthogonalization(dir_obj, dir_eqs, dir_ineqs)

        # Determine position of border of halfspaces and hyperplanes
        pos_obj, pos_eqs, pos_ineqs = self.__DetermineConstraintBorders(len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs)

        # Translate to orthogonal basis
        self.pos_obj_o = TranslateToNewBasis(pos_obj, self.ortho_basis)
        self.dir_obj_o = TranslateToNewBasis([dir_obj], self.ortho_basis)
        self.pos_eqs_o = TranslateToNewBasis(pos_eqs, self.ortho_basis)
        self.dir_eqs_o = TranslateToNewBasis(dir_eqs, self.ortho_basis)
        self.pos_ineqs_o = TranslateToNewBasis(pos_ineqs, self.ortho_basis)
        self.dir_ineqs_o = TranslateToNewBasis(dir_ineqs, self.ortho_basis)

    # --------------------------------------------------------------------------
    def RunProjection(self, len_obj, threshold, nargout):

        # Adjust halfspaces according input
        pos_obj_o, pos_eqs_o, pos_ineqs_o, len_obj, len_eqs, len_ineqs = self.__AdjustHalfSpaces(len_obj, threshold)

        pos_HS_o = HorzCat(pos_obj_o, pos_ineqs_o)
        dir_HS_o = HorzCat(dir_obj_o, dir_ineqs_o)

        # Project to adjusted halfspaces (suboptimization)
        dX_o, subopt_itr, error, exit_code = self.__ProjectToHalfSpaces(pos_HS_o, dir_HS_o)

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
        else:
            return norm_dX, dX, is_projection_sucessfull

    # --------------------------------------------------------------------------
    def __AdjustHalfSpaces(self, len_obj, threshold):

        # Adjusting objective halfspace according input length
        adj_len_obj = len_obj
        adj_pos_obj_o = ScalarVectorProduct(adj_len_obj/self.reference_len_obj, self.pos_obj_o)

        # Adjusting all half spaces according threshold
        if threshold<adj_len_obj:
            adj_len_obj = threshold
            adj_pos_obj_o = ScalarVectorProduct(adj_len_obj/self.reference_len_obj, self.pos_obj_o)

        adj_len_eqs = self.reference_len_eqs
        adj_len_ineqs = self.reference_len_ineqs

        for itr in range(self.num_eqs):
            adj_len_eqs[itr] = min(max(self.reference_len_eqs[itr],-threshold),threshold)
            scaling = adj_len_eqsr[itr]/self.reference_len_eqs[itr]
            adj_pos_eqs_o[itr] =  ScalarVectorProduct(scaling, self.pos_eqs_o[itr])

        for itr in range(self.num_ineqs):
            adj_len_ineqs[itr] = min(self.reference_len_ineqs[itr],threshold)
            scaling = adj_len_ineqs[itr]/self.reference_len_ineqs[itr]
            adj_pos_ineqs_o[itr] =  ScalarVectorProduct(scaling, self.pos_ineqs_o[itr])

        return adj_pos_obj_o, adj_pos_eqs_o, adj_pos_ineqs_o, adj_len_obj, adj_len_eqs, adj_len_ineqs

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
            norm2_v = Norm2(V[0])
            if norm2_v>1e-10:
                B.append( ScalarVectorProduct(1/norm2_v,v) )
            else:
                print("Zero basis vector after Gram-Schmidt orthogonalization --> Problem set has too many response functions or involves only parallel vectors or is at the unconstrained optimum.")
                B.append(v)

        return B

    # --------------------------------------------------------------------------
    def __DetermineConstraintBorders(self, len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs):
        pos_obj = []
        pos_eq = []
        pos_ineq = []

        pos_obj.append(ScalarVectorProduct(-len_obj,dir_obj))

        for i in range(self.num_eqs):
            pos_eq.append(ScalarVectorProduct(-len_eqs[i],dir_eqs[i]))

        for i in range(self.num_ineqs):
            pos_ineq.append(ScalarVectorProduct(-len_ineqs[i],dir_ineqs[i]))

        return pos_obj, pos_eq, pos_ineq

    # --------------------------------------------------------------------------
    def __ProjectToHalfSpaces(self, hsPos_o, hs_Dir_o):
        A = Trans(hs_Dir_o)
        b = [ Dot(hsPos_o[i],hs_Dir_o[i]) for i in range(RowSize(A)) ]

        max_itr = self.algorithm_settings["subopt_max_itr"].GetInt()
        tolerance = self.algorithm_settings["subopt_tolerance"].GetDouble()

        dX_o, subopt_itr, error, exit_code = QuadProg(A, b, max_itr, tolerance)

        return dX_o, subopt_itr, error, exit_code

# ==============================================================================




    # # --------------------------------------------------------------------------
    # def __CheckConvergence(self, val_eqs, val_ineqs):

    #     step_length_tolerance = self.algorithm_settings["step_length_tolerance"].GetDouble()
    #     feasibility_frequency = self.algorithm_settings["feasibility_frequency"].GetDouble()
    #     max_iterations = self.algorithm_settings["max_iterations"].GetInt()

    #     # Check termination conditions
    #     if self.opt_iteration == max_iterations:
    #         self.algo_is_in_termination_phase = True
    #         print("> Algorithm: Begin enforcing feasibility for termination (max number of optimization iterations reached)!")

    #     if self.opt_iteration > 1 and self.value_history["step_length"][-1] < step_length_tolerance:
    #         self.algo_is_in_termination_phase = True
    #         print("> Algorithm: Begin enforcing feasibility for termination (step length < tolerance)!")

    #     # Check feasibility
    #     if self.algo_is_in_termination_phase:
    #         self.algo_is_in_feasibility_phase = True

    #     if self.opt_iteration == self.next_iteration_to_enforce_feasibility:
    #         self.algo_is_in_feasibility_phase = True
    #         self.reduceOscillationFactor = 1
    #         print("> Algorithm: Begin enforcing feasibility according specified frequency!")

    #     feasibility_is_enforced = self.__IsFeasibilityEnforced(val_eqs, val_ineqs)

    #     if feasibility_is_enforced and self.algo_is_in_feasibility_phase:
    #         print("> Algorithm: Feasibility was enforced")
    #         self.algo_is_in_feasibility_phase = False
    #         self.next_iteration_to_enforce_feasibility = iterCounter + feasibility_frequency

    #     if feasibility_is_enforced and not self.algo_is_in_feasibility_phase:
    #         print("> Algorithm: Feasibility is spontaneously enforced.")

    #     if self.algo_is_in_termination_phase and feasibility_is_enforced:
    #         print("> Algorithm: Converged!")
    #         return True
    #     else:
    #         False

    # # --------------------------------------------------------------------------
    # def __IsFeasibilityEnforced(self, val_eqs, val_ineqs):
    #     # Return false if any equality constraint is violated
    #     for itr, eq in enumerate(self.equality_constraints):
    #         eq_value = val_eqs[itr]
    #         eq_id = eq["settings"]["identifier"].GetString()
    #         eq_tolerance = eq["settings"]["tolerance"].GetDouble()
    #         eq_reference = self.communicator.getReferenceValue(eq_id)
    #         if abs(eq_value) > eq_tolerance*eq_reference:
    #             return False

    #     # Return false if any inequality constraint is violated
    #     for val_ineq in val_ineqs:
    #         if val_ineq > 0:
    #             return False

    #     return True

    # # --------------------------------------------------------------------------
    # def __ApplyStepLengthRule(self, val_obj, len_eqs, len_ineqs):

    #     max_step_length = self.algorithm_settings["max_step_length"].GetDouble()
    #     step_length_decrease_factor = self.algorithm_settings["step_length_decrease_factor"].GetDouble()
    #     step_length_increase_factor = self.algorithm_settings["step_length_increase_factor"].GetDouble()

    #     step_length = max_step_length

    #     # dont change step lenght if final feasibility is enforced
    #     if self.algo_is_in_feasibility_phase:
    #         self.last_itr_with_step_length_update = self.opt_iteration
    #         return step_length

    #     # Wait with an updat of the step length for min three iterations
    #     if self.opt_iteration<3 or self.opt_iteration<self.last_itr_with_step_length_update+2:
    #         return step_length

    #     # Decrease step length if necessary conditions are met
    #     if self.__IsStepLengthToBeDecreased(val_obj, len_eqs, len_ineqs):
    #         self.last_itr_with_step_length_update = self.opt_iteration
    #         print("> Algorithm: Decreasing step length!")
    #         return step_length / step_length_decrease_factor

    #     # Increase slowly step length if necessary conditions are met
    #     if self.__IsStepLengthToBeIncreased(val_obj, len_eqs, len_ineqs):
    #         self.last_itr_with_step_length_update = self.opt_iteration
    #         print("> Algorithm: Increasing step length")
    #         return min(max_step_length, step_length * step_length_increase_factor)

    #     #do nothing
    #     return step_length

    # # --------------------------------------------------------------------------
    # def __IsStepLengthToBeDecreased(self, val_obj, len_eqs, len_ineqs):

    #     # Check if enough values are already computed
    #     if self.opt_iteration < 3:
    #         return False

    #     v2, v1, v0 = self.value_history["val_obj"][-2:] + [val_obj]

    #     # Check if objective oscillates
    #     if (v1 < v2 and v1 < v0) and v0 > (v2-0.1*(v2-v1)):
    #         objective_oscillates = True
    #     else:
    #         objective_oscillates = False

    #     # Objective does not improve
    #     if v0>v2:
    #         no_progress_in_obj = True
    #     else:
    #         no_progress_in_obj = False

    #     # No progress in the distance to the feasible domain
    #     s_old = self.value_history["step_length"][-2]

    #     old_len_eqs = self.value_history["len_eqs"][-2]
    #     old_len_ineqs = self.value_history["len_ineqs"][-2]

    #     dist_feasible_domoain = self.__SimpleDistanceToFeasibleDomain(len_eqs,len_ineqs)
    #     dist_feasible_domoain_old = self.__SimpleDistanceToFeasibleDomain(old_len_eqs , old_len_ineqs)

    #     if dist_feasible_domoain > (dist_feasible_domoain_old - 0.2*s_old):
    #         no_progress_towards_feasible_domain = True
    #     else:
    #         no_progress_towards_feasible_domain = False

    #     # Evaluate individual conditions
    #     if (no_progress_in_obj or objective_oscillates) and no_progress_towards_feasible_domain:
    #         return True
    #     else:
    #         return False

    # # --------------------------------------------------------------------------
    # def __SimpleDistanceToFeasibleDomain(self, len_equality, len_inequality):
    #     return max( [max(l,0) for l in len_inequality] + [abs(l) for l in len_equality] + [0] ) # add [0] for case with no constraint

    # # --------------------------------------------------------------------------
    # def __IsStepLengthToBeIncreased(self, val_obj, len_eqs, len_ineqs):
    #     # Progress of objective function is monotonic over a few iterations
    #     number_of_smooth_steps = 5
    #     if self.opt_iteration < number_of_smooth_steps:
    #         return False
    #     else:
    #         val_hist = self.value_history["val_obj"][-number_of_smooth_steps+1:] + [val_obj]

    #         # Check if values are changing monotonically
    #         is_oscillating = []
    #         for i in range(number_of_smooth_steps-2):
    #             is_oscillating.append( (val_hist[i+1]-val_hist[i])*(val_hist[i+2]-val_hist[i+1])<0 )

    #         if sum(is_oscillating)==0: # no oscillation at all
    #             return True
    #         else:
    #             return False
