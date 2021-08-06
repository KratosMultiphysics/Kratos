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
from .algorithm_base import OptimizationAlgorithm
from . import mapper_factory
from . import data_logger_factory
from .custom_timer import Timer
from .custom_variable_utilities import WriteDictionaryDataOnNodalVariable

import math
#import matplotlib.pyplot as plt

# ==============================================================================
class AlgorithmSequentialQuadraticProgrammingWithLineSearch(OptimizationAlgorithm):
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
        self.constraint_laplace_multipliers = {}
        self.constraint_penalty_multipliers = {}
        for itr, constraint in enumerate(self.constraints):
            self.constraint_gradient_variables.update({
                constraint["identifier"].GetString() : {
                    "gradient": KM.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX"),
                    "mapped_gradient": KM.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")
                }
            })
            self.constraint_laplace_multipliers.update({
                constraint["identifier"].GetString() : {
                    "0": 0.0,
                    "-1": 0.0
                }
            })
            self.constraint_penalty_multipliers.update({
                constraint["identifier"].GetString() : 0.0
            })
        self.quadratic_approximation = []
        self.initial_curvature = 1.0
        self.previous_laplace_gradient = KM.Matrix()
        self.last_iteration_count = 0
        self.n_penalty_evaluations = 0
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
            
            self.__computePenaltyFunction()
            
            self.__logCurrentOptimizationStep()
            
            self.__updateQuadraticApproximation()

            self.__computeShapeUpdate()
            
            self.__performLineSearch()

            KM.Logger.Print("")
            KM.Logger.PrintInfo("ShapeOpt", "Time needed for current optimization step = ", timer.GetLapTime(), "s")
            KM.Logger.PrintInfo("ShapeOpt", "Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

            if self.__isAlgorithmConverged():
                break
            #else:
            #    self.__determineAbsoluteChanges()

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
    def __analyzeShapeWithoutGradients(self):
        self.communicator.initializeCommunication()
        self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())

        for constraint in self.constraints:
            con_id =  constraint["identifier"].GetString()
            self.communicator.requestValueOf(con_id)

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.optimization_iteration, self.communicator)

    # --------------------------------------------------------------------------
    def __updateQuadraticApproximation(self):
        if self.optimization_iteration == 1:
            objective_id = self.objectives[0]["identifier"].GetString()
            objective_value = self.communicator.getStandardizedValue(objective_id)
            
            self.initial_curvature = objective_value / (self.step_size)**2
            
        if self.optimization_iteration == 2:
            
            laplace_gradient = KM.Matrix()
            self.optimization_utilities.AssembleVectorMatrix(self.design_surface, laplace_gradient, KSO.DF1DX_MAPPED)
            
            g_gradient = KM.Matrix()
            
            for constraint in self.constraints:
                identifier = constraint["identifier"].GetString()
                
                self.optimization_utilities.AssembleVectorMatrix(self.design_surface, g_gradient, self.constraint_gradient_variables[identifier]["mapped_gradient"])
                
                lap = self.constraint_laplace_multipliers[identifier]["0"]
                self.constraint_penalty_multipliers[identifier] = max(abs(lap),0.5*(abs(lap)+self.constraint_penalty_multipliers[identifier]))
                
                #if self.communicator.getStandardizedValue(identifier) > 0 :
                #   laplace_gradient += g_gradient * self.constraint_penalty_multipliers[identifier]
                
                laplace_gradient += g_gradient*self.constraint_laplace_multipliers[identifier]["0"]
                
            self.previous_laplace_gradient = laplace_gradient
            
        elif self.optimization_iteration > 2:
            
            p = KM.Matrix()
            self.optimization_utilities.AssembleVectorMatrix(self.design_surface, p, KSO.SEARCH_DIRECTION)
            p_vec = KM.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, p_vec, KSO.SEARCH_DIRECTION)
            
            laplace_gradient = KM.Matrix()
            self.optimization_utilities.AssembleVectorMatrix(self.design_surface, laplace_gradient, KSO.DF1DX_MAPPED)
            
            g_gradient = KM.Matrix()
            
            for constraint in self.constraints:
                identifier = constraint["identifier"].GetString()
                
                self.optimization_utilities.AssembleVectorMatrix(self.design_surface, g_gradient, self.constraint_gradient_variables[identifier]["mapped_gradient"])
                
                lap = self.constraint_laplace_multipliers[identifier]["0"]
                self.constraint_penalty_multipliers[identifier] = max(abs(lap),0.5*(abs(lap)+self.constraint_penalty_multipliers[identifier]))
                
                #if self.communicator.getStandardizedValue(identifier) > 0 :
                #    laplace_gradient += g_gradient * self.constraint_penalty_multipliers[identifier]
                
                laplace_gradient += g_gradient*self.constraint_laplace_multipliers[identifier]["0"]
            
            y = laplace_gradient - self.previous_laplace_gradient
            
            self.previous_laplace_gradient = laplace_gradient
            
            
            pBp_mat = KM.Matrix()
            
            Bp = KM.Matrix()
            Bp = self.__multiplyQuadraticApproximation(p)
            
            self.optimization_utilities.MatrixScalarProduct(pBp_mat, p, Bp)
            pBp = pBp_mat[0,0]
            
            py_mat = KM.Matrix()
            self.optimization_utilities.MatrixScalarProduct(py_mat, p, y)
            py = py_mat[0,0]
            
            if py < 0.2*pBp:
                omega = 0.8*pBp/(pBp-py)
            else:
                omega = 1.0;
            eta = y*omega + Bp*(1-omega)
            
            peta_mat = KM.Matrix()
            self.optimization_utilities.MatrixScalarProduct(peta_mat, p, eta)
            peta = peta_mat[0,0]
            
            eta *= 1.0/math.sqrt(peta)
            self.quadratic_approximation.append(eta)
                      
            Bp *= 1.0/math.sqrt(pBp)
            self.quadratic_approximation.append(Bp)
            
    # --------------------------------------------------------------------------
    def __multiplyQuadraticApproximation(self, rightVector):
        
        result = KM.Matrix()
        result.Resize(rightVector.Size1(),1)
        result.fill(0.0) 
        
        result += rightVector * self.initial_curvature
        
        for i in range(len(self.quadratic_approximation)):
            currentVector = self.quadratic_approximation[i]
            
            factor_mat = KM.Matrix()
            self.optimization_utilities.MatrixScalarProduct(factor_mat, currentVector, rightVector)
            factor = factor_mat[0,0]
            
            if i%2 == 0:
                # add eta terms
                result += currentVector * factor
            else:
                # subtract Bp terms
                result -= currentVector * factor
            
        return result
        
        
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
        self.model_part_controller.DampNodalVariableIfSpecified(KSO.SHAPE_UPDATE)

    # --------------------------------------------------------------------------
    def __computeControlPointUpdate(self):
        g, nabla_g = self.__getConstraints()
        
        objective_id = self.objectives[0]["identifier"].GetString()
        f = self.communicator.getStandardizedValue(objective_id)

        KM.Logger.PrintInfo("ShapeOpt", "Assemble vector of objective gradient.")
        nabla_f = KM.Matrix()
        self.optimization_utilities.AssembleVectorMatrix(self.design_surface, nabla_f, KSO.DF1DX_MAPPED)



        settings = KM.Parameters('{ "solver_type" : "LinearSolversApplication.dense_col_piv_householder_qr" }')
        solver = dense_linear_solver_factory.ConstructSolver(settings)

        KM.Logger.PrintInfo("ShapeOpt", "Calculate projected search direction and correction.")
        
        # ---------------------------------------------------------------------------------
        
        S = KM.Matrix()
        S.Resize(nabla_f.Size1(),1)
        S.fill(0.0)
        
        S_int = KM.Matrix()
        S_int.Resize(nabla_f.Size1(),1)
        S_int.fill(0.0)
        
        alpha_opt = 0.0
        
        #used for debugging
        constraint_history = []
        objective_history = []
        for i in range(len(nabla_g)):
            constraint_history.append([])
        
        N = KM.Matrix()
        activeConstraintMap = []
        
        self.last_iteration_count = 10000
        for it in range(10000):
            print("Iteration: ", it)
            activeConstraintValues = []
            activeConstraintGradient = []
            
            S = S +  S_int * alpha_opt
            
            g_int = []
            for i in range(len(nabla_g)):
                current_nabla_g = nabla_g[i]
                
                current_nabla_g_vec = KM.Vector()
                self.optimization_utilities.AssembleVectorMatrixtoVector(current_nabla_g_vec,current_nabla_g)
                g_norm = current_nabla_g_vec.norm_inf()
                
                g_int_mat = KM.Matrix()
                self.optimization_utilities.MatrixScalarProduct(g_int_mat, current_nabla_g, S)
                g_int.append(g_int_mat[0,0])
                
                if g[i] < 0:
                    g_int[i] += g[i];
                else:
                    g_int[i] += 0.95*g[i];
                
                constraint_history[i].append(g_int[i])

                if g_int[i] > -1e-10 and g_norm > 1e-10:
                    activeConstraintMap.append(i)
                    activeConstraintValues.append(g_int[i])
                    activeConstraintGradient.append(current_nabla_g)
            

            S_int = KM.Matrix()
            C_int = KM.Vector()
            
            BS = KM.Matrix()
            BS = self.__multiplyQuadraticApproximation(S)

            f_int_mat = KM.Matrix()
            self.optimization_utilities.MatrixScalarProduct(f_int_mat, nabla_f, S)
            f_int = f + f_int_mat[0,0]
            SBS_mat = KM.Matrix()
            self.optimization_utilities.MatrixScalarProduct(SBS_mat, S, BS)
            f_int += 0.5 * SBS_mat[0,0]
            objective_history.append(f_int)
            
            
            if len(activeConstraintMap) > 0:
                nabla_f_int = KM.Matrix()
                nabla_f_int = nabla_f + BS

                self.optimization_utilities.AssembleVectorMatrixtoMatrix(self.design_surface, N, activeConstraintGradient)
                
                self.optimization_utilities.CalculateProjectedSearchDirectionAndCorrectionMatrix(
				nabla_f_int,
				N,
				activeConstraintValues,
				solver,
				S_int,
				C_int
                )
                
                c_norm = C_int.norm_inf()
                if c_norm > 1e-10:
                    C_mat = KM.Matrix()
                    self.optimization_utilities.AssembleVectortoVectorMatrix(C_mat,C_int)
                    S_int = C_mat
            else:
                S_int = nabla_f*(-1.0) - BS
            
            S_int_vec = KM.Vector()
            self.optimization_utilities.AssembleVectorMatrixtoVector(S_int_vec,S_int)
            s_norm = S_int_vec.norm_inf()
            
            dSBdS_mat = KM.Matrix()
            
            BdS = KM.Matrix()
            BdS = self.__multiplyQuadraticApproximation(S_int)
            
            self.optimization_utilities.MatrixScalarProduct(dSBdS_mat, S_int, BdS)
            dSBdS = dSBdS_mat[0,0]
	    	
            FdS_mat = KM.Matrix()
            self.optimization_utilities.MatrixScalarProduct(FdS_mat, S_int, nabla_f)
            FdS = FdS_mat[0,0]
	    	
            dSBS_mat = KM.Matrix()
            self.optimization_utilities.MatrixScalarProduct(dSBS_mat, S_int, BS)
            dSBS = dSBS_mat[0,0]
            #print("objective gradient approximation along seach direction: ",FdS + dSBS)
     
            alpha_opt = 0.0
            
            if dSBdS < 0.0:
                self.last_iteration_count = it
                print("error negative curvature")
                break
	    	
            alpha_opt = (-1)*(FdS + dSBS)/dSBdS
            if alpha_opt < 0:
                alpha_opt = 0.0
	    	
            allConstraintsFullfilled = True
            for i in range(len(nabla_g)):
                current_nabla_g = nabla_g[i]
                
                dgdS_mat = KM.Matrix()
                self.optimization_utilities.MatrixScalarProduct(dgdS_mat, current_nabla_g, S_int)
                dgdS = dgdS_mat[0,0]
                            
                if g_int[i] > 1e-10:
                    allConstraintsFullfilled = False
 	    		
                if abs(dgdS) > 1e-10:

                    current_alpha = (-1)*(g_int[i])/dgdS

                    if current_alpha > 0 and (current_alpha < alpha_opt or alpha_opt == 0.0):
                        alpha_opt = current_alpha

            if allConstraintsFullfilled and FdS + dSBS > -1e-12:
                self.last_iteration_count = it
                break
            if s_norm < 1e-10:
                self.last_iteration_count = it
                break
            if alpha_opt * s_norm < 1e-10:
                self.last_iteration_count = it
                break
        
        #plt.figure()
        #plt.plot(objective_history)
        #plt.show()
        #for i in range(len(nabla_g)):
            #if len(constraint_history[i]) > 0:
                #plt.figure()
                #plt.plot(constraint_history[i])
                #plt.show()
        
        lambda_mat = KM.Matrix()
        
        if len(activeConstraintMap) > 0:
            self.optimization_utilities.CalculateLaplaceMultipliers(lambda_mat, N, nabla_f+BS, solver)
            
            j = 0
            for i, constraint in enumerate(self.constraints):
                identifier = constraint["identifier"].GetString()
                
                if i == activeConstraintMap[j]:
                    self.constraint_laplace_multipliers[identifier]["-1"] = self.constraint_laplace_multipliers[identifier]["0"]
                    self.constraint_laplace_multipliers[identifier]["0"] = lambda_mat[j,0]
                    
                    j += 1
                    if j >= len(activeConstraintMap):
                        break
                else:
                    self.constraint_laplace_multipliers[identifier]["-1"] = self.constraint_laplace_multipliers[identifier]["0"]
                    self.constraint_laplace_multipliers[identifier]["0"] = 0.0

        else:
            for constraint in self.constraints:
                identifier = constraint["identifier"].GetString()
                
                self.constraint_laplace_multipliers[identifier]["-1"] = self.constraint_laplace_multipliers[identifier]["0"]
                self.constraint_laplace_multipliers[identifier]["0"] = 0.0
        

        self.optimization_utilities.AssignMatrixToVariable(self.design_surface, S, KSO.SEARCH_DIRECTION)
        S_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.SEARCH_DIRECTION)
        if S_norm > self.step_size:
            S *= self.step_size/S_norm
            self.optimization_utilities.AssignMatrixToVariable(self.design_surface, S, KSO.SEARCH_DIRECTION)
        self.optimization_utilities.AssignMatrixToVariable(self.design_surface, S, KSO.CORRECTION)
        self.optimization_utilities.AssignMatrixToVariable(self.design_surface, S, KSO.CONTROL_POINT_UPDATE)

    # --------------------------------------------------------------------------
    def __getConstraints(self):
        active_constraint_values = []
        active_constraint_variables = []

        for constraint in self.constraints:
            g_a_variable_matrix = KM.Matrix()
            identifier = constraint["identifier"].GetString()
            constraint_value = self.communicator.getStandardizedValue(identifier)
            active_constraint_values.append(constraint_value)
            g_a_variable = self.constraint_gradient_variables[identifier]["mapped_gradient"]
            self.optimization_utilities.AssembleVectorMatrix(self.design_surface, g_a_variable_matrix, g_a_variable)
            active_constraint_variables.append(
                g_a_variable_matrix)

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
    def __performLineSearch(self):
        S = KM.Matrix()
        self.optimization_utilities.AssembleVectorMatrix(self.design_surface, S, KSO.SEARCH_DIRECTION)
        
        S_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.SEARCH_DIRECTION)
        
        initialValue = self.__computePenaltyFunction()
        
        self.n_penalty_evaluations = 1
        alpha_upper = 1.0
        
        self.__initializeNewShape()
        self.__analyzeShapeWithoutGradients()
        
        upperValue = self.__computePenaltyFunction()
        noUpperBound = False
        
        if upperValue < initialValue:
            oneStep = False
            
            middleValue = upperValue + 1.0
            while middleValue > upperValue:
                
                alpha_mid = alpha_upper
                alpha_upper = 1.61803 * alpha_upper
                

                self.n_penalty_evaluations += 1
                
                self.optimization_utilities.AssignMatrixToVariable(self.design_surface, S*(alpha_upper - alpha_mid), KSO.CONTROL_POINT_UPDATE)
                
                # ensure that the incremental update is equalto the complete update
                C = KM.Matrix()
                self.optimization_utilities.AssembleVectorMatrix(self.design_surface, C, KSO.CORRECTION)
                C = C + S*(alpha_upper - alpha_mid)
                self.optimization_utilities.AssignMatrixToVariable(self.design_surface, C, KSO.CORRECTION)
                
                self.mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
                self.model_part_controller.DampNodalVariableIfSpecified(KSO.SHAPE_UPDATE)
                self.__determineAbsoluteChanges()
                self.__initializeNewShape()
                self.__analyzeShapeWithoutGradients()
                
                middleValue = upperValue
                upperValue = self.__computePenaltyFunction()
                
                if alpha_upper*S_norm >= self.step_size:
                    noUpperBound = True
                    break
        else:

            oneStep = True
            self.n_penalty_evaluations += 1
            
            alpha_mid = 1.0/1.61803
            
            self.optimization_utilities.AssignMatrixToVariable(self.design_surface, S*(alpha_mid - alpha_upper), KSO.CONTROL_POINT_UPDATE)
            
            # ensure that the incremental update is equalto the complete update
            C = KM.Matrix()
            self.optimization_utilities.AssembleVectorMatrix(self.design_surface, C, KSO.CORRECTION)
            C = C + S*(alpha_mid - alpha_upper)
            self.optimization_utilities.AssignMatrixToVariable(self.design_surface, C, KSO.CORRECTION)
            
            self.mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
            self.model_part_controller.DampNodalVariableIfSpecified(KSO.SHAPE_UPDATE)
            self.__determineAbsoluteChanges()
            self.__initializeNewShape()
            self.__analyzeShapeWithoutGradients()
            
            middleValue = self.__computePenaltyFunction()
            
        if not noUpperBound:   
            
            a2 = (((upperValue - initialValue)/alpha_upper)-((middleValue-initialValue)/alpha_mid))/(alpha_upper-alpha_mid)
            a1 = (middleValue - initialValue)/alpha_mid - a2 * alpha_mid
            
            alpha_opt = - a1/(2*a2)
            
            if oneStep:
                self.optimization_utilities.AssignMatrixToVariable(self.design_surface, S*(alpha_opt - alpha_mid), KSO.CONTROL_POINT_UPDATE)
                self.optimization_utilities.AssignMatrixToVariable(self.design_surface, S*alpha_opt, KSO.SEARCH_DIRECTION)
                
                C = KM.Matrix()
                self.optimization_utilities.AssembleVectorMatrix(self.design_surface, C, KSO.CORRECTION)
                C = C + S*(alpha_opt - alpha_mid)
                self.optimization_utilities.AssignMatrixToVariable(self.design_surface, C, KSO.CORRECTION)
                
                self.mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
                self.model_part_controller.DampNodalVariableIfSpecified(KSO.SHAPE_UPDATE)
                self.__determineAbsoluteChanges()
            else:
                self.optimization_utilities.AssignMatrixToVariable(self.design_surface, S*(alpha_opt - alpha_upper), KSO.CONTROL_POINT_UPDATE)
                self.optimization_utilities.AssignMatrixToVariable(self.design_surface, S*alpha_opt, KSO.SEARCH_DIRECTION)
                
                C = KM.Matrix()
                self.optimization_utilities.AssembleVectorMatrix(self.design_surface, C, KSO.CORRECTION)
                C = C + S*(alpha_opt - alpha_upper)
                self.optimization_utilities.AssignMatrixToVariable(self.design_surface, C, KSO.CORRECTION)
                
                self.mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
                self.model_part_controller.DampNodalVariableIfSpecified(KSO.SHAPE_UPDATE)
                self.__determineAbsoluteChanges()
        else:
            self.optimization_utilities.AssignMatrixToVariable(self.design_surface, S*0.0, KSO.CONTROL_POINT_UPDATE)
            self.optimization_utilities.AssignMatrixToVariable(self.design_surface, S*alpha_upper, KSO.SEARCH_DIRECTION)
            
            C = KM.Matrix()
            self.optimization_utilities.AssembleVectorMatrix(self.design_surface, C, KSO.CORRECTION)
            C = C + S*0.0
            self.optimization_utilities.AssignMatrixToVariable(self.design_surface, C, KSO.CORRECTION)
            
            self.mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
            self.model_part_controller.DampNodalVariableIfSpecified(KSO.SHAPE_UPDATE)
            self.__determineAbsoluteChanges()
            
        
    # --------------------------------------------------------------------------
    def __computePenaltyFunction(self):
        objective_id = self.objectives[0]["identifier"].GetString()
        penaltyValue = self.communicator.getStandardizedValue(objective_id)
        #print("Starting penalty value computation:")
        #print("objective value: ",penaltyValue)
        for constraint in self.constraints:
                identifier = constraint["identifier"].GetString()
                constraint_value = self.communicator.getStandardizedValue(identifier)
                #print("constraint value: ", constraint_value)
                if constraint_value > 0 :
                    penaltyValue += constraint_value * self.constraint_penalty_multipliers[identifier]
                    #print("penalty value: ", penaltyValue)
        #print("Finished penalty value computation:")
        return penaltyValue
    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.step_size
        additional_values_to_log["inf_norm_s"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.SEARCH_DIRECTION)
        additional_values_to_log["inf_norm_c"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.CORRECTION)
        additional_values_to_log["iteration_count"] = self.last_iteration_count
        additional_values_to_log["line_search_function_calls"] = self.n_penalty_evaluations
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
            print("rel objective change: ",relative_change_of_objective_value)
            print("tolerance: ",self.relative_tolerance)
            if abs(relative_change_of_objective_value) < self.relative_tolerance:
                KM.Logger.Print("")
                KM.Logger.PrintInfo("ShapeOpt", "Optimization problem converged within a relative objective tolerance of ",self.relative_tolerance,"%.")
                return True

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges(self):
        self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.CONTROL_POINT_UPDATE, KSO.CONTROL_POINT_CHANGE)
        self.optimization_utilities.AddFirstVariableToSecondVariable(self.design_surface, KSO.SHAPE_UPDATE, KSO.SHAPE_CHANGE)

# ==============================================================================
