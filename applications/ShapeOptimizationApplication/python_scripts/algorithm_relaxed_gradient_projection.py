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


# ==============================================================================
class AlgorithmRelaxedGradientProjection(OptimizationAlgorithm):
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
        self.constraint_buffer_variables = {}
        self.constraint_laplace_multipliers = {}
        for itr, constraint in enumerate(self.constraints):
            constraint_id = constraint["identifier"].GetString()
            self.constraint_gradient_variables.update({
                constraint_id : {
                    "gradient": KM.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX"),
                    "mapped_gradient": KM.KratosGlobals.GetVariable("DC"+str(itr+1)+"DX_MAPPED")
                }
            })
            self.constraint_buffer_variables.update({
                constraint_id : {
                    "buffer_value": 0.0,
                    "buffer_value-1": 0.0,
                    "buffer_size": 1e-12,
                    "buffer_size_factor": 2.0,
                    "central_buffer_value": 0.0,
                    "lower_buffer_value": - 1e-12,
                    "upper_buffer_value": 1e-12,
                    "constraint_value-1": 0.0,
                    "constraint_value-2": 0.0,
                    "constraint_value-3": 0.0,
                    "max_constraint_change": 0.0
                }
            })
            self.constraint_laplace_multipliers.update({
                constraint_id : 0.0
            })

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
            
            #not required
            self.__computeLaplaceMultipliers()
            
            self.__computeBufferValue()

            self.__computeShapeUpdate()
            
            self.__normOutputGradients()
            
            self.__logCurrentOptimizationStep()
            
            self.__updateBufferZone()
            


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
    def __computeLaplaceMultipliers(self):
        
        active_constraint_gradient = []
        active_constraint_map = []
        
        for i, constraint in enumerate(self.constraints):
            identifier = constraint["identifier"].GetString()
            constraint_value = self.communicator.getStandardizedValue(identifier)
            
            if constraint["type"].GetString() == "=":
                active_constraint_gradient.append(self.constraint_gradient_variables[identifier]["mapped_gradient"])
                active_constraint_map.append(i)
            elif constraint_value >= -1e-12:
                active_constraint_gradient.append(self.constraint_gradient_variables[identifier]["mapped_gradient"])
                active_constraint_map.append(i)
            
        if len(active_constraint_map) > 0:
            N = KM.Matrix()
            self.optimization_utilities.AssembleMatrix(self.design_surface, N, active_constraint_gradient)
            
            nabla_f = KM.Matrix()
            self.optimization_utilities.AssembleVectorMatrix(self.design_surface, nabla_f, KSO.DF1DX_MAPPED)
            
            settings = KM.Parameters('{ "solver_type" : "LinearSolversApplication.dense_col_piv_householder_qr" }')
            solver = dense_linear_solver_factory.ConstructSolver(settings)
            
            lambda_mat = KM.Matrix()
            
            self.optimization_utilities.CalculateLaplaceMultipliers(lambda_mat, N, nabla_f, solver)
            
            j = 0
            for i, constraint in enumerate(self.constraints):
                identifier = constraint["identifier"].GetString()
                if i == active_constraint_map[j]:
                    self.constraint_laplace_multipliers[identifier] = lambda_mat[j,0]
                    j += 1
                    if j >= len(active_constraint_map):
                        break
                else:
                    self.constraint_laplace_multipliers[identifier] = 0.0
        else:
            for constraint in self.constraints:
                identifier = constraint["identifier"].GetString()
                self.constraint_laplace_multipliers[identifier] = 0.0
    
    # --------------------------------------------------------------------------
    def __computeBufferValue(self):
    	for constraint in self.constraints:
            identifier = constraint["identifier"].GetString()
            constraint_value = self.communicator.getStandardizedValue(identifier)
            
            self.constraint_buffer_variables[identifier]["buffer_value-1"] = self.constraint_buffer_variables[identifier]["buffer_value"]
            
            if self.optimization_iteration > 1:
                
                if abs(constraint_value - self.constraint_buffer_variables[identifier]["constraint_value-1"]) > self.constraint_buffer_variables[identifier]["max_constraint_change"]:
                    self.constraint_buffer_variables[identifier]["max_constraint_change"] = abs(constraint_value - self.constraint_buffer_variables[identifier]["constraint_value-1"])
                    
                self.constraint_buffer_variables[identifier]["buffer_size"] = max(self.constraint_buffer_variables[identifier]["buffer_size_factor"] * self.constraint_buffer_variables[identifier]["max_constraint_change"], 1e-12)
            
            self.constraint_buffer_variables[identifier]["lower_buffer_value"] = self.constraint_buffer_variables[identifier]["central_buffer_value"] - self.constraint_buffer_variables[identifier]["buffer_size"]
            self.constraint_buffer_variables[identifier]["upper_buffer_value"] = self.constraint_buffer_variables[identifier]["central_buffer_value"] + self.constraint_buffer_variables[identifier]["buffer_size"]
            
            if self.__isConstraintActive(constraint):                
                if constraint["type"].GetString() == "=":
                    self.constraint_buffer_variables[identifier]["buffer_value"] = 1 - abs(constraint_value) / self.constraint_buffer_variables[identifier]["buffer_size"]
                else:
                    self.constraint_buffer_variables[identifier]["buffer_value"] = (constraint_value - self.constraint_buffer_variables[identifier]["lower_buffer_value"]) / self.constraint_buffer_variables[identifier]["buffer_size"]
            else:
                self.constraint_buffer_variables[identifier]["buffer_value"] = 0.0
            
            #print("DEBUG =========================================================")
            #print("constraint: ", identifier)
            #print("buffer value: ", self.constraint_buffer_variables[identifier]["buffer_value"])
            #print("buffer size: ", self.constraint_buffer_variables[identifier]["buffer_size"])
            #print("buffer size factor: ", self.constraint_buffer_variables[identifier]["buffer_size_factor"])
            #print("central buffer value: ", self.constraint_buffer_variables[identifier]["central_buffer_value"])
            #print("lower buffer value: ", self.constraint_buffer_variables[identifier]["lower_buffer_value"])
            #print("END ===========================================================")
                
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
        """adapted from https://msulaiman.org/onewebmedia/GradProj_2.pdf"""
        g_a, g_a_variables, relaxation_coefficients, correction_coefficients = self.__getActiveConstraints()
        
        #print("DEBUG: =========================================== ", relaxation_coefficients)
        #print("DEBUG: =========================================== ", correction_coefficients)
        KM.Logger.PrintInfo("ShapeOpt", "Assemble vector of objective gradient.")
        nabla_f = KM.Vector()
        s = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, nabla_f, KSO.DF1DX_MAPPED)
        f_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.DF1DX_MAPPED)
        if abs(f_norm) > 1e-10:
            nabla_f *= 1.0/f_norm

        if len(g_a) == 0:
            KM.Logger.PrintInfo("ShapeOpt", "No constraints active, use negative objective gradient as search direction.")
            s = nabla_f * (-1.0)
            s *= self.step_size / f_norm
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.SEARCH_DIRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, [0.0]*len(s), KSO.CORRECTION)
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.CONTROL_POINT_UPDATE)
            return

        omega_r = KM.Matrix()
        self.optimization_utilities.AssembleBufferMatrix(omega_r, relaxation_coefficients)
        
        omega_c = KM.Vector()
        self.optimization_utilities.AssembleBufferVector(omega_c, correction_coefficients)

        KM.Logger.PrintInfo("ShapeOpt", "Assemble matrix of constraint gradient.")
        N = KM.Matrix()
        self.optimization_utilities.AssembleVectorstoMatrix(self.design_surface, N, g_a_variables)  # TODO check if gradients are 0.0! - in cpp

        settings = KM.Parameters('{ "solver_type" : "LinearSolversApplication.dense_col_piv_householder_qr" }')
        solver = dense_linear_solver_factory.ConstructSolver(settings)
        
        c = KM.Vector()
        
        KM.Logger.PrintInfo("ShapeOpt", "Calculate projected search direction and correction.")
        self.optimization_utilities.CalculateRelaxedProjectedSearchDirectionAndCorrection(
            nabla_f,
            N,
            omega_r,
            omega_c,
            solver,
            s,
            c)
        
        
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.SEARCH_DIRECTION)
        s_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.SEARCH_DIRECTION)
        s *= 1.0 / s_norm
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, s+c, KSO.CONTROL_POINT_UPDATE)
        step_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.CONTROL_POINT_UPDATE)
        
        if abs(step_norm) > 1e-10:
            s *= self.step_size / step_norm
            c *= self.step_size / step_norm
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, s+c, KSO.CONTROL_POINT_UPDATE)
        
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, s, KSO.SEARCH_DIRECTION)
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, c, KSO.CORRECTION)
        

    # --------------------------------------------------------------------------
    def __getActiveConstraints(self):
        active_constraint_values = []
        active_constraint_variables = []
        active_relaxation_coefficient = []
        active_correction_coefficient = []
        
        for constraint in self.constraints:
            if self.__isConstraintActive(constraint):
                identifier = constraint["identifier"].GetString()
                constraint_value = self.communicator.getStandardizedValue(identifier)
                buffer_value = self.constraint_buffer_variables[identifier]["buffer_value"]
                
                active_constraint_values.append(constraint_value)
                g_a_variable = self.constraint_gradient_variables[identifier]["mapped_gradient"]
                g_a_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, g_a_variable)
                g_a_variable_vector = KM.Vector()
                self.optimization_utilities.AssembleVector(self.design_surface, g_a_variable_vector, g_a_variable)
                
                if abs(g_a_norm) > 1e-10:
                    g_a_variable_vector *= 1.0/g_a_norm
                
                active_constraint_variables.append(g_a_variable_vector)
                active_relaxation_coefficient.append(min(buffer_value,1.0))
                
                maxBuffer = 2
                if buffer_value > 1.0:
                    if buffer_value < maxBuffer:
                        active_correction_coefficient.append(2*(buffer_value - 1))
                    else:
                        active_correction_coefficient.append(2*maxBuffer)
                else:
                    active_correction_coefficient.append(0.0)
                

        return active_constraint_values, active_constraint_variables, active_relaxation_coefficient, active_correction_coefficient
    
    # --------------------------------------------------------------------------
    def __isConstraintActive(self, constraint):
        identifier = constraint["identifier"].GetString()
        constraint_value = self.communicator.getStandardizedValue(identifier)
        if constraint["type"].GetString() == "=":
            return True
        elif constraint_value >= self.constraint_buffer_variables[identifier]["lower_buffer_value"]:
            return True
        else:
            return False
    # --------------------------------------------------------------------------
    def __updateBufferZone(self):
        for constraint in self.constraints:
            
            identifier = constraint["identifier"].GetString()
            constraint_value = self.communicator.getStandardizedValue(identifier)
            
            if self.optimization_iteration > 3:
                delta_g_1 = constraint_value - self.constraint_buffer_variables[identifier]["constraint_value-1"]
                delta_g_2 = self.constraint_buffer_variables[identifier]["constraint_value-1"] -self.constraint_buffer_variables[identifier]["constraint_value-2"]
                delta_g_3 = self.constraint_buffer_variables[identifier]["constraint_value-2"] -self.constraint_buffer_variables[identifier]["constraint_value-3"]
                
                if delta_g_1*delta_g_2 < 0 and delta_g_2*delta_g_3 < 0:
                    self.constraint_buffer_variables[identifier]["buffer_size_factor"] += abs(self.constraint_buffer_variables[identifier]["buffer_value"]-self.constraint_buffer_variables[identifier]["buffer_value-1"])

            if self.optimization_iteration > 1:
                delta_g = constraint_value - self.constraint_buffer_variables[identifier]["constraint_value-1"]
                
                if delta_g >= 0.0 and constraint_value > 0.0 and self.constraint_buffer_variables[identifier]["constraint_value-1"] > 0:
                    self.constraint_buffer_variables[identifier]["central_buffer_value"] -= self.constraint_buffer_variables[identifier]["constraint_value-1"]
                elif delta_g <= 0.0 and constraint_value < 0.0 and self.constraint_buffer_variables[identifier]["constraint_value-1"] < 0:
                    self.constraint_buffer_variables[identifier]["central_buffer_value"] -= self.constraint_buffer_variables[identifier]["constraint_value-1"]
                    self.constraint_buffer_variables[identifier]["central_buffer_value"] = min(self.constraint_buffer_variables[identifier]["central_buffer_value"],0.0)   
            
            self.constraint_buffer_variables[identifier]["constraint_value-3"] = self.constraint_buffer_variables[identifier]["constraint_value-2"]
            self.constraint_buffer_variables[identifier]["constraint_value-2"] = self.constraint_buffer_variables[identifier]["constraint_value-1"]
            self.constraint_buffer_variables[identifier]["constraint_value-1"] = constraint_value
    # --------------------------------------------------------------------------
    
    def __normOutputGradients(self):
        for itr, constraint in enumerate(self.constraints):
            identifier = constraint["identifier"].GetString()
            gradient = self.constraint_gradient_variables[identifier]["mapped_gradient"]
            
            gradient_vector = KM.Vector()
            self.optimization_utilities.AssembleVector(self.design_surface, gradient_vector, gradient)
            grad_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, gradient)
            
            if abs(grad_norm) > 1e-10:
                gradient_vector *= 1.0/grad_norm
                
            self.optimization_utilities.AssignVectorToVariable(self.design_surface, gradient_vector, gradient)
            
        gradient_vector = KM.Vector()
        self.optimization_utilities.AssembleVector(self.design_surface, gradient_vector, KSO.DF1DX_MAPPED)
        grad_norm = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.DF1DX_MAPPED)
        
        if abs(grad_norm) > 1e-10:
            gradient_vector *= 1.0/grad_norm
        
        self.optimization_utilities.AssignVectorToVariable(self.design_surface, gradient_vector, KSO.DF1DX_MAPPED)
    
    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.step_size
        additional_values_to_log["inf_norm_s"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.SEARCH_DIRECTION)
        additional_values_to_log["inf_norm_c"] = self.optimization_utilities.ComputeMaxNormOfNodalVariable(self.design_surface, KSO.CORRECTION)
        itr = 0
        for constraint in self.constraints:
            identifier = constraint["identifier"].GetString()
            additional_values_to_log["c"+str(itr+1)+"_buffer_value"] = self.constraint_buffer_variables[identifier]["buffer_value"]
            additional_values_to_log["c"+str(itr+1)+"_buffer_size"] = self.constraint_buffer_variables[identifier]["buffer_size"]
            additional_values_to_log["c"+str(itr+1)+"_buffer_size_factor"] = self.constraint_buffer_variables[identifier]["buffer_size_factor"]
            additional_values_to_log["c"+str(itr+1)+"_central_buffer_value"] = self.constraint_buffer_variables[identifier]["central_buffer_value"]
            additional_values_to_log["c"+str(itr+1)+"_lower_buffer_value"] = self.constraint_buffer_variables[identifier]["lower_buffer_value"]
            additional_values_to_log["c"+str(itr+1)+"_upper_buffer_value"] = self.constraint_buffer_variables[identifier]["upper_buffer_value"]
            additional_values_to_log["c"+str(itr+1)+"_laplace_multiplier"] = self.constraint_laplace_multipliers[identifier]
            itr += 1
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
