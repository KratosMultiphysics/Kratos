# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl
#                   
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication.algorithms.algorithm_base import OptimizationAlgorithm
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
from KratosMultiphysics.ShapeOptimizationApplication.loggers import data_logger_factory
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_variable_utilities import WriteDictionaryDataOnNodalVariable
import math
# ==============================================================================
class AlgorithmSteepestDescentImplicitVM(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = KM.Parameters("""
        {
            "name"               : "steepest_descent_implicit_vm",
            "max_iterations"     : 100,
            "relative_tolerance" : 1e-3,
            "gradient_tolerance" : 1e-5,
            "line_search" : {
                "line_search_type"           : "manual_stepping",
                "normalize_search_direction" : true,
                "step_size"                  : 1.0,
                "estimation_tolerance"       : 0.1,
                "increase_factor"            : 1.1,
                "max_increase_factor"        : 10.0
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

        self.previos_objective_value = None

        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()
        self.gradient_tolerance = self.algorithm_settings["gradient_tolerance"].GetDouble()
        self.line_search_type = self.algorithm_settings["line_search"]["line_search_type"].GetString()
        self.estimation_tolerance = self.algorithm_settings["line_search"]["estimation_tolerance"].GetDouble()
        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.increase_factor = self.algorithm_settings["line_search"]["increase_factor"].GetDouble()
        self.max_step_size = self.step_size*self.algorithm_settings["line_search"]["max_increase_factor"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.SEARCH_DIRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.HELMHOLTZ_VARS)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.HELMHOLTZ_SOURCE)
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.CONTROL_POINT)



    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Steepest descent algorithm only supports one objective function!")
        if self.constraints.size() > 0:
            raise RuntimeError("Steepest descent algorithm does not allow for any constraints!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.model_part_controller.Initialize()

        self.design_surface = self.model_part_controller.GetDesignSurface()

        # here we add DOFs needed by implicit vertex-morphing.
        for node in self.optimization_model_part.Nodes:
            node.AddDof(KSO.HELMHOLTZ_VARS_X)
            node.AddDof(KSO.HELMHOLTZ_VARS_Y)
            node.AddDof(KSO.HELMHOLTZ_VARS_Z)

        # here we fixed the damping/fixed regions
        for key, value in self.model_part_controller.damping_regions.items():
            for node in value.Nodes:
                node.Fix(KSO.HELMHOLTZ_VARS_X)
                node.Fix(KSO.HELMHOLTZ_VARS_Y)            
                node.Fix(KSO.HELMHOLTZ_VARS_Z)   
                node.SetSolutionStepValue(KSO.HELMHOLTZ_VARS_X,0,node.X0)
                node.SetSolutionStepValue(KSO.HELMHOLTZ_VARS_Y,0,node.Y0)            
                node.SetSolutionStepValue(KSO.HELMHOLTZ_VARS_Z,0,node.Z0)                             

        self.analyzer.InitializeBeforeOptimizationLoop()

        self.mapper = mapper_factory.CreateMapper(self.optimization_model_part, self.design_surface, self.mapper_settings)
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
            KM.Logger.PrintInfo("ShapeOpt", "",timer.GetTimeStamp(), ": Starting optimization iteration ",self.optimization_iteration)
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

    # --------------------------------------------------------------------------
    def __analyzeShape(self):
        self.communicator.initializeCommunication()
        self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
        self.communicator.requestGradientOf(self.objectives[0]["identifier"].GetString())

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.optimization_iteration, self.communicator)

        objGradientDict = self.communicator.getStandardizedGradient(self.objectives[0]["identifier"].GetString())
        WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, KSO.DF1DX)

        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ComputeUnitSurfaceNormals()
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(KSO.DF1DX)

        if self.objectives[0]["deintegrate_area"].GetBool():
            self.model_part_controller.ComputeUnitSurfaceNormals()
            self.model_part_controller.AreaDeintegrateNodalVariable(KSO.DF1DX)            

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):
        self.mapper.Update()
        self.mapper.InverseMap(KSO.DF1DX, KSO.DF1DX_MAPPED)

        # self.optimization_utilities.ComputeSearchDirectionSteepestDescent(self.design_surface)
        # normalize = self.algorithm_settings["line_search"]["normalize_search_direction"].GetBool()
        # self.optimization_utilities.ComputeControlPointUpdate(self.design_surface, self.step_size, normalize)

        self.mapper.Map(KSO.DF1DX_MAPPED, KSO.SHAPE_UPDATE)
        max_norm = 0
        for node in self.optimization_model_part.Nodes:
            node_shape_update = node.GetSolutionStepValue(KSO.SHAPE_UPDATE)
            norm = math.sqrt(node_shape_update[0]**2 + node_shape_update[1]**2 + node_shape_update[2]**2)
            if norm > max_norm :
                max_norm = norm

        for node in self.optimization_model_part.Nodes:
            node_shape_update = node.GetSolutionStepValue(KSO.SHAPE_UPDATE)  
            node.X -= self.step_size * (node_shape_update[0]/max_norm)
            node.X0 = node.X
            node.Y -= self.step_size * (node_shape_update[1]/max_norm)
            node.Y0 = node.Y
            node.Z -= self.step_size * (node_shape_update[2]/max_norm)   
            node.Z0 = node.Z 


    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        self.previos_objective_value = self.communicator.getStandardizedValue(self.objectives[0]["identifier"].GetString())
        self.norm_objective_gradient = self.optimization_utilities.ComputeL2NormOfNodalVariable(self.design_surface, KSO.DF1DX_MAPPED)

        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.step_size
        additional_values_to_log["norm_objective_gradient"] = self.norm_objective_gradient
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

            # Check gradient norm
            if self.optimization_iteration == 2:
                self.initial_norm_objective_gradient = self.norm_objective_gradient
            else:
                if self.norm_objective_gradient < self.gradient_tolerance*self.initial_norm_objective_gradient:
                    KM.Logger.Print("")
                    KM.Logger.PrintInfo("ShapeOpt", "Optimization problem converged as gradient norm reached specified tolerance of ",self.gradient_tolerance)
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
