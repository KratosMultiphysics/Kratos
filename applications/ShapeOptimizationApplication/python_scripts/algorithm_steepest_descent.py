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
import mapper_factory
import data_logger_factory
from custom_timer import Timer
from custom_variable_utilities import WriteDictionaryDataOnNodalVariable

# ==============================================================================
class AlgorithmSteepestDescent(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = Parameters("""
        {
            "name"               : "steepest_descent",
            "max_iterations"     : 100,
            "relative_tolerance" : 1e-3,
            "gradient_tolerance" : 1e-5,
            "line_search" : {
                "line_search_type"           : "manual_stepping",
                "normalize_search_direction" : true,
                "step_size"                  : 1.0,
                "approximation_tolerance"    : 0.1,
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

        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(SEARCH_DIRECTION)

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Steepest descent algorithm only supports one objective function!")
        if self.constraints.size() > 0:
            raise RuntimeError("Steepest descent algorithm does not allow for any constraints!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.model_part_controller.ImportOptimizationModelPart()
        self.model_part_controller.InitializeMeshController()

        self.line_search_type = self.algorithm_settings["line_search"]["line_search_type"].GetString()
        self.approximation_tolerance = self.algorithm_settings["line_search"]["approximation_tolerance"].GetDouble()
        self.initial_step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.increase_factor = self.algorithm_settings["line_search"]["increase_factor"].GetDouble()
        self.max_increase_factor = self.algorithm_settings["line_search"]["max_increase_factor"].GetDouble()
        self.analyzer.InitializeBeforeOptimizationLoop()

        self.design_surface = self.model_part_controller.GetDesignSurface()

        self.mapper = mapper_factory.CreateMapper(self.design_surface, self.design_surface, self.mapper_settings)
        self.mapper.Initialize()

        self.data_logger = data_logger_factory.CreateDataLogger(self.model_part_controller, self.communicator, self.optimization_settings)
        self.data_logger.InitializeDataLogging()

        self.optimization_utilities = OptimizationUtilities(self.design_surface, self.optimization_settings)

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        for self.optimization_iteration in range(1,self.max_iterations):
            print("\n>===================================================================")
            print("> ",timer.GetTimeStamp(),": Starting optimization iteration ",self.optimization_iteration)
            print(">===================================================================\n")

            timer.StartNewLap()

            self.__initializeNewShape()

            self.__analyzeShape()

            self.__adjustStepSize()

            self.__computeShapeUpdate()

            self.__logCurrentOptimizationStep()

            print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
            print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

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
        self.model_part_controller.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
        self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyzeShape(self):
        self.communicator.initializeCommunication()
        self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
        self.communicator.requestGradientOf(self.objectives[0]["identifier"].GetString())

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.design_surface, self.optimization_iteration, self.communicator)

        objGradientDict = self.communicator.getStandardizedGradient(self.objectives[0]["identifier"].GetString())
        WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, DF1DX)

        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ComputeUnitSurfaceNormals()
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(DF1DX)

        self.model_part_controller.DampNodalVariableIfSpecified(DF1DX)

    # --------------------------------------------------------------------------
    def __adjustStepSize(self):
        if self.line_search_type == "manual_stepping":
            return

        if self.optimization_iteration > 1:
            current_step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()

            # Compare actual and estimated improvement using linear information
            dfd0 = 0.0
            for node in self.DesignSurface.Nodes:
                vec1 = node.GetSolutionStepValue(SEARCH_DIRECTION)
                vec2 = node.GetSolutionStepValue(DF1DX_MAPPED)
                dfd0 = dfd0 + vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]

            fa = self.Communicator.getStandardizedValue(self.only_obj["identifier"].GetString())
            f0 = self.previos_objective_value

            df_actual = fa - f0
            df_estimated = current_step_size*dfd0

            # Adjust step size if necessary
            if fa < f0:
                estimation_error = (df_actual-df_estimated)/df_actual

                # Increase step size if linear approximation matches the actual improvement within a specified tolerance
                if estimation_error < self.approximation_tolerance:
                    new_step_size = current_step_size * self.increase_factor
                    new_step_size = min(new_step_size, self.initial_step_size*self.max_increase_factor)

                # Leave step size unchanged if a nonliner change in the objective is observed but still a descent direction is obtained
                else:
                    new_step_size = current_step_size
            else:
                # Search approximation of optimal step using interpolation
                a = current_step_size
                corrected_step_size = - 0.5 * dfd0 * a**2 / (fa - f0 - dfd0 * a )

                # Starting from the new design, and assuming an opposite gradient direction, the step size to the approximated optimum behaves reciprocal
                new_step_size = current_step_size-corrected_step_size

            self.algorithm_settings["line_search"]["step_size"].SetDouble(new_step_size)

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):
        self.mapper.Update()
        self.mapper.InverseMap(DF1DX, DF1DX_MAPPED)

        self.optimization_utilities.ComputeSearchDirectionSteepestDescent()
        self.optimization_utilities.ComputeControlPointUpdate()

        self.mapper.Map(CONTROL_POINT_UPDATE, SHAPE_UPDATE)
        self.model_part_controller.DampNodalVariableIfSpecified(SHAPE_UPDATE)

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        self.previos_objective_value = self.Communicator.getStandardizedValue(self.only_obj["identifier"].GetString())
        self.norm_obj_gradient = self.OptimizationUtilities.ComputeL2NormOfNodalVariable(DF1DX_MAPPED)

        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        additional_values_to_log["norm_obj_gradient"] = self.norm_obj_gradient
        self.data_logger.LogCurrentValues(self.optimization_iteration, additional_values_to_log)
        self.data_logger.LogCurrentDesign(self.optimization_iteration)

    # --------------------------------------------------------------------------
    def __isAlgorithmConverged(self):

        if self.optimization_iteration > 1 :
            # Check if maximum iterations were reached
            if self.optimization_iteration == self.max_iterations:
                print("\n> Maximal iterations of optimization problem reached!")
                return True

            # Check gradient norm
            if self.optimization_iteration == 2:
                self.initial_norm_obj_gradient = self.norm_obj_gradient
            else:
                if self.norm_obj_gradient < self.gradient_tolerance*self.initial_norm_obj_gradient:
                    print("\n> Optimization problem converged as gradient norm reached specified tolerance of ",self.gradient_tolerance)
                    return True

            # Check for relative tolerance
            relativeChangeOfObjectiveValue = self.data_logger.GetValue("rel_change_obj", self.optimization_iteration)
            if abs(relativeChangeOfObjectiveValue) < self.relative_tolerance:
                print("\n> Optimization problem converged within a relative objective tolerance of ",self.relative_tolerance,"%.")
                return True

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges(self):
        self.optimization_utilities.AddFirstVariableToSecondVariable(CONTROL_POINT_UPDATE, CONTROL_POINT_CHANGE)
        self.optimization_utilities.AddFirstVariableToSecondVariable(SHAPE_UPDATE, SHAPE_CHANGE)

# ==============================================================================
