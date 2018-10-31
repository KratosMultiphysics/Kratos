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
class AlgorithmPenalizedProjection(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = Parameters("""
        {
            "name"                    : "penalized_projection",
            "correction_scaling"      : 1.0,
            "use_adaptive_correction" : true,
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

        self.objectives = optimization_settings["objectives"]
        self.constraints = optimization_settings["constraints"]

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(SEARCH_DIRECTION)

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Penalized projection algorithm only supports one objective function!")
        if self.constraints.size() == 0:
            raise RuntimeError("Penalized projection algorithm requires definition of a constraint!")
        if self.constraints.size() > 1:
            raise RuntimeError("Penalized projection algorithm only supports one constraint!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.model_part_controller.ImportOptimizationModelPart()
        self.model_part_controller.InitializeMeshController()

        self.analyzer.InitializeBeforeOptimizationLoop()

        self.design_surface = self.model_part_controller.GetDesignSurface()

        self.mapper = mapper_factory.CreateMapper(self.design_surface, self.design_surface, self.mapper_settings)
        self.mapper.Initialize()

        self.data_logger = data_logger_factory.CreateDataLogger(self.model_part_controller, self.communicator, self.optimization_settings)
        self.data_logger.InitializeDataLogging()

        self.optimization_utilities = OptimizationUtilities(self.design_surface, self.optimization_settings)

        self.only_obj = self.objectives[0]
        self.only_con = self.constraints[0]

        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        for self.optimization_iteration in range(1,self.max_iterations):
            print("\n>===================================================================")
            print("> ",timer.GetTimeStamp(),": Starting optimization iteration ", self.optimization_iteration)
            print(">===================================================================\n")

            timer.StartNewLap()

            self.__initializeNewShape()

            self.__analyzeShape()

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
        self.model_part_controller.InitializeNewOptimizationStep(self.optimization_iteration)
        self.model_part_controller.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
        self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyzeShape(self):
        self.communicator.initializeCommunication()
        self.communicator.requestValueOf(self.only_obj["identifier"].GetString())
        self.communicator.requestGradientOf(self.only_obj["identifier"].GetString())
        self.communicator.requestValueOf(self.only_con["identifier"].GetString())
        self.communicator.requestGradientOf(self.only_con["identifier"].GetString())

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.design_surface, self.optimization_iteration, self.communicator)

        objGradientDict = self.communicator.getStandardizedGradient(self.only_obj["identifier"].GetString())
        conGradientDict = self.communicator.getStandardizedGradient(self.only_con["identifier"].GetString())

        WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, DF1DX)
        WriteDictionaryDataOnNodalVariable(conGradientDict, self.optimization_model_part, DC1DX)

        if self.only_obj["project_gradient_on_surface_normals"].GetBool() or self.only_con["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ComputeUnitSurfaceNormals()

        if self.only_obj["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(DF1DX)

        if self.only_con["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(DC1DX)

        self.model_part_controller.DampNodalVariableIfSpecified(DF1DX)
        self.model_part_controller.DampNodalVariableIfSpecified(DC1DX)

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):
        self.mapper.Update()
        self.mapper.InverseMap(DF1DX, DF1DX_MAPPED)
        self.mapper.InverseMap(DC1DX, DC1DX_MAPPED)

        constraint_value = self.communicator.getStandardizedValue(self.only_con["identifier"].GetString())
        if self.__isConstraintActive(constraint_value):
            self.optimization_utilities.ComputeProjectedSearchDirection()
            self.optimization_utilities.CorrectProjectedSearchDirection(constraint_value)
        else:
            self.optimization_utilities.ComputeSearchDirectionSteepestDescent()
        self.optimization_utilities.ComputeControlPointUpdate()

        self.mapper.Map(CONTROL_POINT_UPDATE, SHAPE_UPDATE)
        self.model_part_controller.DampNodalVariableIfSpecified(SHAPE_UPDATE)

    # --------------------------------------------------------------------------
    def __isConstraintActive(self, constraintValue):
        if self.only_con["type"].GetString() == "=":
            return True
        elif constraintValue > 0:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        additional_values_to_log["correction_scaling"] = self.algorithm_settings["correction_scaling"].GetDouble()
        self.data_logger.LogCurrentValues(self.optimization_iteration, additional_values_to_log)
        self.data_logger.LogCurrentDesign(self.optimization_iteration)

    # --------------------------------------------------------------------------
    def __isAlgorithmConverged(self):

        if self.optimization_iteration > 1 :

            # Check if maximum iterations were reached
            if self.optimization_iteration == self.max_iterations:
                print("\n> Maximal iterations of optimization problem reached!")
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
