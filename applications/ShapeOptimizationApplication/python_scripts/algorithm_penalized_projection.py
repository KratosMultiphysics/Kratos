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

# Additional imports
from .algorithm_base import OptimizationAlgorithm
from . import mapper_factory
from . import data_logger_factory
from .custom_timer import Timer
from .custom_variable_utilities import WriteDictionaryDataOnNodalVariable

# ==============================================================================
class AlgorithmPenalizedProjection(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = KM.Parameters("""
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

        self.design_surface = None
        self.mapper = None
        self.data_logger = None
        self.optimization_utilities = None

        self.objectives = optimization_settings["objectives"]
        self.constraints = optimization_settings["constraints"]

        self.step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relative_tolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.optimization_model_part.AddNodalSolutionStepVariable(KSO.SEARCH_DIRECTION)

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
        self.model_part_controller.UpdateTimeStep(self.optimization_iteration)
        self.model_part_controller.UpdateMeshAccordingInputVariable(KSO.SHAPE_UPDATE)
        self.model_part_controller.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyzeShape(self):
        self.communicator.initializeCommunication()
        self.communicator.requestValueOf(self.objectives[0]["identifier"].GetString())
        self.communicator.requestGradientOf(self.objectives[0]["identifier"].GetString())
        self.communicator.requestValueOf(self.constraints[0]["identifier"].GetString())
        self.communicator.requestGradientOf(self.constraints[0]["identifier"].GetString())

        self.analyzer.AnalyzeDesignAndReportToCommunicator(self.optimization_model_part, self.optimization_iteration, self.communicator)

        objGradientDict = self.communicator.getStandardizedGradient(self.objectives[0]["identifier"].GetString())
        conGradientDict = self.communicator.getStandardizedGradient(self.constraints[0]["identifier"].GetString())

        WriteDictionaryDataOnNodalVariable(objGradientDict, self.optimization_model_part, KSO.DF1DX)
        WriteDictionaryDataOnNodalVariable(conGradientDict, self.optimization_model_part, KSO.DC1DX)

        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool() or self.constraints[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ComputeUnitSurfaceNormals()

        if self.objectives[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(KSO.DF1DX)

        if self.constraints[0]["project_gradient_on_surface_normals"].GetBool():
            self.model_part_controller.ProjectNodalVariableOnUnitSurfaceNormals(KSO.DC1DX)

        self.model_part_controller.DampNodalVariableIfSpecified(KSO.DF1DX)
        self.model_part_controller.DampNodalVariableIfSpecified(KSO.DC1DX)

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):
        self.mapper.Update()
        self.mapper.InverseMap(KSO.DF1DX, KSO.DF1DX_MAPPED)
        self.mapper.InverseMap(KSO.DC1DX, KSO.DC1DX_MAPPED)

        constraint_value = self.communicator.getStandardizedValue(self.constraints[0]["identifier"].GetString())
        if self.__isConstraintActive(constraint_value):
            self.optimization_utilities.ComputeProjectedSearchDirection()
            self.optimization_utilities.CorrectProjectedSearchDirection(constraint_value)
        else:
            self.optimization_utilities.ComputeSearchDirectionSteepestDescent()
        self.optimization_utilities.ComputeControlPointUpdate(self.step_size)

        self.mapper.Map(KSO.CONTROL_POINT_UPDATE, KSO.SHAPE_UPDATE)
        self.model_part_controller.DampNodalVariableIfSpecified(KSO.SHAPE_UPDATE)

    # --------------------------------------------------------------------------
    def __isConstraintActive(self, constraintValue):
        if self.constraints[0]["type"].GetString() == "=":
            return True
        elif constraintValue > 0:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        additional_values_to_log = {}
        additional_values_to_log["correction_scaling"] = self.optimization_utilities.GetCorrectionScaling()
        additional_values_to_log["step_size"] = self.step_size
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
            relativeChangeOfObjectiveValue = self.data_logger.GetValues("rel_change_objective")[self.optimization_iteration]
            if abs(relativeChangeOfObjectiveValue) < self.relative_tolerance:
                print("\n> Optimization problem converged within a relative objective tolerance of ",self.relative_tolerance,"%.")
                return True

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges(self):
        self.optimization_utilities.AddFirstVariableToSecondVariable(KSO.CONTROL_POINT_UPDATE, KSO.CONTROL_POINT_CHANGE)
        self.optimization_utilities.AddFirstVariableToSecondVariable(KSO.SHAPE_UPDATE, KSO.SHAPE_CHANGE)

# ==============================================================================
