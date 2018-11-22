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
    def __init__(self, OptimizationSettings, Analyzer, Communicator, ModelPartController):
        default_algorithm_settings = Parameters("""
        {
            "name"               : "steepest_descent",
            "max_iterations"     : 100,
            "relative_tolerance" : 1e-3,
            "line_search" : {
                "line_search_type"           : "manual_stepping",
                "normalize_search_direction" : true,
                "step_size"                  : 1.0
            }
        }""")
        self.algorithm_settings =  OptimizationSettings["optimization_algorithm"]
        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

        self.Analyzer = Analyzer
        self.Communicator = Communicator
        self.ModelPartController = ModelPartController

        self.objectives = OptimizationSettings["objectives"]
        self.constraints = OptimizationSettings["constraints"]

        self.OptimizationModelPart = ModelPartController.GetOptimizationModelPart()
        self.DesignSurface = ModelPartController.GetDesignSurface()

        self.Mapper = mapper_factory.CreateMapper(self.DesignSurface, self.DesignSurface, OptimizationSettings["design_variables"]["filter"])
        self.DataLogger = data_logger_factory.CreateDataLogger(ModelPartController, Communicator, OptimizationSettings)

        self.OptimizationUtilities = OptimizationUtilities(self.DesignSurface, OptimizationSettings)


    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if self.objectives.size() > 1:
            raise RuntimeError("Steepest descent algorithm only supports one objective function!")
        if self.constraints.size() > 0:
            raise RuntimeError("Steepest descent algorithm does not allow for any constraints!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.only_obj = self.objectives[0]

        self.maxIterations = self.algorithm_settings["max_iterations"].GetInt() + 1
        self.relativeTolerance = self.algorithm_settings["relative_tolerance"].GetDouble()

        self.ModelPartController.InitializeMeshController()
        self.Mapper.Initialize()
        self.Analyzer.InitializeBeforeOptimizationLoop()
        self.DataLogger.InitializeDataLogging()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        for self.optimization_iteration in range(1,self.maxIterations):
            print("\n>===================================================================")
            print("> ",timer.GetTimeStamp(),": Starting optimization iteration ",self.optimization_iteration)
            print(">===================================================================\n")

            timer.StartNewLap()

            self.__initializeNewShape()

            self.__analyzeShape()

            self.__computeShapeUpdate()

            self.__performLineSearchUsingThreePoints()

            self.__logCurrentOptimizationStep()

            print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
            print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

            if self.__isAlgorithmConverged():
                break
            else:
                self.__determineAbsoluteChanges()

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.DataLogger.FinalizeDataLogging()
        self.Analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __initializeNewShape(self):
        self.ModelPartController.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
        self.ModelPartController.SetReferenceMeshToMesh()

    # --------------------------------------------------------------------------
    def __analyzeShape(self):
        self.Communicator.initializeCommunication()
        self.Communicator.requestValueOf(self.only_obj["identifier"].GetString())
        self.Communicator.requestGradientOf(self.only_obj["identifier"].GetString())

        self.Analyzer.AnalyzeDesignAndReportToCommunicator(self.DesignSurface, self.optimization_iteration, self.Communicator)

        objGradientDict = self.Communicator.getStandardizedGradient(self.only_obj["identifier"].GetString())
        WriteDictionaryDataOnNodalVariable(objGradientDict, self.OptimizationModelPart, DF1DX)

        if self.only_obj["project_gradient_on_surface_normals"].GetBool():
            self.ModelPartController.ComputeUnitSurfaceNormals()
            self.ModelPartController.ProjectNodalVariableOnUnitSurfaceNormals(DF1DX)

        self.ModelPartController.DampNodalVariableIfSpecified(DF1DX)

    # --------------------------------------------------------------------------
    def __computeShapeUpdate(self):
        self.Mapper.Update()
        self.Mapper.InverseMap(DF1DX, DF1DX_MAPPED)

        self.OptimizationUtilities.ComputeSearchDirectionSteepestDescent()
        self.OptimizationUtilities.ComputeControlPointUpdate()

        self.Mapper.Map(CONTROL_POINT_UPDATE, SHAPE_UPDATE)

        self.ModelPartController.DampNodalVariableIfSpecified(SHAPE_UPDATE)

    # --------------------------------------------------------------------------
    def __performLineSearchUsingGradientAndTwoPoints(self):
        current_step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()

        f0 = self.Communicator.getStandardizedValue(self.only_obj["identifier"].GetString())

        dfda = 0.0
        for node in self.DesignSurface.Nodes:
            vec1 = 1/current_step_size * node.GetSolutionStepValue(SHAPE_UPDATE)
            vec2 = node.GetSolutionStepValue(DF1DX)
            dfda = dfda + vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]

        old_node_coordinates = []
        for node in self.OptimizationModelPart.Nodes:
            old_node_coordinates.append(node.X0)
            old_node_coordinates.append(node.Y0)
            old_node_coordinates.append(node.Z0)

        self.ModelPartController.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
        self.ModelPartController.SetReferenceMeshToMesh()

        self.Communicator.initializeCommunication()
        self.Communicator.requestValueOf(self.only_obj["identifier"].GetString())
        self.Analyzer.AnalyzeDesignAndReportToCommunicator(self.DesignSurface, self.optimization_iteration, self.Communicator)

        fa = self.Communicator.getStandardizedValue(self.only_obj["identifier"].GetString())

        is_decrease_condition_fullfilled = fa < f0 + 1e-4*current_step_size*dfda

        if is_decrease_condition_fullfilled:
            return
        else:
            a = current_step_size

            # aparab = (fa - f0 - dfda * a ) / a**2
            # bparab = dfda
            # cparab = f0

            a_optimized = - 0.5 * dfda * a**2 / (fa - f0 - dfda * a )

            # Update shape update and reset additional mesh motion
            for node in self.DesignSurface.Nodes:
                corrected_update = a_optimized*node.GetSolutionStepValue(SHAPE_UPDATE)
                node.SetSolutionStepValue(SHAPE_UPDATE, corrected_update)

        for counter, node in enumerate(self.OptimizationModelPart.Nodes):
            node.X = old_node_coordinates[3*counter+0]
            node.Y = old_node_coordinates[3*counter+1]
            node.Z = old_node_coordinates[3*counter+2]

            node.X0 = node.X
            node.Y0 = node.Y
            node.Z0 = node.Z

            # Update step size
            # new_step_size = a_optimized * current_step_size
            self.algorithm_settings["line_search"]["step_size"].SetDouble(a_optimized)

        self.Communicator.reportValue(self.only_obj["identifier"].GetString(), f0)

        # # Go half way back
        # for node in self.DesignSurface.Nodes:
        #     node.SetSolutionStepValue(SHAPE_UPDATE,-0.5*node.GetSolutionStepValue(SHAPE_UPDATE))

        # self.ModelPartController.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
        # self.ModelPartController.SetReferenceMeshToMesh()

        # self.Communicator.initializeCommunication()
        # self.Communicator.requestValueOf(self.only_obj["identifier"].GetString())
        # self.Analyzer.AnalyzeDesignAndReportToCommunicator(self.DesignSurface, self.optimization_iteration, self.Communicator)

        # f2 = self.Communicator.getStandardizedGradient(self.only_obj["identifier"].GetString())

    # --------------------------------------------------------------------------
    def __performLineSearchUsingThreePoints(self):
        current_step_size = self.algorithm_settings["line_search"]["step_size"].GetDouble()

        old_node_coordinates = []
        for node in self.OptimizationModelPart.Nodes:
            old_node_coordinates.append(node.X0)
            old_node_coordinates.append(node.Y0)
            old_node_coordinates.append(node.Z0)

        f0 = self.Communicator.getStandardizedValue(self.only_obj["identifier"].GetString())

        self.ModelPartController.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
        self.ModelPartController.SetReferenceMeshToMesh()

        self.Communicator.initializeCommunication()
        self.Communicator.requestValueOf(self.only_obj["identifier"].GetString())
        self.Analyzer.AnalyzeDesignAndReportToCommunicator(self.DesignSurface, self.optimization_iteration, self.Communicator)

        f1 = self.Communicator.getStandardizedValue(self.only_obj["identifier"].GetString())

        self.ModelPartController.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
        self.ModelPartController.SetReferenceMeshToMesh()

        self.Communicator.initializeCommunication()
        self.Communicator.requestValueOf(self.only_obj["identifier"].GetString())
        self.Analyzer.AnalyzeDesignAndReportToCommunicator(self.DesignSurface, self.optimization_iteration, self.Communicator)

        f2 = self.Communicator.getStandardizedValue(self.only_obj["identifier"].GetString())

        # Compute vertex of parabula
        x1 = 0
        x2 = 1
        x3 = 2
        y1 = f0
        y2 = f1
        y3 = f2

        denom = (x1 - x2)*(x1 - x3)*(x2 - x3)
        A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
        B = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3)) / denom
        C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom

        a_optimized = -B/(2.0*A)

        if a_optimized < 0:
            if f2<f1:
                a_optimized = 2
            else:
                a_optimized = 1
        else:
            a_optimized = min(2,a_optimized)

        # Update shape update and reset additional mesh motion
        for node in self.DesignSurface.Nodes:
            corrected_update = a_optimized*node.GetSolutionStepValue(SHAPE_UPDATE)
            node.SetSolutionStepValue(SHAPE_UPDATE, corrected_update)

        for counter, node in enumerate(self.OptimizationModelPart.Nodes):
            node.X = old_node_coordinates[3*counter+0]
            node.Y = old_node_coordinates[3*counter+1]
            node.Z = old_node_coordinates[3*counter+2]

            node.X0 = node.X
            node.Y0 = node.Y
            node.Z0 = node.Z

        # Update step size
        self.algorithm_settings["line_search"]["step_size"].SetDouble(a_optimized*current_step_size)

        self.Communicator.reportValue(self.only_obj["identifier"].GetString(), f0)

    # --------------------------------------------------------------------------
    def __logCurrentOptimizationStep(self):
        additional_values_to_log = {}
        additional_values_to_log["step_size"] = self.algorithm_settings["line_search"]["step_size"].GetDouble()
        self.DataLogger.LogCurrentValues(self.optimization_iteration, additional_values_to_log)
        self.DataLogger.LogCurrentDesign(self.optimization_iteration)

    # --------------------------------------------------------------------------
    def __isAlgorithmConverged(self):

        if self.optimization_iteration > 1 :

            # Check if maximum iterations were reached
            if self.optimization_iteration == self.maxIterations:
                print("\n> Maximal iterations of optimization problem reached!")
                return True

            # Check for relative tolerance
            relativeChangeOfObjectiveValue = self.DataLogger.GetValue("rel_change_obj", self.optimization_iteration)
            if abs(relativeChangeOfObjectiveValue) < self.relativeTolerance:
                print("\n> Optimization problem converged within a relative objective tolerance of ",self.relativeTolerance,"%.")
                return True

    # --------------------------------------------------------------------------
    def __determineAbsoluteChanges(self):
        self.OptimizationUtilities.AddFirstVariableToSecondVariable(CONTROL_POINT_UPDATE, CONTROL_POINT_CHANGE)
        self.OptimizationUtilities.AddFirstVariableToSecondVariable(SHAPE_UPDATE, SHAPE_CHANGE)

# ==============================================================================
