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

# additional imports
from custom_timer import Timer
from analyzer_empty import EmptyAnalyzer
import model_part_controller_factory
import analyzer_factory
import communicator_factory
import algorithm_factory

# ==============================================================================
def CreateOptimizer(parameters, optimization_mdpa, external_analyzer=EmptyAnalyzer()):
    optimization_settings = parameters["optimization_settings"]

    model_part_controller = model_part_controller_factory.CreateController(optimization_settings, optimization_mdpa)

    analyzer = analyzer_factory.CreateAnalyzer(parameters, model_part_controller, external_analyzer)

    communicator = communicator_factory.CreateCommunicator(optimization_settings)

    if optimization_settings["design_variables"]["type"].GetString() == "vertex_morphing":
        return VertexMorphingMethod(optimization_settings, model_part_controller, analyzer, communicator)
    else:
        raise NameError("The following type of design variables is not supported by the optimizer: " + variable_type)

# ==============================================================================
class VertexMorphingMethod:
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, model_part_controller, analyzer, communicator):
        self.optimization_settings = optimization_settings
        self.model_part_controller = model_part_controller
        self.analyzer = analyzer
        self.communicator = communicator

        self.__AddNodalVariablesNeededForOptimization()

    # --------------------------------------------------------------------------
    def __AddNodalVariablesNeededForOptimization(self):
        optimization_mdpa = self.model_part_controller.GetOptimizationModelPart()
        optimization_mdpa.AddNodalSolutionStepVariable(NORMAL)
        optimization_mdpa.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)
        optimization_mdpa.AddNodalSolutionStepVariable(OBJECTIVE_SENSITIVITY)
        optimization_mdpa.AddNodalSolutionStepVariable(OBJECTIVE_SURFACE_SENSITIVITY)
        optimization_mdpa.AddNodalSolutionStepVariable(MAPPED_OBJECTIVE_SENSITIVITY)
        optimization_mdpa.AddNodalSolutionStepVariable(CONSTRAINT_SENSITIVITY)
        optimization_mdpa.AddNodalSolutionStepVariable(CONSTRAINT_SURFACE_SENSITIVITY)
        optimization_mdpa.AddNodalSolutionStepVariable(MAPPED_CONSTRAINT_SENSITIVITY)
        optimization_mdpa.AddNodalSolutionStepVariable(CONTROL_POINT_UPDATE)
        optimization_mdpa.AddNodalSolutionStepVariable(CONTROL_POINT_CHANGE)
        optimization_mdpa.AddNodalSolutionStepVariable(SEARCH_DIRECTION)
        optimization_mdpa.AddNodalSolutionStepVariable(SHAPE_UPDATE)
        optimization_mdpa.AddNodalSolutionStepVariable(SHAPE_CHANGE)
        optimization_mdpa.AddNodalSolutionStepVariable(MESH_CHANGE)

    # --------------------------------------------------------------------------
    def Optimize(self):
        algorithm_name = self.optimization_settings["optimization_algorithm"]["name"].GetString()

        print("\n> ==============================================================================================================")
        print("> ", Timer().GetTimeStamp(),": Starting optimization using the following algorithm: ", algorithm_name)
        print("> ==============================================================================================================\n")

        if self.model_part_controller.IsOptimizationModelPartAlreadyImported():
            print("> Skipping import of optimization model part as already done by another application. ")
        else:
            self.model_part_controller.ImportOptimizationModelPart()

        algorithm = algorithm_factory.CreateAlgorithm(self.optimization_settings,
                                                      self.model_part_controller,
                                                      self.analyzer,
                                                      self.communicator)

        algorithm.InitializeOptimizationLoop()
        algorithm.RunOptimizationLoop()
        algorithm.FinalizeOptimizationLoop()

        print("\n> ==============================================================================================================")
        print("> Finished optimization                                                                                           ")
        print("> ==============================================================================================================\n")

# ==============================================================================