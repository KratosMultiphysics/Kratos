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
from analyzer_empty import EmptyAnalyzer

# ==============================================================================
def CreateOptimizer(project_parameters, optimization_mdpa, external_analyzer=EmptyAnalyzer()):
    variable_type = project_parameters["optimization_settings"]["design_variables"]["type"].GetString()

    if variable_type == "vertex_morphing":
        return VertexMorphingMethod(project_parameters, optimization_mdpa, external_analyzer)
    else:
        raise NameError("The following type of design variables is not supported by the optimizer: " + variable_type)

# ==============================================================================
class VertexMorphingMethod:
    # --------------------------------------------------------------------------
    def __init__(self, project_parameters, optimization_mdpa, external_analyzer):
        self.project_parameters = project_parameters
        self.optimization_mdpa = optimization_mdpa
        self.external_analyzer = external_analyzer

        self.__addNodalVariablesNeededForOptimization()

    # --------------------------------------------------------------------------
    def __addNodalVariablesNeededForOptimization(self):
        self.optimization_mdpa.AddNodalSolutionStepVariable(NORMAL)
        self.optimization_mdpa.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)
        self.optimization_mdpa.AddNodalSolutionStepVariable(OBJECTIVE_SENSITIVITY)
        self.optimization_mdpa.AddNodalSolutionStepVariable(OBJECTIVE_SURFACE_SENSITIVITY)
        self.optimization_mdpa.AddNodalSolutionStepVariable(MAPPED_OBJECTIVE_SENSITIVITY)
        self.optimization_mdpa.AddNodalSolutionStepVariable(CONSTRAINT_SENSITIVITY)
        self.optimization_mdpa.AddNodalSolutionStepVariable(CONSTRAINT_SURFACE_SENSITIVITY)
        self.optimization_mdpa.AddNodalSolutionStepVariable(MAPPED_CONSTRAINT_SENSITIVITY)
        self.optimization_mdpa.AddNodalSolutionStepVariable(CONTROL_POINT_UPDATE)
        self.optimization_mdpa.AddNodalSolutionStepVariable(CONTROL_POINT_CHANGE)
        self.optimization_mdpa.AddNodalSolutionStepVariable(SEARCH_DIRECTION)
        self.optimization_mdpa.AddNodalSolutionStepVariable(SHAPE_UPDATE)
        self.optimization_mdpa.AddNodalSolutionStepVariable(SHAPE_CHANGE)
        self.optimization_mdpa.AddNodalSolutionStepVariable(MESH_CHANGE)

    # --------------------------------------------------------------------------
    def Optimize(self):
        optimization_settings = self.project_parameters["optimization_settings"]
        algorithm_name = optimization_settings["optimization_algorithm"]["name"].GetString()

        from custom_timer import Timer
        print("\n> ==============================================================================================================")
        print("> ", Timer().GetTimeStamp(),": Starting optimization using the following algorithm: ", algorithm_name)
        print("> ==============================================================================================================\n")

        import model_part_controller_factory
        mdpa_controller = model_part_controller_factory.CreateController(optimization_settings, self.optimization_mdpa)

        import analyzer_factory
        analyzer = analyzer_factory.CreateAnalyzer(self.project_parameters, mdpa_controller, self.external_analyzer)

        import communicator_factory
        communicator = communicator_factory.CreateCommunicator(optimization_settings)

        if mdpa_controller.IsOptimizationModelPartAlreadyImported():
            print("> Skipping import of optimization model part as already done by another application. ")
        else:
            mdpa_controller.ImportOptimizationModelPart()

        import algorithm_factory
        algorithm = algorithm_factory.CreateAlgorithm(optimization_settings, mdpa_controller, analyzer, communicator)

        algorithm.InitializeOptimizationLoop()
        algorithm.RunOptimizationLoop()
        algorithm.FinalizeOptimizationLoop()

        print("\n> ==============================================================================================================")
        print("> Finished optimization                                                                                           ")
        print("> ==============================================================================================================\n")

# ==============================================================================