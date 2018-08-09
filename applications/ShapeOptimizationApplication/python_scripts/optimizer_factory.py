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
def CreateOptimizer(optimization_settings, optimization_mdpa, external_analyzer=EmptyAnalyzer()):

    model_part_controller = model_part_controller_factory.CreateController(optimization_settings, optimization_mdpa)

    analyzer = analyzer_factory.CreateAnalyzer(optimization_settings, model_part_controller, external_analyzer)

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

        self.__AddNodalVariablesNeededForOptimization( model_part_controller.GetOptimizationModelPart() )

    # --------------------------------------------------------------------------
    def __AddNodalVariablesNeededForOptimization(self, optimization_mdpa):
        number_of_objectives = self.optimization_settings["objectives"].size()
        number_of_constraints = self.optimization_settings["constraints"].size()

        for itr in range(1,number_of_objectives+1):
            nodal_variable = KratosGlobals.GetVariable("DF"+str(itr)+"DX")
            optimization_mdpa.AddNodalSolutionStepVariable(nodal_variable)
            nodal_variable = KratosGlobals.GetVariable("DF"+str(itr)+"DX_MAPPED")
            optimization_mdpa.AddNodalSolutionStepVariable(nodal_variable)

        for itr in range(1,number_of_constraints+1):
            nodal_variable = KratosGlobals.GetVariable("DC"+str(itr)+"DX")
            optimization_mdpa.AddNodalSolutionStepVariable(nodal_variable)
            nodal_variable = KratosGlobals.GetVariable("DC"+str(itr)+"DX_MAPPED")
            optimization_mdpa.AddNodalSolutionStepVariable(nodal_variable)

        optimization_mdpa.AddNodalSolutionStepVariable(CONTROL_POINT_UPDATE)
        optimization_mdpa.AddNodalSolutionStepVariable(CONTROL_POINT_CHANGE)
        optimization_mdpa.AddNodalSolutionStepVariable(SEARCH_DIRECTION)
        optimization_mdpa.AddNodalSolutionStepVariable(SHAPE_UPDATE)
        optimization_mdpa.AddNodalSolutionStepVariable(SHAPE_CHANGE)
        optimization_mdpa.AddNodalSolutionStepVariable(MESH_CHANGE)
        optimization_mdpa.AddNodalSolutionStepVariable(NORMAL)
        optimization_mdpa.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)

    # --------------------------------------------------------------------------
    def Optimize(self):
        algorithm_name = self.optimization_settings["optimization_algorithm"]["name"].GetString()

        print("\n> ==============================================================================================================")
        print("> ", Timer().GetTimeStamp(),": Starting optimization using the following algorithm: ", algorithm_name)
        print("> ==============================================================================================================\n")

        self.model_part_controller.ImportOptimizationModelPart()

        algorithm = algorithm_factory.CreateOptimizationAlgorithm(self.optimization_settings,
                                                                  self.analyzer,
                                                                  self.communicator,
                                                                  self.model_part_controller)

        algorithm.CheckApplicability()
        algorithm.InitializeOptimizationLoop()
        algorithm.RunOptimizationLoop()
        algorithm.FinalizeOptimizationLoop()

        print("\n> ==============================================================================================================")
        print("> Finished optimization                                                                                           ")
        print("> ==============================================================================================================\n")

# ==============================================================================