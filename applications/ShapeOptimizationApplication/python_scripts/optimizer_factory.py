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
def CreateOptimizer(optimization_settings, model, external_analyzer=EmptyAnalyzer()):

    default_settings = Parameters("""
    {
        "model_settings" : { },
        "objectives" : [ ],
        "constraints" : [ ],
        "design_variables" : { },
        "optimization_algorithm" : { },
        "output" : { }
    }""")

    for key in default_settings.keys():
        if not optimization_settings.Has(key):
            raise RuntimeError("CreateOptimizer: Required setting '{}' missing in 'optimization_settings'!".format(key))

    optimization_settings.ValidateAndAssignDefaults(default_settings)

    model_part_controller = model_part_controller_factory.CreateController(optimization_settings["model_settings"], model)

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

        self.__AddNodalVariablesNeededForOptimization()

    # --------------------------------------------------------------------------
    def __AddNodalVariablesNeededForOptimization(self):
        model_part = self.model_part_controller.GetOptimizationModelPart()
        number_of_objectives = self.optimization_settings["objectives"].size()
        number_of_constraints = self.optimization_settings["constraints"].size()

        for itr in range(1,number_of_objectives+1):
            nodal_variable = KratosGlobals.GetVariable("DF"+str(itr)+"DX")
            model_part.AddNodalSolutionStepVariable(nodal_variable)
            nodal_variable = KratosGlobals.GetVariable("DF"+str(itr)+"DX_MAPPED")
            model_part.AddNodalSolutionStepVariable(nodal_variable)

        for itr in range(1,number_of_constraints+1):
            nodal_variable = KratosGlobals.GetVariable("DC"+str(itr)+"DX")
            model_part.AddNodalSolutionStepVariable(nodal_variable)
            nodal_variable = KratosGlobals.GetVariable("DC"+str(itr)+"DX_MAPPED")
            model_part.AddNodalSolutionStepVariable(nodal_variable)

        model_part.AddNodalSolutionStepVariable(CONTROL_POINT_UPDATE)
        model_part.AddNodalSolutionStepVariable(CONTROL_POINT_CHANGE)
        model_part.AddNodalSolutionStepVariable(SEARCH_DIRECTION)
        model_part.AddNodalSolutionStepVariable(SHAPE_UPDATE)
        model_part.AddNodalSolutionStepVariable(SHAPE_CHANGE)
        model_part.AddNodalSolutionStepVariable(MESH_CHANGE)
        model_part.AddNodalSolutionStepVariable(NORMAL)
        model_part.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)

        # For bead optimization
        model_part.AddNodalSolutionStepVariable(ALPHA)
        model_part.AddNodalSolutionStepVariable(ALPHA_MAPPED)
        model_part.AddNodalSolutionStepVariable(DF1DALPHA)
        model_part.AddNodalSolutionStepVariable(DF1DALPHA_MAPPED)
        model_part.AddNodalSolutionStepVariable(DPDALPHA)
        model_part.AddNodalSolutionStepVariable(DPDALPHA_MAPPED)
        model_part.AddNodalSolutionStepVariable(DLDALPHA)

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