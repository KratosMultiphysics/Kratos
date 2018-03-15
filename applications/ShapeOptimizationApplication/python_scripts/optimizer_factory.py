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

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# ==============================================================================
def CreateOptimizer( project_parameters, optimization_model_part ):
    variable_type = project_parameters["optimization_settings"]["design_variables"]["type"].GetString()

    if variable_type == "vertex_morphing":
        return VertexMorphingMethod( project_parameters, optimization_model_part )
    else:
        raise NameError("The following type of design variables is not supported by the optimizer: " + variable_type)

# ==============================================================================
class VertexMorphingMethod:
    # --------------------------------------------------------------------------
    def __init__( self, project_parameters, optimization_model_part ):
        self.project_parameters = project_parameters
        self.optimization_model_part = optimization_model_part
        self.external_analyzer = None

        self.__addNodalVariablesNeededForOptimization()

    # --------------------------------------------------------------------------
    def __addNodalVariablesNeededForOptimization( self ):
        self.optimization_model_part.AddNodalSolutionStepVariable(NORMAL)
        self.optimization_model_part.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)
        self.optimization_model_part.AddNodalSolutionStepVariable(OBJECTIVE_SENSITIVITY)
        self.optimization_model_part.AddNodalSolutionStepVariable(OBJECTIVE_SURFACE_SENSITIVITY)
        self.optimization_model_part.AddNodalSolutionStepVariable(MAPPED_OBJECTIVE_SENSITIVITY)
        self.optimization_model_part.AddNodalSolutionStepVariable(CONSTRAINT_SENSITIVITY)
        self.optimization_model_part.AddNodalSolutionStepVariable(CONSTRAINT_SURFACE_SENSITIVITY)
        self.optimization_model_part.AddNodalSolutionStepVariable(MAPPED_CONSTRAINT_SENSITIVITY)
        self.optimization_model_part.AddNodalSolutionStepVariable(CONTROL_POINT_UPDATE)
        self.optimization_model_part.AddNodalSolutionStepVariable(CONTROL_POINT_CHANGE)
        self.optimization_model_part.AddNodalSolutionStepVariable(SEARCH_DIRECTION)
        self.optimization_model_part.AddNodalSolutionStepVariable(SHAPE_UPDATE)
        self.optimization_model_part.AddNodalSolutionStepVariable(SHAPE_CHANGE)
        self.optimization_model_part.AddNodalSolutionStepVariable(MESH_CHANGE)

    # --------------------------------------------------------------------------
    def importExternalAnalyzer( self, new_analyzer ):
        self.external_analyzer = new_analyzer

    # --------------------------------------------------------------------------
    def optimize( self ):
        from custom_timer import Timer
        algorithm_name = self.project_parameters["optimization_settings"]["optimization_algorithm"]["name"].GetString()

        print("\n> ==============================================================================================================")
        print("> ", Timer().GetTimeStamp(),": Starting optimization using the following algorithm: ", algorithm_name )
        print("> ==============================================================================================================\n")

        import model_part_controller_factory
        mdpa_controller = model_part_controller_factory.CreateController( self.project_parameters["optimization_settings"],
                                                                          self.optimization_model_part )

        import analyzer_factory
        analyzer = analyzer_factory.CreateAnalyzer( self.project_parameters, mdpa_controller, self.external_analyzer )

        import communicator_factory
        communicator = communicator_factory.CreateCommunicator( self.project_parameters["optimization_settings"] )

        if mdpa_controller.IsOptimizationModelPartAlreadyImported():
            print("> Skipping import of optimization model part as already done by another application. ")
        else:
            mdpa_controller.ImportOptimizationModelPart()

        import algorithm_factory
        algorithm = algorithm_factory.CreateAlgorithm( self.project_parameters["optimization_settings"],
                                                       mdpa_controller,
                                                       analyzer,
                                                       communicator )

        algorithm.InitializeOptimizationLoop()
        algorithm.RunOptimizationLoop()
        algorithm.FinalizeOptimizationLoop()

        print("\n> ==============================================================================================================")
        print("> Finished optimization                                                                                           ")
        print("> ==============================================================================================================\n")

# ==============================================================================