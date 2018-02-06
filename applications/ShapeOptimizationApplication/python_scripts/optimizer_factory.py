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

# Additional imports
import timer_factory
import algorithm_factory
# import analyzer_factory
import communicator_factory
import model_part_controller_factory

# ==============================================================================
def CreateOptimizer( OptimizationModelPart, OptimizationSettings ):
    design_variables_type = OptimizationSettings["design_variables"]["design_variables_type"].GetString()
    if design_variables_type == "vertex_morphing":
        return VertexMorphingMethod( OptimizationModelPart, OptimizationSettings )
    else:
        raise NameError("The following design variables type is not supported by the optimizer (name may be misspelled): " + design_variables_type)              

# ==============================================================================
class VertexMorphingMethod:
    # --------------------------------------------------------------------------
    def __init__( self, OptimizationModelPart, OptimizationSettings ):
        self.OptimizationModelPart = OptimizationModelPart
        self.OptimizationSettings = OptimizationSettings

        self.ModelPartController = model_part_controller_factory.CreateController( OptimizationModelPart, OptimizationSettings )        
        self.Communicator = communicator_factory.CreateCommunicator( OptimizationSettings )

        self.__addNodalVariablesNeededForOptimization()

    # --------------------------------------------------------------------------
    def __addNodalVariablesNeededForOptimization( self ):
        self.OptimizationModelPart.AddNodalSolutionStepVariable(NORMAL)
        self.OptimizationModelPart.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)
        self.OptimizationModelPart.AddNodalSolutionStepVariable(OBJECTIVE_SENSITIVITY)
        self.OptimizationModelPart.AddNodalSolutionStepVariable(OBJECTIVE_SURFACE_SENSITIVITY)
        self.OptimizationModelPart.AddNodalSolutionStepVariable(MAPPED_OBJECTIVE_SENSITIVITY)
        self.OptimizationModelPart.AddNodalSolutionStepVariable(CONSTRAINT_SENSITIVITY) 
        self.OptimizationModelPart.AddNodalSolutionStepVariable(CONSTRAINT_SURFACE_SENSITIVITY)
        self.OptimizationModelPart.AddNodalSolutionStepVariable(MAPPED_CONSTRAINT_SENSITIVITY) 
        self.OptimizationModelPart.AddNodalSolutionStepVariable(CONTROL_POINT_UPDATE)
        self.OptimizationModelPart.AddNodalSolutionStepVariable(CONTROL_POINT_CHANGE)  
        self.OptimizationModelPart.AddNodalSolutionStepVariable(SEARCH_DIRECTION) 
        self.OptimizationModelPart.AddNodalSolutionStepVariable(SHAPE_UPDATE) 
        self.OptimizationModelPart.AddNodalSolutionStepVariable(SHAPE_CHANGE)
        self.OptimizationModelPart.AddNodalSolutionStepVariable(MESH_CHANGE)   

    # --------------------------------------------------------------------------
    def importModelPart( self ):
        model_part_io = ModelPartIO( self.OptimizationSettings["design_variables"]["optimization_model_part_name"].GetString() )
        model_part_io.ReadModelPart( self.OptimizationModelPart )
        self.OptimizationModelPart.SetBufferSize( 1 )
        self.OptimizationModelPart.ProcessInfo.SetValue( DOMAIN_SIZE, self.OptimizationSettings["design_variables"]["domain_size"].GetInt() )

    # --------------------------------------------------------------------------
    def importAnalyzer( self, newAnalyzer ): 
        self.Analyzer = newAnalyzer

    # --------------------------------------------------------------------------
    def optimize( self ):
        timer = timer_factory.CreateTimer()
        algorithmName = self.OptimizationSettings["optimization_algorithm"]["name"].GetString()

        print("\n> ==============================================================================================================")
        print("> ",timer.GetTimeStamp(),": Starting optimization using the following algorithm: ", algorithmName)
        print("> ==============================================================================================================\n")

        algorithm = algorithm_factory.CreateAlgorithm( self.ModelPartController,
                                                       self.Analyzer, 
                                                       self.Communicator,
                                                       self.OptimizationSettings )

        algorithm.execute()       

        print("\n> ==============================================================================================================")
        print("> Finished optimization                                                                                           ")
        print("> ==============================================================================================================\n")                

# ==============================================================================