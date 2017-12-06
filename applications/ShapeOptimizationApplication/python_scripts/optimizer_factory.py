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
        self.__addVariablesNeededForOptimization( OptimizationModelPart )

    # --------------------------------------------------------------------------
    def __addVariablesNeededForOptimization( self, OptimizationModelPart ):
        OptimizationModelPart.AddNodalSolutionStepVariable(NORMAL)
        OptimizationModelPart.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)
        OptimizationModelPart.AddNodalSolutionStepVariable(OBJECTIVE_SENSITIVITY)
        OptimizationModelPart.AddNodalSolutionStepVariable(OBJECTIVE_SURFACE_SENSITIVITY)
        OptimizationModelPart.AddNodalSolutionStepVariable(MAPPED_OBJECTIVE_SENSITIVITY)
        OptimizationModelPart.AddNodalSolutionStepVariable(CONSTRAINT_SENSITIVITY) 
        OptimizationModelPart.AddNodalSolutionStepVariable(CONSTRAINT_SURFACE_SENSITIVITY)
        OptimizationModelPart.AddNodalSolutionStepVariable(MAPPED_CONSTRAINT_SENSITIVITY) 
        OptimizationModelPart.AddNodalSolutionStepVariable(CONTROL_POINT_UPDATE)
        OptimizationModelPart.AddNodalSolutionStepVariable(CONTROL_POINT_CHANGE)  
        OptimizationModelPart.AddNodalSolutionStepVariable(SEARCH_DIRECTION) 
        OptimizationModelPart.AddNodalSolutionStepVariable(SHAPE_UPDATE) 
        OptimizationModelPart.AddNodalSolutionStepVariable(SHAPE_CHANGE)
        OptimizationModelPart.AddNodalSolutionStepVariable(MESH_CHANGE)        

    # --------------------------------------------------------------------------
    def importModelPart( self ):
        model_part_io = ModelPartIO( self.OptimizationSettings["design_variables"]["optimization_model_part_name"].GetString() )
        model_part_io.ReadModelPart( self.OptimizationModelPart )
        buffer_size = 1
        self.OptimizationModelPart.SetBufferSize( buffer_size )
        self.OptimizationModelPart.ProcessInfo.SetValue( DOMAIN_SIZE, self.OptimizationSettings["design_variables"]["domain_size"].GetInt() )

    # --------------------------------------------------------------------------
    def importAnalyzer( self, newAnalyzer ): 
        self.analyzer = newAnalyzer

    # --------------------------------------------------------------------------
    def optimize( self ):
        
        timer = timer_factory.CreateTimer()
        algorithmName = self.OptimizationSettings["optimization_algorithm"]["name"].GetString()

        print("\n> ==============================================================================================================")
        print("> ",timer.GetTimeStamp(),": Starting optimization using the following algorithm: ", algorithmName)
        print("> ==============================================================================================================\n")
    
        algorithm = algorithm_factory.CreateAlgorithm( self.OptimizationModelPart, self.analyzer, self.OptimizationSettings )
        algorithm.execute()       

        print("\n> ==============================================================================================================")
        print("> Finished optimization                                                                                           ")
        print("> ==============================================================================================================\n")            

# ==============================================================================