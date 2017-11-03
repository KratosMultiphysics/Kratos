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

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# Import logger base classes
from design_logger_base import DesignLogger

# ==============================================================================
class DesignLoggerUNV( DesignLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, OptimizationModelPart, DesignSurface, OptimizationSettings ):
        WriteCompleteOptimizationModelPart = self.OutputSettings["output_complete_optimization_model_part"].GetBool()
        if WriteCompleteOptimizationModelPart:
            self.UNVIO = UniversalFileIO( OptimizationModelPart, OptimizationSettings )
        else:
            self.UNVIO = UniversalFileIO( DesignSurface, OptimizationSettings )
            
    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        self.UNVIO.InitializeLogging()

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, OptimizationIteration ):
        self.UNVIO.LogNodalResults( OptimizationIteration )

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):      
        pass       

# ==============================================================================
