# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
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
    def __init__( self, designSurface, optimizationSettings ):
        self.designSurface = designSurface
        self.optimizationSettings = optimizationSettings
        self.unvIO = UniversalFileIO( designSurface, optimizationSettings )                

    # --------------------------------------------------------------------------
    def initializeLogging( self ):
        iteratorForInitialDesign = 0
        self.unvIO.initializeLogging()
        self.unvIO.logCurrentDesign( iteratorForInitialDesign )

    # --------------------------------------------------------------------------
    def logCurrentDesign( self, optimizationIteration ):
        self.unvIO.logCurrentDesign( optimizationIteration )

    # --------------------------------------------------------------------------
    def finalizeLogging( self ):      
        pass       

# ==============================================================================
