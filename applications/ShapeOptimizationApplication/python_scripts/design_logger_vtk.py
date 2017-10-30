# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#                   Sferza Massimo, https://github.com/IIIaxS
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
class DesignLoggerVTK( DesignLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, InputModelPart, OptimizationSettings ):
        self.OptimizationSettings = OptimizationSettings
        self.VtkIO = VTKFileIO( InputModelPart, OptimizationSettings )                

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        self.VtkIO.InitializeLogging()

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, optimizationIteration ):
        self.VtkIO.LogNodalResults( optimizationIteration )

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):      
        pass       

# ==============================================================================
