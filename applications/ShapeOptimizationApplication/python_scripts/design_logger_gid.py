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

# For GID output
from gid_output import GiDOutput

# Import logger base classes
from design_logger_base import DesignLogger

# ==============================================================================
class DesignLoggerGID( DesignLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, OptimizationModelPart, DesignSurface, OptimizationSettings ):
        self.OptimizationModelPart = OptimizationModelPart
        self.OptimizationSettings = OptimizationSettings
        self.NodalResults = self.__CreateListOfNodalResults( self.OptimizationSettings["output"] )
        self.gaussPointResults = []
        self.GidIO = self.__CreateGiDIO( OptimizationSettings )

    # ---------------------------------------------------------------------------
    def __CreateListOfNodalResults( self, outputSettings ):
        ListOfNodalResults = []
        for NodalResultNumber in range(outputSettings["nodal_results"].size()):
            ListOfNodalResults.append( outputSettings["nodal_results"][NodalResultNumber].GetString() )
        return ListOfNodalResults          

    # --------------------------------------------------------------------------
    def __CreateGiDIO( self, OptimizationSettings ):
        ResultsDirectory = OptimizationSettings["output"]["output_directory"].GetString()
        DesignHistoryFilename = OptimizationSettings["output"]["design_history_filename"].GetString()
        DesignHistoryFilenameWithPath =  ResultsDirectory+"/"+DesignHistoryFilename
        GidIO = GiDOutput( DesignHistoryFilenameWithPath,
                           OptimizationSettings["output"]["output_format"]["VolumeOutput"].GetBool(),
                           OptimizationSettings["output"]["output_format"]["GiDPostMode"].GetString(),
                           OptimizationSettings["output"]["output_format"]["GiDMultiFileFlag"].GetString(),
                           OptimizationSettings["output"]["output_format"]["GiDWriteMeshFlag"].GetBool(),
                           OptimizationSettings["output"]["output_format"]["GiDWriteConditionsFlag"].GetBool() )
        return GidIO            

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        iteratorForInitialDesign = 0
        self.GidIO.initialize_results( self.OptimizationModelPart )

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, optimizationIteration ):
        self.GidIO.write_results( optimizationIteration, self.OptimizationModelPart, self.NodalResults, self.gaussPointResults )         

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):      
        self.GidIO.finalize_results()        

# ==============================================================================
