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

# For GID output
from gid_output import GiDOutput

# Import logger base classes
from design_logger_base import DesignLogger

# ==============================================================================
class DesignLoggerGID( DesignLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, designSurface, optimizationSettings ):
        self.designSurface = designSurface
        self.optimizationSettings = optimizationSettings
        self.nodalResults = self.createListOfNodalResults( self.optimizationSettings["output"] )
        self.gaussPointResults = []
        self.gidIO = self.createGiDIO( optimizationSettings )

    # ---------------------------------------------------------------------------
    def createListOfNodalResults( self, outputSettings ):
        listOfNodalResults = []
        for nodalResultNumber in range(outputSettings["nodal_results"].size()):
            listOfNodalResults.append( outputSettings["nodal_results"][nodalResultNumber].GetString() )
        return listOfNodalResults          

    # --------------------------------------------------------------------------
    def createGiDIO( self, optimizationSettings ):
        resultsDirectory = optimizationSettings["output"]["output_directory"].GetString()
        designHistoryFilename = optimizationSettings["output"]["design_history_filename"].GetString()
        designHistoryFilenameWithPath =  resultsDirectory+"/"+designHistoryFilename
        gidIO = GiDOutput( designHistoryFilenameWithPath,
                           optimizationSettings["output"]["output_format"]["VolumeOutput"].GetBool(),
                           optimizationSettings["output"]["output_format"]["GiDPostMode"].GetString(),
                           optimizationSettings["output"]["output_format"]["GiDMultiFileFlag"].GetString(),
                           optimizationSettings["output"]["output_format"]["GiDWriteMeshFlag"].GetBool(),
                           optimizationSettings["output"]["output_format"]["GiDWriteConditionsFlag"].GetBool() )
        return gidIO            

    # --------------------------------------------------------------------------
    def initializeLogging( self ):
        iteratorForInitialDesign = 0
        self.gidIO.initialize_results( self.designSurface )
        self.gidIO.write_results( iteratorForInitialDesign, self.designSurface, self.nodalResults, self.gaussPointResults )           

    # --------------------------------------------------------------------------
    def logCurrentDesign( self, optimizationIteration ):
        self.gidIO.write_results( optimizationIteration, self.designSurface, self.nodalResults, self.gaussPointResults )         

    # --------------------------------------------------------------------------
    def finalizeLogging( self ):      
        self.gidIO.finalize_results()        

# ==============================================================================
