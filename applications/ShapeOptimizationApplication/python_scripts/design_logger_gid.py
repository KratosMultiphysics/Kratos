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
from gid_output_process import GiDOutputProcess

# Import logger base classes
from design_logger_base import DesignLogger

# ==============================================================================
class DesignLoggerGID( DesignLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, OptimizationModelPart, DesignSurface, OptimizationSettings ):
        self.OptimizationModelPart = OptimizationModelPart
        self.DesignSurface = DesignSurface
        self.OptimizationSettings = OptimizationSettings
        self.WriteCompleteOptimizationModelPart = OptimizationSettings["output"]["output_complete_optimization_model_part"].GetBool()
        
        self.__ModifySettingsToStickToDefaultGiDOutputProcess()
        self.__CreateGiDIO()

    # --------------------------------------------------------------------------
    def __ModifySettingsToStickToDefaultGiDOutputProcess( self ):      
        self.__AddNodalResultsToGidConfiguration()
        self.__EnsureConditionsAreWrittenIfOnlyDesignSurfaceShallBeOutput()

    # --------------------------------------------------------------------------
    def __AddNodalResultsToGidConfiguration( self ):
        NodalResults = self.OptimizationSettings["output"]["nodal_results"]        
        self.OptimizationSettings["output"]["output_format"]["gid_configuration"]["result_file_configuration"].AddValue("nodal_results", NodalResults)
        
    # --------------------------------------------------------------------------
    def __EnsureConditionsAreWrittenIfOnlyDesignSurfaceShallBeOutput( self ):  
        GidConfig = self.OptimizationSettings["output"]["output_format"]["gid_configuration"]
        
        ConditionSettingsAlreadyGiven = GidConfig["result_file_configuration"]["gidpost_flags"].Has("WriteConditionsFlag")
        if not ConditionSettingsAlreadyGiven:           
            GidConfig["result_file_configuration"]["gidpost_flags"].AddEmptyValue("WriteConditionsFlag")

        if self.WriteCompleteOptimizationModelPart:
            GidConfig["result_file_configuration"]["gidpost_flags"]["WriteConditionsFlag"].SetString("WriteElementsOnly")       
        else:     
            GidConfig["result_file_configuration"]["gidpost_flags"]["WriteConditionsFlag"].SetString("WriteConditions")

    # --------------------------------------------------------------------------
    def __CreateGiDIO( self ):
        ResultsDirectory = self.OptimizationSettings["output"]["output_directory"].GetString()
        DesignHistoryFilename = self.OptimizationSettings["output"]["design_history_filename"].GetString()
        DesignHistoryFilenameWithPath =  ResultsDirectory+"/"+DesignHistoryFilename
        
        GidConfig = self.OptimizationSettings["output"]["output_format"]["gid_configuration"]

        if self.WriteCompleteOptimizationModelPart:
            self.GidIO = GiDOutputProcess(self.OptimizationModelPart, DesignHistoryFilenameWithPath, GidConfig)
        else:
            self.GidIO = GiDOutputProcess(self.DesignSurface, DesignHistoryFilenameWithPath, GidConfig)

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        self.GidIO.ExecuteInitialize()
        self.GidIO.ExecuteBeforeSolutionLoop()        

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, optimizationIteration ):
        self.GidIO.ExecuteInitializeSolutionStep()
        if(self.GidIO.IsOutputStep()):
            self.GidIO.PrintOutput()
        self.GidIO.ExecuteFinalizeSolutionStep()
                    
    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):      
        self.GidIO.ExecuteFinalize()

# ==============================================================================
