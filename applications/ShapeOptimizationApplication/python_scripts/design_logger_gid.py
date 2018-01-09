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
        self.OutputSettings = OptimizationSettings["output"]
        
        self.__DetermineOutputMode()
        self.__CreateGiDIO()

    # --------------------------------------------------------------------------
    def __DetermineOutputMode( self ):
        OutputMode = self.OutputSettings["design_output_mode"].GetString()
        
        self.WriteDesignSurface = False
        self.WriteOptimizationModelPart = False

        if OutputMode == "WriteDesignSurface":
            self.WriteDesignSurface = True
        elif OutputMode == "WriteOptimizationModelPart":
            if self.OptimizationModelPart.NumberOfElements() == 0:
                raise NameError("Output of optimization model part in Gid-format requires definition of elements. No elements are given in current mdpa! You may change the design output mode.")              
            self.WriteOptimizationModelPart = True
        else:
            raise NameError("The following design output mode is not defined within a GiD output (name may be misspelled): " + OutputMode)              

    # --------------------------------------------------------------------------
    def __CreateGiDIO( self ):
        self.__ModifySettingsToMatchDefaultGiDOutputProcess()        

        GidConfig = self.OutputSettings["output_format"]["gid_configuration"]
        ResultsDirectory = self.OutputSettings["output_directory"].GetString()
        DesignHistoryFilename = self.OutputSettings["design_history_filename"].GetString()
        DesignHistoryFilenameWithPath =  ResultsDirectory+"/"+DesignHistoryFilename

        if self.WriteDesignSurface:     
            self.GidIO = GiDOutputProcess(self.DesignSurface, DesignHistoryFilenameWithPath, GidConfig)
        elif self.WriteOptimizationModelPart:
            self.GidIO = GiDOutputProcess(self.OptimizationModelPart, DesignHistoryFilenameWithPath, GidConfig)

    # --------------------------------------------------------------------------
    def __ModifySettingsToMatchDefaultGiDOutputProcess( self ):      
        self.__AddNodalResultsToGidConfiguration()
        self.__SetConditionsFlagAccordingOutputMode()

    # --------------------------------------------------------------------------
    def __AddNodalResultsToGidConfiguration( self ):
        NodalResults = self.OutputSettings["nodal_results"]        
        self.OutputSettings["output_format"]["gid_configuration"]["result_file_configuration"].AddValue("nodal_results", NodalResults)
        
    # --------------------------------------------------------------------------
    def __SetConditionsFlagAccordingOutputMode( self ):  
        GidConfig = self.OutputSettings["output_format"]["gid_configuration"]
        
        if not GidConfig["result_file_configuration"]["gidpost_flags"].Has("WriteConditionsFlag"):           
            GidConfig["result_file_configuration"]["gidpost_flags"].AddEmptyValue("WriteConditionsFlag")

        if self.WriteDesignSurface:     
            GidConfig["result_file_configuration"]["gidpost_flags"]["WriteConditionsFlag"].SetString("WriteConditions")
        elif self.WriteOptimizationModelPart:
            GidConfig["result_file_configuration"]["gidpost_flags"]["WriteConditionsFlag"].SetString("WriteElementsOnly")       

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        self.GidIO.ExecuteInitialize()
        self.GidIO.ExecuteBeforeSolutionLoop()        

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, optimizationIteration ):
        OriginalTime = self.OptimizationModelPart.ProcessInfo[TIME]
        self.OptimizationModelPart.ProcessInfo[TIME] = optimizationIteration

        self.GidIO.ExecuteInitializeSolutionStep()
        if(self.GidIO.IsOutputStep()):
            self.GidIO.PrintOutput()
        self.GidIO.ExecuteFinalizeSolutionStep()

        self.OptimizationModelPart.ProcessInfo[TIME] = OriginalTime
                    
    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):      
        self.GidIO.ExecuteFinalize()

# ==============================================================================
