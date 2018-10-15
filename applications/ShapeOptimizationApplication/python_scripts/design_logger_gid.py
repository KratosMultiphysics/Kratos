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

# Kratos Core and Apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# For GID output
from gid_output_process import GiDOutputProcess

# Import logger base classes
from design_logger_base import DesignLogger

# ==============================================================================
class DesignLoggerGID( DesignLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, model_part_controller, optimization_settings ):
        self.output_settings = optimization_settings["output"]
        default_gid_settings = Parameters("""
        {
            "name"              : "gid",
            "gid_configuration" : {
                "result_file_configuration" : {
                    "gidpost_flags"         : {
                        "GiDPostMode"       : "GiD_PostBinary"
                    },
                    "output_frequency": 1.0
                }
            }
        }""")
        self.output_settings["output_format"].RecursivelyValidateAndAssignDefaults(default_gid_settings)

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surface = model_part_controller.GetDesignSurface()

        self.__DetermineOutputMode()
        self.__CreateGiDIO()

    # --------------------------------------------------------------------------
    def __DetermineOutputMode( self ):
        output_mode = self.output_settings["design_output_mode"].GetString()

        self.write_design_surface = False
        self.write_optimization_model_part = False

        if output_mode == "WriteDesignSurface":
            self.write_design_surface = True
        elif output_mode == "WriteOptimizationModelPart":
            if self.optimization_model_part.NumberOfElements() == 0:
                raise NameError("Output of optimization model part in Gid-format requires definition of elements. No elements are given in current mdpa! You may change the design output mode.")
            self.write_optimization_model_part = True
        else:
            raise NameError("The following design output mode is not defined within a GiD output (name may be misspelled): " + output_mode)

    # --------------------------------------------------------------------------
    def __CreateGiDIO( self ):
        self.__ModifySettingsToMatchDefaultGiDOutputProcess()

        gid_config = self.output_settings["output_format"]["gid_configuration"]
        results_directory = self.output_settings["output_directory"].GetString()
        design_history_filename = self.output_settings["design_history_filename"].GetString()
        design_history_filename_with_path =  results_directory+"/"+design_history_filename

        if self.write_design_surface:
            self.gid_io = GiDOutputProcess(self.design_surface, design_history_filename_with_path, gid_config)
        elif self.write_optimization_model_part:
            self.gid_io = GiDOutputProcess(self.optimization_model_part, design_history_filename_with_path, gid_config)

    # --------------------------------------------------------------------------
    def __ModifySettingsToMatchDefaultGiDOutputProcess( self ):
        self.__AddNodalResultsToGidConfiguration()
        self.__SetConditionsFlagAccordingOutputMode()

    # --------------------------------------------------------------------------
    def __AddNodalResultsToGidConfiguration( self ):
        nodal_results = self.output_settings["nodal_results"]
        self.output_settings["output_format"]["gid_configuration"]["result_file_configuration"].AddValue("nodal_results", nodal_results)

    # --------------------------------------------------------------------------
    def __SetConditionsFlagAccordingOutputMode( self ):
        gid_config = self.output_settings["output_format"]["gid_configuration"]

        if not gid_config["result_file_configuration"]["gidpost_flags"].Has("WriteConditionsFlag"):
            gid_config["result_file_configuration"]["gidpost_flags"].AddEmptyValue("WriteConditionsFlag")

        if self.write_design_surface:
            gid_config["result_file_configuration"]["gidpost_flags"]["WriteConditionsFlag"].SetString("WriteConditions")
        elif self.write_optimization_model_part:
            gid_config["result_file_configuration"]["gidpost_flags"]["WriteConditionsFlag"].SetString("WriteElementsOnly")

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        self.gid_io.ExecuteInitialize()
        self.gid_io.ExecuteBeforeSolutionLoop()

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, optimizationIteration ):
        OriginalTime = self.optimization_model_part.ProcessInfo[TIME]
        self.optimization_model_part.ProcessInfo[TIME] = optimizationIteration

        self.gid_io.ExecuteInitializeSolutionStep()
        if(self.gid_io.IsOutputStep()):
            self.gid_io.PrintOutput()
        self.gid_io.ExecuteFinalizeSolutionStep()

        self.optimization_model_part.ProcessInfo[TIME] = OriginalTime

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        self.gid_io.ExecuteFinalize()

# ==============================================================================
