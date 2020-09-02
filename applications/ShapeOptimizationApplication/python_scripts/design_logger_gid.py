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
import KratosMultiphysics as KM

# For GID output
from KratosMultiphysics.gid_output_process import GiDOutputProcess

# Import logger base classes
from .design_logger_base import DesignLogger

# ==============================================================================
class DesignLoggerGID( DesignLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, model_part_controller, optimization_settings ):
        self.output_settings = optimization_settings["output"]
        minimal_gid_parameters = KM.Parameters("""
        {
            "name"           : "gid",
            "gid_parameters" : { }
        }""")

        self.output_settings["output_format"].ValidateAndAssignDefaults(minimal_gid_parameters)
        self.output_settings["output_format"]["gid_parameters"].RecursivelyValidateAndAssignDefaults(GiDOutputProcess.defaults)

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surface = model_part_controller.GetDesignSurface()

        self.write_design_surface = False
        self.write_optimization_model_part = False
        self.design_history_filename = None

        self.__DetermineOutputMode()
        self.__ModifySettingsToMatchDefaultGiDOutputProcess()
        self.__CreateGiDIO()

    # --------------------------------------------------------------------------
    def __DetermineOutputMode( self ):
        output_mode = self.output_settings["design_output_mode"].GetString()

        if output_mode == "write_design_surface":
            self.write_design_surface = True
            self.design_history_filename = self.design_surface.Name
        elif output_mode == "write_optimization_model_part":
            if self.optimization_model_part.NumberOfElements() == 0:
                raise NameError("Output of optimization model part in Gid-format requires definition of elements. No elements are given in current mdpa! You may change the design output mode.")
            self.write_optimization_model_part = True
            self.design_history_filename = self.optimization_model_part.Name
        else:
            raise NameError("The following design output mode is not defined within a GiD output (name may be misspelled): " + output_mode)

    # --------------------------------------------------------------------------
    def __ModifySettingsToMatchDefaultGiDOutputProcess( self ):
        gid_parameters = self.output_settings["output_format"]["gid_parameters"]

        # Add nodal results
        self.output_settings["output_format"]["gid_parameters"]["result_file_configuration"]["nodal_results"] = self.output_settings["nodal_results"]

        # Set condition flag
        if not gid_parameters["result_file_configuration"]["gidpost_flags"].Has("WriteConditionsFlag"):
            gid_parameters["result_file_configuration"]["gidpost_flags"].AddEmptyValue("WriteConditionsFlag")

        if self.write_design_surface:
            gid_parameters["result_file_configuration"]["gidpost_flags"]["WriteConditionsFlag"].SetString("WriteConditions")
        elif self.write_optimization_model_part:
            gid_parameters["result_file_configuration"]["gidpost_flags"]["WriteConditionsFlag"].SetString("WriteElementsOnly")

    # --------------------------------------------------------------------------
    def __CreateGiDIO( self ):

        gid_config = self.output_settings["output_format"]["gid_parameters"]
        results_directory = self.output_settings["output_directory"].GetString()
        design_history_file_path =  results_directory+"/"+self.design_history_filename

        if self.write_design_surface:
            self.gid_io = GiDOutputProcess(self.design_surface, design_history_file_path, gid_config)
        elif self.write_optimization_model_part:
            self.gid_io = GiDOutputProcess(self.optimization_model_part, design_history_file_path, gid_config)

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        self.gid_io.ExecuteInitialize()
        self.gid_io.ExecuteBeforeSolutionLoop()

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, optimizationIteration ):
        OriginalTime = self.optimization_model_part.ProcessInfo[KM.TIME]
        self.optimization_model_part.ProcessInfo[KM.TIME] = optimizationIteration

        self.gid_io.ExecuteInitializeSolutionStep()
        if(self.gid_io.IsOutputStep()):
            self.gid_io.PrintOutput()
        self.gid_io.ExecuteFinalizeSolutionStep()

        self.optimization_model_part.ProcessInfo[KM.TIME] = OriginalTime

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        self.gid_io.ExecuteFinalize()

# ==============================================================================