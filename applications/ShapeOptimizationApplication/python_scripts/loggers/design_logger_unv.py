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
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# Import logger base classes
from KratosMultiphysics.ShapeOptimizationApplication.loggers.design_logger_base import DesignLogger

# ==============================================================================
class DesignLoggerUNV( DesignLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, model_part_controller, optimization_settings ):
        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surface = model_part_controller.GetDesignSurface()
        self.output_settings = optimization_settings["output"]

        self.write_design_surface = False
        self.write_optimization_model_part = False
        self.design_history_filename = None

        self.__DetermineOutputMode()
        self.__CreateUNVIO()

    # --------------------------------------------------------------------------
    def __DetermineOutputMode( self ):
        output_mode = self.output_settings["design_output_mode"].GetString()

        if output_mode == "write_design_surface":
            self.write_design_surface = True
            self.design_history_filename = self.design_surface.Name
        elif output_mode == "write_optimization_model_part":
            if self.optimization_model_part.NumberOfElements() == 0:
                raise NameError("Output of optimization model part in UNV-format requires definition of elements. No elements are given in current mdpa! You may change the design output mode.")
            self.write_optimization_model_part = True
            self.design_history_filename = self.optimization_model_part.Name
        else:
            raise NameError("The following design output mode is not defined within a UNV output (name may be misspelled): " + OutputMode)

    # --------------------------------------------------------------------------
    def __CreateUNVIO( self ):
        results_directory = self.output_settings["output_directory"].GetString()
        design_history_file_path =  results_directory+"/"+self.design_history_filename
        nodal_results = self.output_settings["nodal_results"]

        if self.write_design_surface:
            self.UNVIO = KSO.UniversalFileIO( self.design_surface, design_history_file_path, "WriteConditionsOnly", nodal_results )
        elif self.write_optimization_model_part:
            self.UNVIO = KSO.UniversalFileIO( self.optimization_model_part, design_history_file_path, "WriteElementsOnly", nodal_results )

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
