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

# Additional imports
import shutil
import os

from KratosMultiphysics.ShapeOptimizationApplication.loggers.design_logger_gid import DesignLoggerGID
from KratosMultiphysics.ShapeOptimizationApplication.loggers.design_logger_unv import DesignLoggerUNV
from KratosMultiphysics.ShapeOptimizationApplication.loggers.design_logger_vtk import DesignLoggerVTK

from KratosMultiphysics.ShapeOptimizationApplication.loggers.value_logger_steepest_descent import ValueLoggerSteepestDescent
from KratosMultiphysics.ShapeOptimizationApplication.loggers.value_logger_penalized_projection import ValueLoggerPenalizedProjection
from KratosMultiphysics.ShapeOptimizationApplication.loggers.value_logger_trust_region import ValueLoggerTrustRegion
from KratosMultiphysics.ShapeOptimizationApplication.loggers.value_logger_bead_optimization import ValueLoggerBeadOptimization
from KratosMultiphysics.ShapeOptimizationApplication.loggers.value_logger_gradient_projection import ValueLoggerGradientProjection

# ==============================================================================
def CreateDataLogger( ModelPartController, Communicator, OptimizationSettings ):
    return DataLogger( ModelPartController, Communicator, OptimizationSettings )

# ==============================================================================
class DataLogger():
    # --------------------------------------------------------------------------
    def __init__( self, ModelPartController, Communicator, OptimizationSettings ):
        self.ModelPartController = ModelPartController
        self.Communicator = Communicator
        self.OptimizationSettings = OptimizationSettings

        default_logger_settings = KM.Parameters("""
        {
            "output_directory"          : "Optimization_Results",
            "optimization_log_filename" : "optimization_log",
            "design_output_mode"        : "write_optimization_model_part",
            "nodal_results"             : [ "SHAPE_CHANGE" ],
            "output_format"             : { "name": "vtk" }
        }""")

        self.OptimizationSettings["output"].ValidateAndAssignDefaults(default_logger_settings)

        self.ValueLogger = self.__CreateValueLogger()
        self.DesignLogger = self.__CreateDesignLogger()

        self.__CreateFolderToStoreOptimizationResults()
        self.__OutputInformationAboutResponseFunctions()

    # -----------------------------------------------------------------------------
    def __CreateValueLogger( self ):
        AlgorithmName = self.OptimizationSettings["optimization_algorithm"]["name"].GetString()
        if AlgorithmName == "steepest_descent":
            return ValueLoggerSteepestDescent( self.Communicator, self.OptimizationSettings )
        elif AlgorithmName == "penalized_projection":
            return ValueLoggerPenalizedProjection( self.Communicator, self.OptimizationSettings )
        elif AlgorithmName == "gradient_projection":
            return ValueLoggerGradientProjection( self.Communicator, self.OptimizationSettings )
        elif AlgorithmName == "trust_region":
            return ValueLoggerTrustRegion( self.Communicator, self.OptimizationSettings )
        elif AlgorithmName == "bead_optimization":
            return ValueLoggerBeadOptimization( self.Communicator, self.OptimizationSettings )
        else:
            raise NameError("The following optimization algorithm not supported by the response logger (name may be misspelled): " + AlgorithmName)

    # -----------------------------------------------------------------------------
    def __CreateDesignLogger( self ):
        valid_output_modes = ["write_design_surface", "write_optimization_model_part", "none"]
        output_mode = self.OptimizationSettings["output"]["design_output_mode"].GetString()

        # backward compatibility
        if output_mode == "WriteDesignSurface":
            KM.Logger.PrintWarning("ShapeOpt", "'design_output_mode': 'WriteDesignSurface' is deprecated and replaced with 'write_design_surface'.")
            self.OptimizationSettings["output"]["design_output_mode"].SetString("write_design_surface")
            output_mode = self.OptimizationSettings["output"]["design_output_mode"].GetString()

        if output_mode == "WriteOptimizationModelPart":
            KM.Logger.PrintWarning("ShapeOpt", "'design_output_mode': 'WriteOptimizationModelPart' is deprecated and replaced with 'write_optimization_model_part'.")
            self.OptimizationSettings["output"]["design_output_mode"].SetString("write_optimization_model_part")
            output_mode = self.OptimizationSettings["output"]["design_output_mode"].GetString()

        if output_mode not in valid_output_modes:
            raise RuntimeError("Invalid 'design_output_mode', available options are: {}".format(valid_output_modes))

        if output_mode == "none":
            KM.Logger.Print("")
            KM.Logger.PrintInfo("ShapeOpt", "No design output will be created because 'design_output_mode' = 'None'.")
            return None

        outputFormatName = self.OptimizationSettings["output"]["output_format"]["name"].GetString()
        if outputFormatName == "gid":
            return DesignLoggerGID( self.ModelPartController, self.OptimizationSettings )
        if outputFormatName == "unv":
            return DesignLoggerUNV( self.ModelPartController, self.OptimizationSettings )
        if outputFormatName == "vtk":
            return DesignLoggerVTK( self.ModelPartController, self.OptimizationSettings )
        else:
            raise NameError("The following output format is not supported by the design logger (name may be misspelled): " + outputFormatName)

    # --------------------------------------------------------------------------
    def __CreateFolderToStoreOptimizationResults ( self ):
        resultsDirectory = self.OptimizationSettings["output"]["output_directory"].GetString()
        if os.path.exists(resultsDirectory):
            shutil.rmtree(resultsDirectory)
        os.makedirs(resultsDirectory)

    # --------------------------------------------------------------------------
    def __OutputInformationAboutResponseFunctions( self ):
        numberOfObjectives = self.OptimizationSettings["objectives"].size()
        numberOfConstraints = self.OptimizationSettings["constraints"].size()

        KM.Logger.Print("")
        KM.Logger.PrintInfo("ShapeOpt", "The following objectives are defined:\n")
        for objectiveNumber in range(numberOfObjectives):
            KM.Logger.Print(self.OptimizationSettings["objectives"][objectiveNumber])

        if numberOfConstraints != 0:
            KM.Logger.PrintInfo("ShapeOpt", "The following constraints are defined:\n")
            for constraintNumber in range(numberOfConstraints):
                KM.Logger.Print(self.OptimizationSettings["constraints"][constraintNumber],"\n")
        else:
            KM.Logger.PrintInfo("ShapeOpt", "No constraints defined.\n")

    # --------------------------------------------------------------------------
    def InitializeDataLogging( self ):
        if self.DesignLogger:
            self.DesignLogger.InitializeLogging()
        self.ValueLogger.InitializeLogging()

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, current_iteration ):
        if self.DesignLogger:
            self.DesignLogger.LogCurrentDesign( current_iteration )

    # --------------------------------------------------------------------------
    def LogCurrentValues( self, current_iteration, additional_values ):
        self.ValueLogger.LogCurrentValues( current_iteration, additional_values )

    # --------------------------------------------------------------------------
    def FinalizeDataLogging( self ):
        if self.DesignLogger:
            self.DesignLogger.FinalizeLogging()
        self.ValueLogger.FinalizeLogging()

    # --------------------------------------------------------------------------
    def GetValues( self, key ):
        return self.ValueLogger.GetValues(key)

# ==============================================================================
