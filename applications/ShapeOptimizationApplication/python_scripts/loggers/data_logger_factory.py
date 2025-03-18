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


# Kratos Core and Apps
import KratosMultiphysics as Kratos

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
from KratosMultiphysics.ShapeOptimizationApplication.loggers.value_logger_relaxed_gradient_projection import ValueLoggerRelaxedGradientProjection
from KratosMultiphysics.ShapeOptimizationApplication.loggers.value_logger_thickness_relaxed_gradient_projection import ValueLoggerThicknessRelaxedGradientProjection
from KratosMultiphysics.ShapeOptimizationApplication.loggers.value_logger_shape_fraction_optimization import ValueLoggerShapeFractionOptimization
from KratosMultiphysics.ShapeOptimizationApplication.loggers.sensitivity_heatmap_logger import (
    SensitivityHeatmapLoggerSteepestDescent,
    SensitivityHeatmapLoggerPenalizedProjection,
    SensitivityHeatmapLoggerGradientProjection,
    SensitivityHeatmapLoggerTrustRegion,
    SensitivityHeatmapLoggerBeadOptimization)

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

        default_logger_settings = Kratos.Parameters("""
        {
            "output_directory"          : "Optimization_Results",
            "optimization_log_filename" : "optimization_log",
            "design_output_mode"        : "write_optimization_model_part",
            "sensitivity_heatmap"       : false,
            "sensitivity_heatmap_settings": {},
            "nodal_results"             : [ "SHAPE_CHANGE" ],
            "condition_results"         : [],
            "output_format"             : { "name": "vtk" }
        }""")

        self.OptimizationSettings["output"].ValidateAndAssignDefaults(default_logger_settings)

        self.ValueLogger = self.__CreateValueLogger()
        self.DesignLogger = self.__CreateDesignLogger()
        self.SensitivityHeatmapLogger = self.__CreateSensitivityHeatmapLogger()

        self.__CreateFolderToStoreOptimizationResults()
        self.__OutputInformationAboutResponseFunctions()

    # -----------------------------------------------------------------------------
    def __CreateValueLogger( self ):
        AlgorithmName = self.OptimizationSettings["optimization_algorithm"]["name"].GetString()
        if AlgorithmName == "steepest_descent":
            return ValueLoggerSteepestDescent( self.Communicator, self.OptimizationSettings )
        elif AlgorithmName == "penalized_projection":
            return ValueLoggerPenalizedProjection( self.Communicator, self.OptimizationSettings )
        elif AlgorithmName == "gradient_projection" or AlgorithmName == "free_thickness_optimization" \
            or AlgorithmName == "free_thickness_optimization_v2" or AlgorithmName == "thickness_optimization" \
            or AlgorithmName == "free_thickness_optimization_v3":
            return ValueLoggerGradientProjection( self.Communicator, self.OptimizationSettings )
        elif AlgorithmName == "trust_region":
            return ValueLoggerTrustRegion( self.Communicator, self.OptimizationSettings )
        elif AlgorithmName == "bead_optimization":
            return ValueLoggerBeadOptimization( self.Communicator, self.OptimizationSettings )
        elif AlgorithmName == "relaxed_gradient_projection" or AlgorithmName == "free_thickness_optimization_v2_rgp" \
            or AlgorithmName == "free_thickness_optimization_v3_rgp":
            return ValueLoggerRelaxedGradientProjection(self.Communicator, self.OptimizationSettings)
        elif AlgorithmName == "free_thickness_rgp":
            return ValueLoggerThicknessRelaxedGradientProjection(self.Communicator, self.OptimizationSettings)
        elif AlgorithmName == "shape_fraction_optimization":
            return ValueLoggerShapeFractionOptimization( self.Communicator, self.OptimizationSettings )
        else:
            raise NameError("The following optimization algorithm not supported by the response logger (name may be misspelled): " + AlgorithmName)

    # -----------------------------------------------------------------------------
    def __CreateDesignLogger( self ):
        valid_output_modes = ["write_design_surface", "write_optimization_model_part", "none"]
        output_mode = self.OptimizationSettings["output"]["design_output_mode"].GetString()

        # backward compatibility
        if output_mode == "WriteDesignSurface":
            Kratos.Logger.PrintWarning("ShapeOpt", "'design_output_mode': 'WriteDesignSurface' is deprecated and replaced with 'write_design_surface'.")
            self.OptimizationSettings["output"]["design_output_mode"].SetString("write_design_surface")
            output_mode = self.OptimizationSettings["output"]["design_output_mode"].GetString()

        if output_mode == "WriteOptimizationModelPart":
            Kratos.Logger.PrintWarning("ShapeOpt", "'design_output_mode': 'WriteOptimizationModelPart' is deprecated and replaced with 'write_optimization_model_part'.")
            self.OptimizationSettings["output"]["design_output_mode"].SetString("write_optimization_model_part")
            output_mode = self.OptimizationSettings["output"]["design_output_mode"].GetString()

        if output_mode not in valid_output_modes:
            raise RuntimeError("Invalid 'design_output_mode', available options are: {}".format(valid_output_modes))

        if output_mode == "none":
            Kratos.Logger.Print("")
            Kratos.Logger.PrintInfo("ShapeOpt", "No design output will be created because 'design_output_mode' = 'None'.")
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

    # -----------------------------------------------------------------------------
    def __CreateSensitivityHeatmapLogger( self ):
        if not self.OptimizationSettings["output"]["sensitivity_heatmap"].GetBool():
            return None

        AlgorithmName = self.OptimizationSettings["optimization_algorithm"]["name"].GetString()
        if AlgorithmName == "steepest_descent":
            return SensitivityHeatmapLoggerSteepestDescent( self.ModelPartController, self.OptimizationSettings )
        elif AlgorithmName == "penalized_projection":
            return SensitivityHeatmapLoggerPenalizedProjection( self.ModelPartController, self.OptimizationSettings )
        elif AlgorithmName == "gradient_projection":
            return SensitivityHeatmapLoggerGradientProjection( self.ModelPartController, self.OptimizationSettings )
        elif AlgorithmName == "trust_region":
            return SensitivityHeatmapLoggerTrustRegion( self.ModelPartController, self.OptimizationSettings )
        elif AlgorithmName == "bead_optimization":
            return SensitivityHeatmapLoggerBeadOptimization( self.ModelPartController, self.OptimizationSettings )
        else:
            raise NameError("The following optimization algorithm is not supported by the sensitivity heatmap logger (name may be misspelled): " + AlgorithmName)

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

        Kratos.Logger.Print("")
        Kratos.Logger.PrintInfo("ShapeOpt", "The following objectives are defined:\n")
        for objectiveNumber in range(numberOfObjectives):
            Kratos.Logger.Print(self.OptimizationSettings["objectives"][objectiveNumber])

        if numberOfConstraints != 0:
            Kratos.Logger.PrintInfo("ShapeOpt", "The following constraints are defined:\n")
            for constraintNumber in range(numberOfConstraints):
                Kratos.Logger.Print(self.OptimizationSettings["constraints"][constraintNumber],"\n")
        else:
            Kratos.Logger.PrintInfo("ShapeOpt", "No constraints defined.\n")

    # --------------------------------------------------------------------------
    def InitializeDataLogging( self ):
        if self.DesignLogger:
            self.DesignLogger.InitializeLogging()
        self.ValueLogger.InitializeLogging()
        if self.SensitivityHeatmapLogger:
            self.SensitivityHeatmapLogger.InitializeLogging()

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, current_iteration ):
        if self.DesignLogger:
            self.DesignLogger.LogCurrentDesign( current_iteration )

    # --------------------------------------------------------------------------
    def LogCurrentValues( self, current_iteration, additional_values ):
        self.ValueLogger.LogCurrentValues( current_iteration, additional_values )

    # --------------------------------------------------------------------------
    def LogSensitivityHeatmap( self, current_iteration, mapper ):
        if self.SensitivityHeatmapLogger:
            self.SensitivityHeatmapLogger.LogSensitivityHeatmap( current_iteration, mapper )

    # --------------------------------------------------------------------------
    def FinalizeDataLogging( self ):
        if self.DesignLogger:
            self.DesignLogger.FinalizeLogging()
        self.ValueLogger.FinalizeLogging()

    # --------------------------------------------------------------------------
    def GetValues( self, key ):
        return self.ValueLogger.GetValues(key)

# ==============================================================================
