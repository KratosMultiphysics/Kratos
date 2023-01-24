# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Schmölz David, https://github.com/dschmoelz
#
# ==============================================================================

# importing the Kratos Library
import KratosMultiphysics as Kratos
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_sensitivity_heatmap import ComputeSensitivityHeatmap

# ==============================================================================
class SensitivityHeatmapLoggerBase():

    # --------------------------------------------------------------------------
    def __init__( self, model_part_controller, optimization_settings ):

        default_heatmap_settings = Kratos.Parameters("""
        {
            "norm_type": "l2",
            "sensitivity_weighting": true,
            "mapping" : true,
            "relaxation_coefficient": 0.5
        }""")

        if not optimization_settings["output"].Has("sensitivity_heatmap_settings"):
            optimization_settings["output"].AddValue("sensitivity_heatmap_settings", default_heatmap_settings)

        self.heatmap_settings = optimization_settings["output"]["sensitivity_heatmap_settings"]

        # TODO: change default to reciprocal
        reciprocal = False
        if self.heatmap_settings["relaxation_coefficient"].IsString():
            if self.heatmap_settings["relaxation_coefficient"].GetString() == "reciprocal":
                reciprocal = True
                self.heatmap_settings["relaxation_coefficient"].SetDouble(-1.0)
            else:
                raise Exception("\"relaxation_coefficient\" either should be double value or \"reciprocal\".")

        self.heatmap_settings.ValidateAndAssignDefaults(default_heatmap_settings)
        if reciprocal:
            self.heatmap_settings["relaxation_coefficient"].SetString("reciprocal")

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surface = model_part_controller.GetDesignSurface()

        self.objectives = optimization_settings["objectives"]
        self.constraints = optimization_settings["constraints"]

        self.design_variable_name = None
        self.design_variable_dimension = None

    # --------------------------------------------------------------------------
    def InitializeLogging ( self ):
        raise RuntimeError("Sensitivity Heatmap logger base class is called. Please check your implementation of the function >> InitializeLogging << .")

    # --------------------------------------------------------------------------
    def LogSensitivityHeatmap(self, optimization_iteration, mapper):
        ComputeSensitivityHeatmap(self.design_surface,
                                  self.objectives, self.constraints,
                                  optimization_iteration, mapper,
                                  self.design_variable_name, self.design_variable_dimension,
                                  self.heatmap_settings)

# ==============================================================================
class SensitivityHeatmapLoggerSteepestDescent(SensitivityHeatmapLoggerBase):

    def InitializeLogging(self):
        self.design_variable_name = "X"
        self.design_variable_dimension = 3

# ==============================================================================
class SensitivityHeatmapLoggerPenalizedProjection(SensitivityHeatmapLoggerBase):

    def InitializeLogging(self):
        self.design_variable_name = "X"
        self.design_variable_dimension = 3

# ==============================================================================
class SensitivityHeatmapLoggerGradientProjection(SensitivityHeatmapLoggerBase):

    def InitializeLogging(self):
        self.design_variable_name = "X"
        self.design_variable_dimension = 3

# ==============================================================================
class SensitivityHeatmapLoggerTrustRegion(SensitivityHeatmapLoggerBase):

    def InitializeLogging(self):
        self.design_variable_name = "X"
        self.design_variable_dimension = 3

# ==============================================================================
class SensitivityHeatmapLoggerBeadOptimization(SensitivityHeatmapLoggerBase):

    def InitializeLogging(self):
        self.design_variable_name = "ALPHA"
        self.design_variable_dimension = 1
