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
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_sensitivity_heatmap import SensitivityHeatmapCalculator

# ==============================================================================
class SensitivityHeatmapLoggerBase:

    # --------------------------------------------------------------------------
    def __init__( self, model_part_controller, optimization_settings ):

        default_heatmap_settings = Kratos.Parameters("""
        {
            "norm_type": "l2",
            "sensitivity_weighting": true,
            "mapping" : true,
            "relaxation_method": "reciprocal",
            "relaxation_coefficient": 0.5
        }""")

        self.heatmap_settings = optimization_settings["output"]["sensitivity_heatmap_settings"]
        if self.heatmap_settings.Has("relaxation_method") and self.heatmap_settings["relaxation_method"].GetString() == "constant":
            if not self.heatmap_settings.Has("relaxation_coefficient"):
                raise RuntimeError("'relaxation_coefficient' is missing for relaxation method 'constant'!")
        self.heatmap_settings.ValidateAndAssignDefaults(default_heatmap_settings)

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surface = model_part_controller.GetDesignSurface()

        self.objectives = optimization_settings["objectives"]
        self.constraints = optimization_settings["constraints"]

        self.sensitivity_heatmap_calculator = None

    # --------------------------------------------------------------------------
    def InitializeLogging ( self ):
        raise RuntimeError("Sensitivity Heatmap logger base class is called. Please check your implementation of the function >> InitializeLogging << .")

    # --------------------------------------------------------------------------
    def LogSensitivityHeatmap(self, optimization_iteration, mapper):
        self.sensitivity_heatmap_calculator.ComputeHeatmaps(optimization_iteration, mapper)

# ==============================================================================
class SensitivityHeatmapLoggerSteepestDescent(SensitivityHeatmapLoggerBase):

    def InitializeLogging(self):
        self.heatmap_settings.AddString("design_variable_name", "X")
        self.heatmap_settings.AddInt("design_variable_dimension", 3)
        self.sensitivity_heatmap_calculator = SensitivityHeatmapCalculator(self.design_surface, self.objectives, self.constraints, self.heatmap_settings)


# ==============================================================================
class SensitivityHeatmapLoggerPenalizedProjection(SensitivityHeatmapLoggerBase):

    def InitializeLogging(self):
        self.heatmap_settings.AddString("design_variable_name", "X")
        self.heatmap_settings.AddInt("design_variable_dimension", 3)
        self.sensitivity_heatmap_calculator = SensitivityHeatmapCalculator(self.design_surface, self.objectives, self.constraints, self.heatmap_settings)

# ==============================================================================
class SensitivityHeatmapLoggerGradientProjection(SensitivityHeatmapLoggerBase):

    def InitializeLogging(self):
        self.heatmap_settings.AddString("design_variable_name", "X")
        self.heatmap_settings.AddInt("design_variable_dimension", 3)
        self.sensitivity_heatmap_calculator = SensitivityHeatmapCalculator(self.design_surface, self.objectives, self.constraints, self.heatmap_settings)

# ==============================================================================
class SensitivityHeatmapLoggerTrustRegion(SensitivityHeatmapLoggerBase):

    def InitializeLogging(self):
        self.heatmap_settings.AddString("design_variable_name", "X")
        self.heatmap_settings.AddInt("design_variable_dimension", 3)
        self.sensitivity_heatmap_calculator = SensitivityHeatmapCalculator(self.design_surface, self.objectives, self.constraints, self.heatmap_settings)

# ==============================================================================
class SensitivityHeatmapLoggerBeadOptimization(SensitivityHeatmapLoggerBase):

    def InitializeLogging(self):
        self.heatmap_settings.AddString("design_variable_name", "ALPHA")
        self.heatmap_settings.AddInt("design_variable_dimension", 1)
        self.sensitivity_heatmap_calculator = SensitivityHeatmapCalculator(self.design_surface, self.objectives, self.constraints, self.heatmap_settings)
