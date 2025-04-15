# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================


# Kratos Core and Apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# additional imports
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer
from KratosMultiphysics.ShapeOptimizationApplication.analyzers.analyzer_empty import EmptyAnalyzer
from KratosMultiphysics.ShapeOptimizationApplication import model_part_controller_factory
from KratosMultiphysics.ShapeOptimizationApplication.analyzers import analyzer_factory
from KratosMultiphysics.ShapeOptimizationApplication import communicator_factory
from KratosMultiphysics.ShapeOptimizationApplication.algorithms import algorithm_factory

## Purely for backward compatibility, should be removed soon.
def CreateOptimizer(optimization_settings,model,external_analyzer=EmptyAnalyzer()):
    return Optimizer(model, optimization_settings, external_analyzer)

def Create( model, optimization_settings, external_analyzer=EmptyAnalyzer()):
    return Optimizer(model, optimization_settings, external_analyzer)

# ==============================================================================
class Optimizer:
    # --------------------------------------------------------------------------
    def __init__(self, model, optimization_settings, external_analyzer=EmptyAnalyzer()):
        self._ValidateSettings(optimization_settings)
        self.optimization_settings = optimization_settings

        self.model_part_controller = model_part_controller_factory.CreateController(optimization_settings["model_settings"], model)
        self.analyzer = analyzer_factory.CreateAnalyzer(optimization_settings, self.model_part_controller, external_analyzer)
        self.communicator = communicator_factory.CreateCommunicator(optimization_settings)

        self.design_variables = optimization_settings["design_variables"].values()
        for design_variable in self.design_variables:
            if design_variable["type"].GetString() not in ["vertex_morphing", "free_thickness", "thickness_parameter", "free_thickness_original_vm"]:
                variable_type = design_variable["type"].GetString()
                raise NameError(f"The following type of design variables is not supported by the optimizer: {variable_type}")

        self.__AddVariablesToBeUsedByAllAlgorithms()
        self.__AddVariablesToBeUsedByDesignVariables()

    # --------------------------------------------------------------------------
    def __AddVariablesToBeUsedByAllAlgorithms(self):
        model_part = self.model_part_controller.GetOptimizationModelPart()
        number_of_objectives = self.optimization_settings["objectives"].size()
        number_of_constraints = self.optimization_settings["constraints"].size()

        nodal_variable = KM.KratosGlobals.GetVariable("DF1DX")
        model_part.AddNodalSolutionStepVariable(nodal_variable)
        nodal_variable = KM.KratosGlobals.GetVariable("DF1DX_MAPPED")
        model_part.AddNodalSolutionStepVariable(nodal_variable)

        for itr in range(1,number_of_constraints+1):
            nodal_variable = KM.KratosGlobals.GetVariable(f"DC{(itr)}DX")
            model_part.AddNodalSolutionStepVariable(nodal_variable)
            nodal_variable = KM.KratosGlobals.GetVariable(f"DC{(itr)}DX_MAPPED")
            model_part.AddNodalSolutionStepVariable(nodal_variable)

        model_part.AddNodalSolutionStepVariable(KSO.CONTROL_POINT_UPDATE)
        model_part.AddNodalSolutionStepVariable(KSO.CONTROL_POINT_CHANGE)
        model_part.AddNodalSolutionStepVariable(KSO.SHAPE_UPDATE)
        model_part.AddNodalSolutionStepVariable(KSO.SHAPE_CHANGE)
        model_part.AddNodalSolutionStepVariable(KSO.MESH_CHANGE)
        model_part.AddNodalSolutionStepVariable(KM.NORMAL)
        model_part.AddNodalSolutionStepVariable(KSO.NORMALIZED_SURFACE_NORMAL)

        # variables required for remeshing
        model_part.AddNodalSolutionStepVariable(KM.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KM.DISTANCE_GRADIENT)

        # thickness variables
        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_SEARCH_DIRECTION)
        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CORRECTION)
        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_PROJECTION)

        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CONTROL)
        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CONTROL_UPDATE)
        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CONTROL_CHANGE)

        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_INITIAL)
        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS)
        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_UPDATE)
        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CHANGE)

        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CHANGE_CONTROL)
        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CHANGE_CONTROL_PROJECTED)
        model_part.AddNodalSolutionStepVariable(KSO.THICKNESS_CONTROL_PROJECTED_UPDATE)

        model_part.AddNodalSolutionStepVariable(KM.NODAL_AREA)
        nodal_variable = KM.KratosGlobals.GetVariable("DF1DT")
        model_part.AddNodalSolutionStepVariable(nodal_variable)
        nodal_variable = KM.KratosGlobals.GetVariable("DF1DT_MAPPED")
        model_part.AddNodalSolutionStepVariable(nodal_variable)
        nodal_variable = KM.KratosGlobals.GetVariable("DF1DT_PROJECTED")
        model_part.AddNodalSolutionStepVariable(nodal_variable)

        for itr in range(1,number_of_constraints+1):
            nodal_variable = KM.KratosGlobals.GetVariable(f"DC{(itr)}DT")
            model_part.AddNodalSolutionStepVariable(nodal_variable)
            nodal_variable = KM.KratosGlobals.GetVariable(f"DC{(itr)}DT_MAPPED")
            model_part.AddNodalSolutionStepVariable(nodal_variable)
            nodal_variable = KM.KratosGlobals.GetVariable(f"DC{(itr)}DT_PROJECTED")
            model_part.AddNodalSolutionStepVariable(nodal_variable)

        # sensitivity heatmap
        if self.optimization_settings["output"].Has("sensitivity_heatmap") and \
            self.optimization_settings["output"]["sensitivity_heatmap"].GetBool():
            self.__AddVariablesToBeUsedBySensitivityHeatmap()

        # water drain response
        if self.optimization_settings["objectives"][0]["identifier"].GetString() == "water_drain":
            model_part.AddNodalSolutionStepVariable(KM.NODAL_AREA)
            model_part.AddNodalSolutionStepVariable(KSO.WATER_LEVEL)
            model_part.AddNodalSolutionStepVariable(KSO.WATER_VOLUMES)


    def __AddVariablesToBeUsedByDesignVariables(self):
        for design_variable in self.design_variables:
            if design_variable.Has("filter"):
                if design_variable["filter"].Has("in_plane_morphing") and \
                    design_variable["filter"]["in_plane_morphing"].GetBool():
                        model_part = self.model_part_controller.GetOptimizationModelPart()
                        model_part.AddNodalSolutionStepVariable(KSO.BACKGROUND_COORDINATE)
                        model_part.AddNodalSolutionStepVariable(KSO.BACKGROUND_NORMAL)
                        model_part.AddNodalSolutionStepVariable(KSO.OUT_OF_PLANE_DELTA)
                if design_variable["filter"].Has("sliding_morphing") and \
                    design_variable["filter"]["sliding_morphing"].GetBool():
                        model_part = self.model_part_controller.GetOptimizationModelPart()
                        model_part.AddNodalSolutionStepVariable(KSO.BACKGROUND_COORDINATE)
                        model_part.AddNodalSolutionStepVariable(KSO.BACKGROUND_NORMAL)
                        model_part.AddNodalSolutionStepVariable(KSO.OUT_OF_PLANE_DELTA)

                if design_variable["filter"]["filter_radius"].IsString() and \
                    design_variable["filter"]["filter_radius"].GetString() == "adaptive":
                    # variables required for adaptive
                    model_part = self.model_part_controller.GetOptimizationModelPart()
                    model_part.AddNodalSolutionStepVariable(KSO.VERTEX_MORPHING_RADIUS)
                    model_part.AddNodalSolutionStepVariable(KSO.VERTEX_MORPHING_RADIUS_RAW)
                    model_part.AddNodalSolutionStepVariable(KSO.GAUSSIAN_CURVATURE)
                    model_part.AddNodalSolutionStepVariable(KSO.MAX_NEIGHBOUR_DISTANCE)

    def __AddVariablesToBeUsedBySensitivityHeatmap(self):
        model_part = self.model_part_controller.GetOptimizationModelPart()

        model_part.AddNodalSolutionStepVariable(KM.NODAL_AREA)

        sensitivity_heatmap_settings = None
        if self.optimization_settings["output"].Has("sensitivity_heatmap_settings"):
            sensitivity_heatmap_settings = self.optimization_settings["output"]["sensitivity_heatmap_settings"]

        weigthing = True
        mapping = True
        if sensitivity_heatmap_settings:
            if sensitivity_heatmap_settings.Has("sensitivity_weighting") and \
                not sensitivity_heatmap_settings["sensitivity_weighting"].GetBool():
                weigthing = False
            if sensitivity_heatmap_settings.Has("mapping") and \
                not sensitivity_heatmap_settings["mapping"].GetBool():
                mapping = False

        if self.optimization_settings["optimization_algorithm"]["name"].GetString() == "bead_optimization":
            model_part.AddNodalSolutionStepVariable(KSO.HEATMAP_DF1DALPHA)
            if weigthing:
                model_part.AddNodalSolutionStepVariable(KSO.DF1DALPHA_WEIGHTED)
                if mapping:
                    model_part.AddNodalSolutionStepVariable(KSO.DF1DALPHA_WEIGHTED_MAPPED)
        else:
            model_part.AddNodalSolutionStepVariable(KSO.HEATMAP_DF1DX)
            if weigthing:
                model_part.AddNodalSolutionStepVariable(KSO.DF1DX_WEIGHTED)
                if mapping:
                    model_part.AddNodalSolutionStepVariable(KSO.DF1DX_WEIGHTED_MAPPED)

            number_of_constraints = self.optimization_settings["constraints"].size()
            if number_of_constraints != 0:
                model_part.AddNodalSolutionStepVariable(KSO.HEATMAP_MAX)
                model_part.AddNodalSolutionStepVariable(KSO.HEATMAP_L2)
            for itr in range(1,number_of_constraints+1):
                nodal_variable = KM.KratosGlobals.GetVariable(f"HEATMAP_DC{(itr)}DX")
                model_part.AddNodalSolutionStepVariable(nodal_variable)
                if weigthing:
                    nodal_variable = KM.KratosGlobals.GetVariable(f"DC{(itr)}DX_WEIGHTED")
                    model_part.AddNodalSolutionStepVariable(nodal_variable)
                    if mapping:
                        nodal_variable = KM.KratosGlobals.GetVariable(f"DC{(itr)}DX_WEIGHTED_MAPPED")
                        model_part.AddNodalSolutionStepVariable(nodal_variable)

    # --------------------------------------------------------------------------
    def Optimize(self):
        algorithm_name = self.optimization_settings["optimization_algorithm"]["name"].GetString()

        KM.Logger.Print("")
        KM.Logger.Print("===============================================================================")
        KM.Logger.PrintInfo("ShapeOpt", Timer().GetTimeStamp(), ": Starting optimization using the following algorithm: ", algorithm_name)
        KM.Logger.Print("===============================================================================\n")

        algorithm = algorithm_factory.CreateOptimizationAlgorithm(self.optimization_settings,
                                                                  self.analyzer,
                                                                  self.communicator,
                                                                  self.model_part_controller)
        algorithm.CheckApplicability()
        algorithm.InitializeOptimizationLoop()
        algorithm.RunOptimizationLoop()
        algorithm.FinalizeOptimizationLoop()

        KM.Logger.Print("")
        KM.Logger.Print("===============================================================================")
        KM.Logger.PrintInfo("ShapeOpt", "Finished optimization")
        KM.Logger.Print("===============================================================================\n")

    # ==============================================================================
    # ------------------------------------------------------------------------------
    # ==============================================================================
    def _ValidateSettings(self, optimization_settings):
        self._ValidateTopLevelSettings(optimization_settings)
        self._ValidateObjectiveSettingsRecursively(optimization_settings["objectives"])
        self._ValidateConstraintSettings(optimization_settings["constraints"])

    # ------------------------------------------------------------------------------
    def _ValidateTopLevelSettings(self, optimization_settings):
        default_settings = KM.Parameters("""
        {
            "model_settings" : { },
            "objectives" : [ ],
            "constraints" : [ ],
            "design_variables" : [ ],
            "optimization_algorithm" : { },
            "output" : { }
        }""")

        for key in default_settings.keys():
            if not optimization_settings.Has(key):
                raise RuntimeError("Optimizer: Required setting '{}' missing in 'optimization_settings'!".format(key))

        # doing some tricks since the type of "design_variables" can be sub parameter or a list
        if optimization_settings["design_variables"].IsSubParameter():
            design_variable = optimization_settings["design_variables"].Clone()
            optimization_settings.RemoveValue("design_variables")
            optimization_settings.AddEmptyList("design_variables")
            optimization_settings["design_variables"].Append(design_variable)

        optimization_settings.ValidateAndAssignDefaults(default_settings)

    # ------------------------------------------------------------------------------
    def _ValidateObjectiveSettingsRecursively(self, objective_settings):
        default_settings = KM.Parameters("""
        {
            "identifier"                          : "NO_IDENTIFIER_SPECIFIED",
            "type"                                : "minimization",
            "scaling_factor"                      : 1.0,
            "analyzer"                            : "external",
            "response_settings"                   : {},
            "is_combined"                         : false,
            "combination_type"                    : "sum",
            "combined_responses"                  : [],
            "weight"                              : 1.0,
            "project_gradient_on_surface_normals" : false
        }""")
        for itr in range(objective_settings.size()):
            objective_settings[itr].ValidateAndAssignDefaults(default_settings)

            if objective_settings[itr]["is_combined"].GetBool():
                self._ValidateObjectiveSettingsRecursively(objective_settings[itr]["combined_responses"])

    # ------------------------------------------------------------------------------
    def _ValidateConstraintSettings(self, constraint_settings):
        default_settings = KM.Parameters("""
        {
            "identifier"                          : "NO_IDENTIFIER_SPECIFIED",
            "type"                                : "<",
            "scaling_factor"                      : 1.0,
            "reference"                           : "initial_value",
            "reference_value"                     : 1.0,
            "analyzer"                            : "external",
            "response_settings"                   : {},
            "project_gradient_on_surface_normals" : false
        }""")
        for itr in range(constraint_settings.size()):
            constraint_settings[itr].ValidateAndAssignDefaults(default_settings)
