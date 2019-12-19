# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO

# additional imports
from .custom_timer import Timer
from .analyzer_empty import EmptyAnalyzer
from . import model_part_controller_factory
from . import analyzer_factory
from . import communicator_factory
from . import algorithm_factory

# ==============================================================================
def CreateOptimizer(optimization_settings, model, external_analyzer=EmptyAnalyzer()):

    _ValidateSettings(optimization_settings)

    model_part_controller = model_part_controller_factory.CreateController(optimization_settings["model_settings"], model)

    analyzer = analyzer_factory.CreateAnalyzer(optimization_settings, model_part_controller, external_analyzer)

    communicator = communicator_factory.CreateCommunicator(optimization_settings)

    if optimization_settings["design_variables"]["type"].GetString() == "vertex_morphing":
        optimizer =  VertexMorphingMethod(optimization_settings, model_part_controller, analyzer, communicator)
    else:
        raise NameError("The following type of design variables is not supported by the optimizer: " + variable_type)

    return optimizer

# ------------------------------------------------------------------------------
def _ValidateSettings(optimization_settings):
    _ValidateTopLevelSettings(optimization_settings)
    _ValidateObjectiveSettingsRecursively(optimization_settings["objectives"])
    _ValidateConstraintSettings(optimization_settings["constraints"])

# ------------------------------------------------------------------------------
def _ValidateTopLevelSettings(optimization_settings):
    default_settings = KM.Parameters("""
    {
        "model_settings" : { },
        "objectives" : [ ],
        "constraints" : [ ],
        "design_variables" : { },
        "optimization_algorithm" : { },
        "output" : { }
    }""")

    for key in default_settings.keys():
        if not optimization_settings.Has(key):
            raise RuntimeError("CreateOptimizer: Required setting '{}' missing in 'optimization_settings'!".format(key))

    optimization_settings.ValidateAndAssignDefaults(default_settings)

# ------------------------------------------------------------------------------
def _ValidateObjectiveSettingsRecursively(objective_settings):
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
            _ValidateObjectiveSettingsRecursively(objective_settings[itr]["combined_responses"])

# ------------------------------------------------------------------------------
def _ValidateConstraintSettings(constraint_settings):
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


# ==============================================================================
class VertexMorphingMethod:
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, model_part_controller, analyzer, communicator):
        self.optimization_settings = optimization_settings
        self.model_part_controller = model_part_controller
        self.analyzer = analyzer
        self.communicator = communicator

        self.__AddVariablesToBeUsedByAllAglorithms()
        self.__AddVariablesToBeUsedByDesignVariables()

    # --------------------------------------------------------------------------
    def __AddVariablesToBeUsedByAllAglorithms(self):
        model_part = self.model_part_controller.GetOptimizationModelPart()
        number_of_objectives = self.optimization_settings["objectives"].size()
        number_of_constraints = self.optimization_settings["constraints"].size()

        nodal_variable = KM.KratosGlobals.GetVariable("DF1DX")
        model_part.AddNodalSolutionStepVariable(nodal_variable)
        nodal_variable = KM.KratosGlobals.GetVariable("DF1DX_MAPPED")
        model_part.AddNodalSolutionStepVariable(nodal_variable)

        for itr in range(1,number_of_constraints+1):
            nodal_variable = KM.KratosGlobals.GetVariable("DC"+str(itr)+"DX")
            model_part.AddNodalSolutionStepVariable(nodal_variable)
            nodal_variable = KM.KratosGlobals.GetVariable("DC"+str(itr)+"DX_MAPPED")
            model_part.AddNodalSolutionStepVariable(nodal_variable)

        model_part.AddNodalSolutionStepVariable(KSO.CONTROL_POINT_UPDATE)
        model_part.AddNodalSolutionStepVariable(KSO.CONTROL_POINT_CHANGE)
        model_part.AddNodalSolutionStepVariable(KSO.SHAPE_UPDATE)
        model_part.AddNodalSolutionStepVariable(KSO.SHAPE_CHANGE)
        model_part.AddNodalSolutionStepVariable(KSO.MESH_CHANGE)
        model_part.AddNodalSolutionStepVariable(KM.NORMAL)
        model_part.AddNodalSolutionStepVariable(KSO.NORMALIZED_SURFACE_NORMAL)

    def __AddVariablesToBeUsedByDesignVariables(self):
        if self.optimization_settings["design_variables"]["filter"].Has("in_plane_morphing") and \
            self.optimization_settings["design_variables"]["filter"]["in_plane_morphing"].GetBool():
                model_part = self.model_part_controller.GetOptimizationModelPart()
                model_part.AddNodalSolutionStepVariable(KSO.BACKGROUND_COORDINATE)
                model_part.AddNodalSolutionStepVariable(KSO.BACKGROUND_NORMAL)
                model_part.AddNodalSolutionStepVariable(KSO.OUT_OF_PLANE_DELTA)

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