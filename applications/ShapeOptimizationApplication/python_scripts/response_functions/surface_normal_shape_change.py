# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

import time as timer
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication as KSO
from KratosMultiphysics import Logger
from KratosMultiphysics.response_functions.response_function_interface import ResponseFunctionInterface

class SurfaceNormalShapeChange(ResponseFunctionInterface):
    """
    This simple geometric response function calculates the sum of all shape updates in
    surface normal direction. Shape update in surface normal direction is considered positive,
    in the opposite direction negative.

    Note that the direction of the surface normals changes in every step,
    this is why the response value is not a meaningful measure but only an indication on how much the shape changed.

    Important settings:
    flip_normal_direction : boolean flag for changing on which side of the surface the update is considered positive.
    """

    def __init__(self, identifier, response_settings, model):
        self.identifier = identifier

        response_settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.response_settings = response_settings
        self.model = model

        self._model_part_name = response_settings["model_part_name"].GetString()
        input_type = response_settings["model_import_settings"]["input_type"].GetString()
        if input_type == "mdpa":
            self.model_part = self.model.CreateModelPart(self._model_part_name)
            domain_size = response_settings["domain_size"].GetInt()
            if domain_size not in [2, 3]:
                raise Exception("SurfaceNormalShapeChange: Invalid 'domain_size': {}".format(domain_size))
            self.model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, domain_size)
        elif input_type == "use_input_model_part":
            self.model_part = None  # will be retrieved in Initialize()
        else:
            raise Exception("SurfaceNormalShapeChange: '{}' model part input type not implemented.".format(input_type))

        self.signed_distances = None
        self.directions = None

        self.previous_value = 0.0
        self.value = 0.0

        self.gradient = {}

        self.flip_normal_direction = self.response_settings["flip_normal_direction"].GetBool()

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "response_type"         : "surface_normal_shape_change",
            "model_part_name"       : "UNKNOWN_NAME",
            "domain_size"           : 3,
            "model_import_settings" : {
                "input_type"        : "use_input_model_part",
                "input_filename"    : "UNKNOWN_NAME"
            },
            "flip_normal_direction" : false
        }""")
        return this_defaults

    def Initialize(self):
        if self.response_settings["model_import_settings"]["input_type"].GetString() == "mdpa":
            file_name = self.response_settings["model_import_settings"]["input_filename"].GetString()
            model_part_io = KM.ModelPartIO(file_name)
            model_part_io.ReadModelPart(self.model_part)
        else:
            self.model_part = self.model.GetModelPart(self._model_part_name)

    def InitializeSolutionStep(self):
        self.previous_value = self.value
        self.value = None
        self.gradient = {}
        KSO.GeometryUtilities(self.model_part).ComputeUnitSurfaceNormals()

    def CalculateValue(self):
        Logger.PrintInfo("\n> Starting primal analysis for response", self.identifier)

        startTime = timer.time()

        value = 0.0
        for node in self.model_part.Nodes:
            shape_update = node.GetSolutionStepValue(KSO.SHAPE_UPDATE)
            normalized_normal = node.GetSolutionStepValue(KSO.NORMALIZED_SURFACE_NORMAL)
            value += shape_update[0] * normalized_normal[0]
            value += shape_update[1] * normalized_normal[1]
            value += shape_update[2] * normalized_normal[2]

        self.value = value + self.previous_value

        Logger.PrintInfo("> Time needed for calculating the response value = ", round(timer.time() - startTime,2), "s")

    def CalculateGradient(self):
        Logger.PrintInfo("\n> Starting gradient calculation for response", self.identifier)

        startTime = timer.time()

        for node in self.model_part.Nodes:
            normalized_normal = node.GetSolutionStepValue(KSO.NORMALIZED_SURFACE_NORMAL)
            self.gradient[node.Id] = normalized_normal

        Logger.PrintInfo("> Time needed for calculating gradients = ", round(timer.time() - startTime,2), "s")

    def GetValue(self):
        return self.value

    def GetNodalGradient(self, variable):
        if not self.gradient:
            raise RuntimeError("Gradient was not calculated")
        if variable != KM.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        return self.gradient
