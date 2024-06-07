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
try:
    import KratosMultiphysics.StructuralMechanicsApplication as KSM
except ImportError:
    KSM = None
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

        self.use_initial_normals = self.response_settings["initial_normals"].GetBool()

        self.max_update = self.response_settings["max_update"].GetDouble()

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
            "flip_normal_direction" : false,
            "initial_normals"       : false,
            "max_update"            : 1.0
        }""")
        return this_defaults

    def Initialize(self):
        if self.response_settings["model_import_settings"]["input_type"].GetString() == "mdpa":
            file_name = self.response_settings["model_import_settings"]["input_filename"].GetString()
            model_part_io = KM.ModelPartIO(file_name)
            model_part_io.ReadModelPart(self.model_part)
        else:
            self.model_part = self.model.GetModelPart(self._model_part_name)

        if self.use_initial_normals:
            KSO.GeometryUtilities(self.model_part).ComputeUnitSurfaceNormals()
            for node in self.model_part.Nodes:
                node_normal = node.GetSolutionStepValue(KSO.NORMALIZED_SURFACE_NORMAL)
                node.SetSolutionStepValue(KSO.INITIAL_NORMALIZED_SURFACE_NORMAL, 0, node_normal)

    def InitializeSolutionStep(self):
        self.previous_value = self.value
        self.value = None
        self.gradient = {}
        if not self.use_initial_normals:
            KSO.GeometryUtilities(self.model_part).ComputeUnitSurfaceNormals()

    def CalculateValue(self):
        Logger.PrintInfo("\n> Starting primal analysis for response", self.identifier)

        startTime = timer.time()

        self.value = 0.0
        value = 0.0
        max_value = 0.0
        self.active_nodes_ids = []
        for node in self.model_part.Nodes:
            value_i = self._CalculateNodeValue(node)
            max_value = max(value_i, max_value)

            if self.use_initial_normals:
                # self.value = max_value
                if value_i > 0:
                    self.value += value_i * value_i
                    self.active_nodes_ids.append(node.Id)
            else:
                value += value_i

        if not self.use_initial_normals:
            self.value = value + self.previous_value

        Logger.PrintInfo("> Time needed for calculating the response value = ", round(timer.time() - startTime,2), "s")

    def CalculateGradient(self):
        Logger.PrintInfo("\n> Starting gradient calculation for response", self.identifier)

        startTime = timer.time()

        for node in self.model_part.Nodes:
            if self.use_initial_normals:
                normalized_normal = node.GetSolutionStepValue(KSO.INITIAL_NORMALIZED_SURFACE_NORMAL)
            else:
                normalized_normal = node.GetSolutionStepValue(KSO.NORMALIZED_SURFACE_NORMAL)

            if self.flip_normal_direction:
                normalized_normal *= -1

            value_i = self._CalculateNodeValue(node)

            if self.use_initial_normals:
                if value_i > 0:
                    self.gradient[node.Id] = 2* value_i * normalized_normal
                else:
                    self.gradient[node.Id] = [0.0, 0.0, 0.0]
            else:
                self.gradient[node.Id] = normalized_normal

        Logger.PrintInfo("> Time needed for calculating gradients = ", round(timer.time() - startTime,2), "s")

    def _CalculateNodeValue(self, node):

        shape_update = node.GetSolutionStepValue(KSO.SHAPE_CHANGE)
        if self.use_initial_normals:
            normalized_normal = node.GetSolutionStepValue(KSO.INITIAL_NORMALIZED_SURFACE_NORMAL)
        else:
            normalized_normal = node.GetSolutionStepValue(KSO.NORMALIZED_SURFACE_NORMAL)

        if self.flip_normal_direction:
            normalized_normal *= -1

        value_i = shape_update[0] * normalized_normal[0]
        value_i += shape_update[1] * normalized_normal[1]
        value_i += shape_update[2] * normalized_normal[2]

        value_i -= self.max_update

        return value_i

    def GetValue(self):
        return self.value

    def GetNodalGradient(self, variable):
        if not self.gradient:
            raise RuntimeError("Gradient was not calculated")
        if variable != KM.SHAPE_SENSITIVITY:
            raise RuntimeError("GetNodalGradient: No gradient for {}!".format(variable.Name))
        return self.gradient

    def GetElementalGradient(self, variable):
        if variable != KSM.THICKNESS_SENSITIVITY:
            raise RuntimeError("GetElementalGradient: No gradient for {}!".format(variable.Name))
        gradient = {}
        for condition in self.model_part.Conditions:
            gradient[condition.Id] = 0.0
        return gradient