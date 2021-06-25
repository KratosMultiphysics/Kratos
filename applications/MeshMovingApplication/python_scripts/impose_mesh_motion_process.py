import KratosMultiphysics
import KratosMultiphysics.MeshMovingApplication as MeshMovingApplication


def Factory(parameters, model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception(
            'expected input shall be a Parameters object, encapsulating a json string')
    return ImposeMeshMotionProcess(model, parameters["Parameters"])


# Full docstring of ImposeMotionProcess (here until codacy multiline docstrings are fixed)
#
# Impose a rotation followed by translation on a ModelPart.
#
# The transformation is equivalent to:
# 1) Translation to the reference frame (offset the origin)
# 2) Specified rotation
# 3) Reverse translation from the reference frame (undo origin offset)
# 4) Specified translation
# Note: angles in radians
#
# The rotation can be defined by either "euler_angles"
# or a "rotation_axis" and "rotation_angle" pair. The following parameters can be
# defined parametrically (see GenericFunctionUtility):
# "euler_angles", "rotation_axis", "reference_point", "rotation_angle", "translation_vector"
#
# Default parameters:
# {
#     "model_part_name"       : "",
#     "interval"              : [0.0, "End"],
#     "rotation_definition"   : "rotation_axis",
#     "euler_angles"          : [0.0, 0.0, 0.0],
#     "rotation_axis"         : [0.0, 0.0, 1.0],
#     "reference_point"       : [0.0, 0.0, 0.0]
#     "rotation_angle"        : 0,
#     "translation_vector"    : [0.0, 0.0, 0.0],
# }
# Note: the euler angles follow the convention specified by @ref{Quaternion} (Z, -X', Z")


class ImposeMeshMotionProcess(KratosMultiphysics.Process):

    def __init__(self,
                 model: KratosMultiphysics.Model,
                 parameters: KratosMultiphysics.Parameters):
        """Impose a rotation followed by translation on a ModelPart."""
        KratosMultiphysics.Process.__init__(self)

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.model_part = model[parameters["model_part_name"].GetString()]

        # Parse interval (in lieu of a python interface for IntervalUtility)
        if parameters["interval"][1].IsString():
            if parameters["interval"][1].GetString() == "End":
                parameters["interval"][1].SetDouble(1e30)
            else:
                raise ValueError("interval end is expected to be either a float or 'End'")
        self.interval = parameters["interval"].GetVector()

        # Determine whether a constant transform will suffice or a parametric one is needed
        requires_parametric_transform = False

        euler_angle_parameters = parameters.GetValue("euler_angles")
        rotation_axis_parameters = parameters.GetValue("rotation_axis")
        rotation_angle_parameters = parameters.GetValue("rotation_angle")
        reference_point_parameters = parameters.GetValue("reference_point")
        translation_vector_parameters = parameters.GetValue("translation_vector")

        if rotation_angle_parameters.IsString():
            requires_parametric_transform = True
        else:
            for i in range(3):
                if euler_angle_parameters.GetArrayItem(i).IsString() or rotation_axis_parameters.GetArrayItem(i).IsString() or reference_point_parameters.GetArrayItem(i).IsString() or translation_vector_parameters.GetArrayItem(i).IsString():
                    requires_parametric_transform = True
                    break

        rotation_definition = parameters["rotation_definition"].GetString()

        if requires_parametric_transform:
            if rotation_definition == "rotation_axis":
                self.transform = MeshMovingApplication.ParametricLinearTransform(
                    rotation_axis_parameters,
                    rotation_angle_parameters,
                    reference_point_parameters,
                    translation_vector_parameters)
            elif rotation_definition == "euler_angles":
                self.transform = MeshMovingApplication.ParametricLinearTransform(
                    euler_angle_parameters,
                    reference_point_parameters,
                    translation_vector_parameters)
            else:
                raise ValueError("Invalid rotation definition '{}'".format(rotation_definition))

        else:
            reference_point = reference_point_parameters.GetVector()
            translation_vector = translation_vector_parameters.GetVector()
            if rotation_definition == "rotation_axis":
                rotation_axis = rotation_axis_parameters.GetVector()
                rotation_angle = rotation_angle_parameters.GetDouble()
                self.transform = MeshMovingApplication.LinearTransform(
                    rotation_axis,
                    rotation_angle,
                    reference_point,
                    translation_vector)
            elif rotation_definition == "euler_angles":
                euler_angles = euler_angle_parameters.GetVector()
                self.transform = MeshMovingApplication.LinearTransform(
                    euler_angles,
                    reference_point,
                    translation_vector)
            else:
                raise ValueError("Invalid rotation definition '{}'".format(rotation_definition))


    def ExecuteInitializeSolutionStep(self):
        MeshMovingApplication.MoveModelPart(self.model_part, self.transform)


    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""{
            "model_part_name"       : "",
            "interval"              : [0.0, "End"],
            "rotation_definition"   : "rotation_axis",
            "euler_angles"          : [0.0, 0.0, 0.0],
            "rotation_axis"         : [0.0, 0.0, 1.0],
            "reference_point"       : [0.0, 0.0, 0.0],
            "rotation_angle"        : 0,
            "translation_vector"    : [0.0, 0.0, 0.0]
        }""")