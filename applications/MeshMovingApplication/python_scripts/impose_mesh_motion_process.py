# --- Core Imports ---
import KratosMultiphysics

# --- MeshMoving Imports ---
import KratosMultiphysics.MeshMovingApplication as MeshMovingApplication


def Factory(parameters, model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception(
            'expected input shall be a Parameters object, encapsulating a json string')
    return ImposeMeshMotionProcess(model, parameters["Parameters"])


class ImposeMeshMotionProcess(KratosMultiphysics.Process):
    """ @brief Impose a rotation followed by translation on a ModelPart.

        @details The transformation is equivalent to:
                 1) Translation to the reference frame (offset the origin)
                 2) Specified rotation
                 3) Reverse translation from the reference frame (undo origin offset)
                 4) Specified translation
                 Note: angles in radians

                 The rotation can be defined by either "euler_angles"
                 or a "rotation_axis" and "rotation_angle" pair. The following parameters can be
                 defined parametrically (see GenericFunctionUtility):
                 "euler_angles", "rotation_axis", "reference_point", "rotation_angle", "translation_vector"

                Default parameters:
                @code
                {
                    "model_part_name"       : "",
                    "interval"              : [0.0, "End"],
                    "rotation_definition"   : "rotation_axis",
                    "euler_angles"          : [0.0, 0.0, 0.0],
                    "rotation_axis"         : [0.0, 0.0, 1.0],
                    "reference_point"       : [0.0, 0.0, 0.0]
                    "rotation_angle"        : 0,
                    "translation_vector"    : [0.0, 0.0, 0.0],
                }
                @endcode

        @note the euler angles follow the convention specified by @ref{Quaternion} (Z, -X', Z")
    """

    def __init__(self,
                 model: KratosMultiphysics.Model,
                 parameters: KratosMultiphysics.Parameters):
        """Impose a rotation followed by translation on a ModelPart."""
        KratosMultiphysics.Process.__init__(self)

        # 'rotation_angle' can either be a float or "End", but the validator can't handle both
        # => convert 'rotation_angle' to the input type
        default_parameters = self.GetDefaultParameters()
        if parameters.Has("rotation_angle") and parameters["rotation_angle"].IsString():
            default_parameters.RemoveValue("rotation_angle")
            default_parameters.AddValue("rotation_angle", parameters["rotation_angle"])

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.model_part = model[parameters["model_part_name"].GetString()]

        # Parse interval
        self.interval_utility = KratosMultiphysics.IntervalUtility(parameters)

        # Determine whether a constant transform will suffice or a parametric one is needed
        requires_parametric_transform = False

        euler_angle_parameters = parameters["euler_angles"]
        rotation_axis_parameters = parameters["rotation_axis"]
        rotation_angle_parameters = parameters["rotation_angle"]
        reference_point_parameters = parameters["reference_point"]
        translation_vector_parameters = parameters["translation_vector"]

        if rotation_angle_parameters.IsString():
            requires_parametric_transform = True
        else:
            for i in range(3):
                if euler_angle_parameters[i].IsString() or rotation_axis_parameters[i].IsString() or reference_point_parameters[i].IsString() or translation_vector_parameters[i].IsString():
                    requires_parametric_transform = True
                    break

        rotation_definition = parameters["rotation_definition"].GetString()

        if requires_parametric_transform:
            if rotation_definition == "rotation_axis":
                self.transform = MeshMovingApplication.ParametricAffineTransform(
                    rotation_axis_parameters,
                    rotation_angle_parameters,
                    reference_point_parameters,
                    translation_vector_parameters)
            elif rotation_definition == "euler_angles":
                self.transform = MeshMovingApplication.ParametricAffineTransform(
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
                self.transform = MeshMovingApplication.AffineTransform(
                    rotation_axis,
                    rotation_angle,
                    reference_point,
                    translation_vector)
            elif rotation_definition == "euler_angles":
                euler_angles = euler_angle_parameters.GetVector()
                self.transform = MeshMovingApplication.AffineTransform(
                    euler_angles,
                    reference_point,
                    translation_vector)
            else:
                raise ValueError("Invalid rotation definition '{}'".format(rotation_definition))


    def ExecuteInitializeSolutionStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval_utility.IsInInterval(time):
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