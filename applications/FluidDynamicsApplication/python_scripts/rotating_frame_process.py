import KratosMultiphysics as KM
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import numpy as np

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return RotatingFrameProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class RotatingFrameProcess(KM.Process):
    """A Kratos process that rotates a given RotatingFrameModelPart 
    around a specified axis of rotation. This process uses the Kratos 
    RotateNodesUtility to rotate the nodes of the RotatingFrameModelPart 
    and assigns the corresponding MeshDisplacement to the rotated nodes. 
    The process also assigns the corresponding rotational mesh velocity to 
    the nodes in a RotatingObjectModelPart. The process takes as input the 
    model part to be rotated, the axis of rotation vector, the final angular 
    velocity in radians per second, and the accelerating time in seconds.

    Public member variables:
    model -- the container of the different model parts.
    settings -- Kratos parameters containing process settings.
    """

    def __init__(self, model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        KM.Process.__init__(self) # calling the baseclass constructor

        default_settings = KM.Parameters("""{
            "rotating_frame_model_part_name": "",
            "rotating_object_model_part_name": "",
            "center_of_rotation": [0.0,0.0,0.0],
            "axis_of_rotation": [1.0,0.0,0.0],
            "target_angular_velocity_radians": 0.0,
            "acceleration_time": 0.0
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        settings.ValidateAndAssignDefaults(default_settings)
        # Alternative:
        # settings.RecursivelyValidateAndAssignDefaults(default_settings)

        #Assign settings
        self.model = model

        # Get the rotating frame model part name
        if not settings["rotating_frame_model_part_name"].GetString():
            raise Exception("\'rotating_frame_model_part_name\' not provided. Please specify the model part to rotate the frame.")
        self.rotating_frame_model_part_name = settings["rotating_frame_model_part_name"].GetString() 
        
        # Get the rotating object model part name
        if not settings["rotating_object_model_part_name"].GetString():
            raise Exception("\'rotating_object_model_part_name\' not provided. Please specify the slave model part.")
        self.rotating_object_model_part_name = settings["rotating_object_model_part_name"].GetString()

        ## Assign rotating frame and object model parts
        self.rotating_frame_model_part = self.model.GetModelPart(self.rotating_frame_model_part_name)
        self.rotating_object_model_part = self.model.GetModelPart(self.rotating_object_model_part_name)

        # Get the center of rotation
        if settings.Has("center_of_rotation"):
            self.center_of_rotation = np.array(settings["center_of_rotation"].GetVector())
        else:
            raise Exception("The center_of_rotation parameter is missing from the settings.")

        # Get the axis of rotation
        if settings.Has("axis_of_rotation"):
            self.axis_of_rotation = np.array(settings["axis_of_rotation"].GetVector())
            if self.axis_of_rotation.size == 0:
                raise Exception("The axis_of_rotation vector is empty.")
            axis_norm = np.linalg.norm(self.axis_of_rotation)
            if not np.isclose(np.linalg.norm(axis_norm), 1.0, rtol=1e-6):
                KM.Logger.PrintWarning("RotatingFrameProcess", "The axis_of_rotation vector is not a unit vector... normalizing the vector.")
                self.axis_of_rotation = self.axis_of_rotation / axis_norm
        else:
            raise Exception("The axis_of_rotation parameter is missing from the settings.")

        # Get the angular velocity in radians
        if settings.Has("target_angular_velocity_radians"):
            self.target_angular_velocity_radians = settings["target_angular_velocity_radians"].GetDouble()
        else:
            raise Exception("The target_angular_velocity_radians parameter is missing from the settings.")

        # Get the acceleration time
        if settings.Has("acceleration_time"):
            self.acceleration_time = settings["acceleration_time"].GetDouble()
            if self.acceleration_time < 0.0:
                raise Exception("The acceleration_time parameter must be non-negative.")
        else:
            raise Exception("The acceleration_time parameter is missing from the settings.")

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # Get Current time
        self.time = self.rotating_frame_model_part.ProcessInfo[KM.TIME]

        # Calculate the angle (theta) for the specified time
        self.__GetThetaAndOmega() 

        # Apply Rotation and Mesh Displacement to nodes.
        KratosCFD.RotatingFrameUtility.ApplyRotationAndMeshDisplacement(self.rotating_frame_model_part, self.axis_of_rotation, self.theta, self.center_of_rotation)

        # Apply Velocity to object.
        KratosCFD.RotatingFrameUtility.ApplyVelocityToRotatingObject(self.rotating_object_model_part, self.axis_of_rotation, self.omega, self.center_of_rotation)

    def __GetThetaAndOmega(self):
        # Calculate the angular acceleration (alpha)
        if np.isclose(self.acceleration_time, 0.0):
            alpha = 0.0
        else:
            alpha = self.target_angular_velocity_radians / self.acceleration_time
        
        # Calculate the angle (theta) based on current time
        if self.time <= self.acceleration_time:
            self.omega = alpha * self.time
            self.theta = 0.5 * alpha * self.time**2
        else:
            self.omega = self.target_angular_velocity_radians
            self.theta = 0.5 * alpha * self.acceleration_time**2 + self.omega * (self.time - self.acceleration_time)