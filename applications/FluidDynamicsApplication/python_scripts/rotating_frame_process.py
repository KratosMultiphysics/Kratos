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

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def Check(self):
        """ This method verifies that the input is correct

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed just before the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

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
        self.ApplyVelocityToRotatingObject(self.rotating_object_model_part, self.axis_of_rotation, self.omega, self.center_of_rotation)
        # KratosCFD.RotatingFrameUtility.ApplyVelocityToRotatingObject(self.rotating_object_model_part, self.omega, self.axis_of_rotation)

    def ApplyRotationAndMeshDisplacement(self, model_part, axis, theta, center):
        # Normalizing the rotation axis
        a = np.cos(theta / 2)
        b, c, d = -axis * np.sin(theta / 2)

        # Creating a quaternion rotation matrix
        rot_matrix = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                            [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                            [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

        for node in model_part.Nodes:
            # Getting the initial coordinates of the node
            point = np.array([node.X0, node.Y0, node.Z0])

            # Shifting the rotation center to the origin
            point -= center

            # Applying the rotation
            rotated_point = np.dot(point, rot_matrix)

            # Shifting the point back and updating the node coordinates
            rotated_point += center
            node.X, node.Y, node.Z = rotated_point
            node.SetSolutionStepValue(KM.MESH_DISPLACEMENT_X, rotated_point[0] - node.X0)
            node.SetSolutionStepValue(KM.MESH_DISPLACEMENT_Y, rotated_point[1] - node.Y0)
            node.SetSolutionStepValue(KM.MESH_DISPLACEMENT_Z, rotated_point[2] - node.Z0)
            node.Fix(KM.MESH_DISPLACEMENT_X)
            node.Fix(KM.MESH_DISPLACEMENT_Y)
            node.Fix(KM.MESH_DISPLACEMENT_Z)
    
    def ApplyVelocityToRotatingObject(self, model_part, axis, omega, center):
        # Normalizing the rotation axis and creating the angular velocity vector
        axis = axis / np.sqrt(np.dot(axis, axis))
        angular_velocity_vector = omega * axis

        for node in model_part.Nodes:
            # Getting the current coordinates of the node
            point = np.array([node.X, node.Y, node.Z])

            # Calculating the position vector (relative to the rotation center)
            position_vector = point - center

            # Computing the velocity due to rotation (v = omega cross r)
            velocity_vector = np.cross(angular_velocity_vector, position_vector)

            # Setting the node's velocity
            node.SetSolutionStepValue(KM.VELOCITY_X, velocity_vector[0])
            node.SetSolutionStepValue(KM.VELOCITY_Y, velocity_vector[1])
            node.SetSolutionStepValue(KM.VELOCITY_Z, velocity_vector[2])

            # Fix the velocity components to ensure they remain constant during the simulation
            node.Fix(KM.VELOCITY_X)
            node.Fix(KM.VELOCITY_Y)
            node.Fix(KM.VELOCITY_Z)

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