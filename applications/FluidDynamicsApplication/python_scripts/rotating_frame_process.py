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
            "angular_velocity_radians": 0.0,
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
        if settings.Has("angular_velocity_radians"):
            self.angular_velocity_radians = settings["angular_velocity_radians"].GetDouble()
        else:
            raise Exception("The angular_velocity_radians parameter is missing from the settings.")

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
        self.get_theta_and_omega() #Potentially python

        # Calculate the rotation matrix.
        # transformation_matrix = self.get_rotation_matrix(self.theta) #Potentially C++

        # Apply Rotation and Mesh Displacement to nodes.
        # self.ApplyRotationAndMeshDisplacement(self.rotating_frame_model_part, transformation_matrix) #Potentially C++
        # KratosCFD.RotatingFrameUtility.ApplyRotationAndMeshDisplacement(self.rotating_frame_model_part, transformation_matrix)
        self.ApplyRotationAndMeshDisplacement(self.rotating_frame_model_part, self.axis_of_rotation, self.theta, self.center_of_rotation)

        # Apply Velocity to object.
        self.set_node_velocities(self.rotating_object_model_part, self.axis_of_rotation, self.omega, self.center_of_rotation)
        # self.ApplyVelocityToRotatingObject(self.rotating_object_model_part, self.omega, self.center_of_rotation, self.axis_of_rotation)
        # KratosCFD.RotatingFrameUtility.ApplyVelocityToRotatingObject(self.rotating_object_model_part, self.omega, self.axis_of_rotation)

    def ApplyVelocityToRotatingObject(self, rotating_object_model_part, omega, center_of_rotation, axis_of_rotation):
        for node in rotating_object_model_part.Nodes:
            relative_coordinates = [node.X, node.Y, node.Z]-center_of_rotation
            radius = np.linalg.norm(relative_coordinates)
            #Apply component by component the corresponding velocity.
            velocity_magnitude = omega * radius
            velocity_direction = np.cross(axis_of_rotation, relative_coordinates)
            velocity_vector = velocity_direction / np.linalg.norm(velocity_direction)
            velocity = velocity_magnitude * velocity_vector
            node.SetSolutionStepValue(KM.VELOCITY_X, velocity[0])
            node.SetSolutionStepValue(KM.VELOCITY_Y, velocity[1])
            node.SetSolutionStepValue(KM.VELOCITY_Z, velocity[2])
            node.Fix(KM.VELOCITY_X)
            node.Fix(KM.VELOCITY_Y)
            node.Fix(KM.VELOCITY_Z)
    
    def calculate_magnitude_difference(self, point, center_of_rotation, axis_of_rotation):
        # Step 1: Calculate the vector from the target point to the given point
        vector_difference = point - center_of_rotation

        # Step 2: Calculate the dot product between the vector difference and the target axis
        dot_product = np.dot(vector_difference, axis_of_rotation)

        # Step 3: Calculate the magnitude of the target axis
        magnitude_axis_of_rotation = np.linalg.norm(axis_of_rotation)

        # Step 4: Divide the dot product by the magnitude of the target axis to obtain the projection
        projection = dot_product / magnitude_axis_of_rotation

        # Step 5: Calculate the magnitude of the projection to obtain the desired magnitude of difference
        magnitude_difference = np.abs(projection)

        return magnitude_difference
    
    # def ApplyVelocityToRotatingObject(self, rotating_object_model_part, omega, center_of_rotation, axis_of_rotation):
    #     for node in rotating_object_model_part.Nodes:
    #         # Calculate the magnitude difference
    #         point = np.array([node.X, node.Y, node.Z])
    #         magnitude_difference = self.calculate_magnitude_difference(point, center_of_rotation, axis_of_rotation)

    #         # Apply component by component the corresponding velocity.
    #         velocity_magnitude = omega * magnitude_difference
    #         velocity_direction = np.cross(axis_of_rotation, point)
    #         velocity_vector = velocity_direction / np.linalg.norm(velocity_direction)
    #         velocity = velocity_magnitude * velocity_vector

    #         # Set the solution step values for velocity components
    #         node.SetSolutionStepValue(KM.VELOCITY_X, velocity[0])
    #         node.SetSolutionStepValue(KM.VELOCITY_Y, velocity[1])
    #         node.SetSolutionStepValue(KM.VELOCITY_Z, velocity[2])

    #         # Fix the velocity components to ensure they remain constant during the simulation
    #         node.Fix(KM.VELOCITY_X)
    #         node.Fix(KM.VELOCITY_Y)
    #         node.Fix(KM.VELOCITY_Z)

    def ApplyRotationAndMeshDisplacement(self, model_part, axis, theta, center):
        # Normalizing the rotation axis
        axis = axis / np.sqrt(np.dot(axis, axis))
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
    
    def set_node_velocities(self, model_part, axis, omega, center):
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

    # def ApplyRotationAndMeshDisplacement(self, rotating_frame_model_part, transformation_matrix):
    #     for node in rotating_frame_model_part.Nodes:
    #         coor = np.array([node.X0, node.Y0, node.Z0, 1.0])  # Convert point to homogeneous coordinates
    #         rotated_coor = transformation_matrix @ coor
    #         rotated_coor = rotated_coor[:3] / rotated_coor[3]  # Convert back to 3D coordinates
    #         node.SetSolutionStepValue(KM.MESH_DISPLACEMENT_X, rotated_coor[0] - node.X0)
    #         node.SetSolutionStepValue(KM.MESH_DISPLACEMENT_Y, rotated_coor[1] - node.Y0)
    #         node.SetSolutionStepValue(KM.MESH_DISPLACEMENT_Z, rotated_coor[2] - node.Z0)
    #         node.Fix(KM.MESH_DISPLACEMENT_X)
    #         node.Fix(KM.MESH_DISPLACEMENT_Y)
    #         node.Fix(KM.MESH_DISPLACEMENT_Z)
    #         node.X = rotated_coor[0]
    #         node.Y = rotated_coor[1]
    #         node.Z = rotated_coor[2]
    
    def get_theta_and_omega(self):
        # Calculate the angular acceleration (alpha)
        if np.isclose(self.acceleration_time, 0.0):
            alpha = 0.0
        else:
            alpha = self.angular_velocity_radians / self.acceleration_time
        
        # Calculate the angle (theta) based on current time
        if self.time <= self.acceleration_time:
            self.omega = alpha * self.time
            self.theta = 0.5 * alpha * self.time**2
        else:
            self.omega = self.angular_velocity_radians
            self.theta = 0.5 * alpha * self.acceleration_time**2 + self.omega * (self.time - self.acceleration_time)

    def get_rotation_matrix(self, theta):
        # Create a rotation matrix for the specified axis of rotation
        axis = self.axis_of_rotation
        c = np.cos(theta)
        s = np.sin(theta)
        C = 1 - c
        x, y, z = axis
        rotation_matrix = np.array([[x*x*C + c, x*y*C - z*s, x*z*C + y*s],
                                    [y*x*C + z*s, y*y*C + c, y*z*C - x*s],
                                    [z*x*C - y*s, z*y*C + x*s, z*z*C + c]])

        # Create translation matrices for moving the center of rotation
        translation_to_origin = np.eye(4)
        translation_to_origin[:3, 3] = -self.center_of_rotation
        translation_back = np.eye(4)
        translation_back[:3, 3] = self.center_of_rotation

        # Create a homogeneous rotation matrix
        homogeneous_rotation = np.eye(4)
        homogeneous_rotation[:3, :3] = rotation_matrix

        # Combine the translation and rotation matrices to get the final transformation matrix
        transformation_matrix = np.dot(np.dot(translation_back, homogeneous_rotation), translation_to_origin)
        return transformation_matrix

    def ExecuteBeforeOutputStep(self):
        """ This method is executed before writing the output (if output
        is being written in this step)

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteAfterOutputStep(self):
        """ This method is executed after writing the output (if output
        is being written in this step)

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalize(self):
        """ This method is executed after the computations, at the end of the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass