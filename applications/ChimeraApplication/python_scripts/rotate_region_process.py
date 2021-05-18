import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera

def Factory(settings, Model):
    if ( not isinstance(settings, KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyRotateRegionProcess(Model, settings["Parameters"])


class ApplyRotateRegionProcess(KratosMultiphysics.Process):
    """This process applies a rotation to a given modelpart or a submodelpart

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """
    def __init__(self, Model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)
        # settings for inlet with interface between fluids and separate velocities
        default_settings = KratosMultiphysics.Parameters("""
        {
                "model_part_name":"",
                "center_of_rotation":[0.0,0.0,0.0],
                "calculate_torque":false,
                "torque_model_part_name":"",
                "moment_of_inertia":0.0,
                "rotational_damping":0.0,
                "angular_velocity_radians":0.0,
                "axis_of_rotation":[0.0,0.0,0.0],
                "is_ale" : false,
                "interval": [0.0, 1e30]
        }
        """)

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # compare against the appropriate default settings
        settings.ValidateAndAssignDefaults(default_settings)
        # checking for empty model part name
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("ApplyRotateRegionProcess: A value (string) for the entry 'model_part_name' must be given in the parameters of the process.")
        # Get the modelpart to rotate
        self.model_part = Model[settings["model_part_name"].GetString()]

        if ( settings["axis_of_rotation"].IsVector() ):
            axis_of_rotation = settings["axis_of_rotation"].GetVector()
            if ( axis_of_rotation[0] == 0.0 and axis_of_rotation[1] == 0.0 and axis_of_rotation[2] == 0.0):
                raise Exception("The values (vector) of the entry 'axis_of_rotation' are all zero. This is not admissible.")

        if (settings["calculate_torque"].GetBool() and settings["angular_velocity_radians"].GetDouble()!=0.0):
            raise Exception("'calculate_torque' is set to true and 'angular_velocity_radians' is not zero. This is not admissible.")

        if(settings["calculate_torque"].GetBool() and settings["moment_of_inertia"].GetDouble() == 0.0):
            KratosMultiphysics.Logger.PrintWarning("RotateRegionProcess", " 'moment_of_inertia' is zero !!")
        if(settings["calculate_torque"].GetBool() and settings["rotational_damping"].GetDouble() == 0.0):
            KratosMultiphysics.Logger.PrintWarning("RotateRegionProcess", " 'rotational_damping' is zero !!")

        # If no torque_model_part_name is specified remove it to avoid later problems
        if (settings["torque_model_part_name"].GetString() == ""):
            settings.RemoveValue("torque_model_part_name")
        settings.RemoveValue("interval")

        # Making the actual process
        self.rotate_region_process = KratosChimera.RotateRegionProcess(self.model_part, settings)

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):
            self.rotate_region_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):
            self.rotate_region_process.ExecuteFinalizeSolutionStep()
