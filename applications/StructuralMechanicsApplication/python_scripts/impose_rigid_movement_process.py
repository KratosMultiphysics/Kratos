import KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeRigidMovementProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ImposeRigidMovementProcess(KratosMultiphysics.Process):
    """This class is used in order to impose a rigid body movement in a certain region of the problem

    This class constructs the model parts containing the constrains that enforce the rigid body movement
    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"      : "please_specify_model_part_name",
            "new_model_part_name"  : "Rigid_Movement_ModelPart",
            "interval"             : [0.0, 1e30],
            "penalty_parameter"    : 1.0e8,
            "penalty_parameter"    : "DISPLACEMENT",
            "reference_node_id"    : 0
        }
        """)

        # Detect "End" as a tag and replace it by a large number
        if(settings.Has("interval")):
            if(settings["interval"][1].IsString()):
                if(settings["interval"][1].GetString() == "End"):
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # The penalty parameter
        self.penalty_parameter = settings["penalty_parameter"].GetDouble()

        # The reference node id
        self.reference_node_id = settings["reference_node_id"].GetInt()

        # We get the corresponding model parts
        self.model_part = Model[settings["model_part_name"].GetString()]
        if (self.model_part.HasSubModelPart(settings["new_model_part_name"].GetString())):
            self.rigid_model_part = self.model_part.GetSubModelPart(settings["new_model_part_name"].GetString())
        else:
            self.rigid_model_part = self.model_part.CreateSubModelPart(settings["new_model_part_name"].GetString())

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        # We activate/deactivate conditions dependeding of interval
        if (self.interval.IsInInterval(current_time)):
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, self.rigid_model_part.MasterSlaveConstraints)
        else:
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, False, self.rigid_model_part.MasterSlaveConstraints)
