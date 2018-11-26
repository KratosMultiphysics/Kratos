import KratosMultiphysics as KratosMultiphysics

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyPeriodicBoundaryConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyPeriodicBoundaryConditionProcess(KratosMultiphysics.Process):
    """This class is used in order to impose periodic boundary condition between two given modelparts which are related
        via a rotation operation or a translation operation.

    This class constructs the model parts containing the constraints that enforce the rigid body movement
    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                        : "This process uses LinearMasterSlaveConstraint in order to impose a periodic boundary condition on the given submodelparts. The process takes the first provided submodelpart as master and the second as slave.",
            "computing_model_part_name"   : "",
            "model_part_name"             : "please_specify_model_part_name",
            "first_model_part_name"       : "please_specify_model_part_name",
            "second_model_part_name"      : "please_specify_model_part_name",
            "interval"                    : [0.0, 1e30],
            "transformation_settings"     : {
                                                "rotation_settings":{
                                                    "center":[0,0,0],
                                                    "axis_of_rotation":[0.0,0.0,0.0],
                                                    "angle_degree":0.0
                                                },
                                                "translation_settings":{
                                                    "dir_of_translation":[0.0,0.0,0.0],
                                                    "magnitude":0.0
                                                }
                                            },
            "search_settings"             : {
                                                "max_results":100000,
                                                "tolerance": 1E-6
                                            },
            "variable_names"              : []
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

        if(settings["variable_names"].size() == 0):
            raise Exception("No variables specified to apply periodic boundary conditions.")

        # The computing model part
        main_model_part_name = settings["model_part_name"].GetString()
        self.main_model_part = Model[main_model_part_name]
        computing_model_part_name = settings["computing_model_part_name"].GetString()
        if(computing_model_part_name != ""):
            self.computing_model_part = Model[main_model_part_name+"."+computing_model_part_name]
        else:
            self.computing_model_part = Model[main_model_part_name]

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # Create the process
        periodic_parameters = KratosMultiphysics.Parameters()
        periodic_parameters.AddValue("variable_names", settings["variable_names"])
        periodic_parameters.AddValue("transformation_settings", settings["transformation_settings"])
        periodic_parameters.AddValue("search_settings", settings["search_settings"])


        master_model_part_name = main_model_part_name+"."+settings["first_model_part_name"].GetString()
        slave_model_part_name = main_model_part_name+"."+settings["second_model_part_name"].GetString()
        self.master_model_part = Model[master_model_part_name]
        self.slave_model_part = Model[slave_model_part_name]

        self.periodic_bc_process = KratosMultiphysics.ApplyPeriodicConditionProcess(self.master_model_part,
                                                                                    self.slave_model_part
                                                                                               , periodic_parameters)

    def ExecuteInitialize(self):

        """
        self.periodic_bc_process.ExecuteInitialize()
        list_constraints = [i.Id for i in self.master_model_part.MasterSlaveConstraints]
        self.computing_model_part.AddMasterSlaveConstraints(list_constraints)

    def ExecuteInitializeSolutionStep(self):


        """
        self.periodic_bc_process.ExecuteInitializeSolutionStep()

        current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

        # We activate/deactivate conditions dependeding of interval
        if (self.interval.IsInInterval(current_time)):
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, self.master_model_part.MasterSlaveConstraints)
        else:
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, False, self.master_model_part.MasterSlaveConstraints)
