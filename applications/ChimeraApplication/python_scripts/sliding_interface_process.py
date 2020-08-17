import KratosMultiphysics
import KratosMultiphysics.ChimeraApplication as KratosChimera

def Factory(settings, Model):
    if ( not isinstance(settings, KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplySlidingInterfaceProcess(Model, settings["Parameters"])

class ApplySlidingInterfaceProcess(KratosMultiphysics.Process):
    """This process applies a sliding interface on a given set of modelparts

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
            "help"                        : "This process uses LinearMasterSlaveConstraint in order to impose a sliding interface condition on the given submodelparts. The process takes the first provided submodelpart as master and the second as slave.",
            "first_model_part_name"       : "master_non_moving_mp",
            "second_model_part_name"      : "slave_moving_mp",
            "model_part_name"             : "FluidModelPart",
            "computing_model_part_name"   : "fluid_computational_model_part",
            "interval"                    : [0.0, 1e30],
            "search_settings"             : {},
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
        master_model_part_name = settings["first_model_part_name"].GetString()
        slave_model_part_name = settings["second_model_part_name"].GetString()
        self.master_model_part = Model[master_model_part_name]
        self.slave_model_part = Model[slave_model_part_name]

        # Create the process
        sliding_parameters = KratosMultiphysics.Parameters()
        sliding_parameters.AddValue("variable_names", settings["variable_names"])
        sliding_parameters.AddValue("search_settings", settings["search_settings"])
        nD = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        ## Here we create the c++ process with all the necessary parameters
        if(nD == 2):
            self.sliding_interface_proc = KratosChimera.SlidingInterfaceProcess2D(self.master_model_part , self.slave_model_part , sliding_parameters)
        elif(nD == 3):
            self.sliding_interface_proc = KratosChimera.SlidingInterfaceProcess3D(self.master_model_part , self.slave_model_part , sliding_parameters)
        else:
            raise Exception("SlidingInterfaceProcess is only for 2 and 3 dimension problems !! Check input.")

    def ExecuteInitializeSolutionStep(self):
        """
            Here the linear master slave constraints are made constraints are made/remade depending on the settings.
        """
        current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        self.sliding_interface_proc.ExecuteInitializeSolutionStep()
        if (self.interval.IsInInterval(current_time)):
            list_constraints = [i.Id for i in self.master_model_part.MasterSlaveConstraints]
            self.computing_model_part.AddMasterSlaveConstraints(list_constraints)


    def ExecuteFinalizeSolutionStep(self):
        current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        if (self.interval.IsInInterval(current_time)):
            self.sliding_interface_proc.ExecuteFinalizeSolutionStep()

    def ExecuteInitialize(self):
        self.sliding_interface_proc.ExecuteInitialize()

    def ExecuteFinalize(self):
        self.sliding_interface_proc.ExecuteFinalize()
