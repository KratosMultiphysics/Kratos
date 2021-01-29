# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetInitialStateProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "Process"

class SetInitialStateProcess(KratosMultiphysics.Process):
    """This process sets a given value for a certain flag in all the nodes of a submodelpart

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

        #The value can be a double or a string (function)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "help"            : "This sets the initial conditions in terms of imposed strain, stress or deformation gradient",
            "mesh_id"         : 0,
            "model_part_name" : "please_specify_model_part_name",
            "dimension"       : 3,
            "imposed_strain"  : [0,0,0,0,0,0],
            "imposed_stress"  : [0,0,0,0,0,0],
            "imposed_deformation_gradient"  : [[1,0,0],[0,1,0],[0,0,1]],
            "interval"        : [0.0, 1e30]
        }
        """
        )

        #assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.dimension = settings["dimension"].GetInt()

        self.imposed_strain = settings["imposed_strain"]
        self.imposed_stress = settings["imposed_stress"]
        self.imposed_deformation_gradient = settings["imposed_deformation_gradient"]


    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

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
        if self.interval.IsInInterval(current_time):
            self.SetInitialState()


    def SetInitialState(self):
        """ This method creates the c++ process and sets each entity

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if self.dimension == 3:
            KratosMultiphysics.SetInitialStateProcess3D(self.model_part,
                self.imposed_strain, self.imposed_stress, self.imposed_deformation_gradient).ExecuteInitializeSolutionStep()
        else:  # 2D case
            KratosMultiphysics.SetInitialStateProcess2D(self.model_part,
                self.imposed_strain, self.imposed_stress, self.imposed_deformation_gradient).ExecuteInitializeSolutionStep()
