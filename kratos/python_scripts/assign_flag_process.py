from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignFlagProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "Process"

class AssignFlagProcess(KratosMultiphysics.Process):
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
            "help"            : "This process sets a given value for a certain flag in all the given entities of a submodelpart",
            "mesh_id"         : 0,
            "model_part_name" : "please_specify_model_part_name",
            "flag_name"       : "SPECIFY_FLAG_NAME",
            "interval"        : [0.0, 1e30],
            "value"           : true,
            "entities"        : []
        }
        """
        )

        #assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        settings.ValidateAndAssignDefaults(default_settings)

        self.flag = KratosMultiphysics.KratosGlobals.GetFlag(settings["flag_name"].GetString())

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.value = settings["value"].GetBool()
        # Define entities
        if (settings["entities"].size() == 0):
           raise Exception("This process requires a list of entities. The options are: nodes, conditions, elements and model_part")
        self.entities = []
        for i in range(settings["entities"].size()):
            self.entities.append(settings["entities"][i].GetString())

        # Construct a variable_utils object to speedup fixing
        self.flag_utils = KratosMultiphysics.VariableUtils()
        self.step_is_active = False

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

        if(self.interval.IsInInterval(current_time)):
            for i in range(len(self.entities)):
                if (self.entities[i] == "nodes"):
                    self.flag_utils.SetFlag(self.flag, self.value, self.model_part.Nodes)
                elif (self.entities[i] == "conditions"):
                    self.flag_utils.SetFlag(self.flag, self.value, self.model_part.Conditions)
                elif (self.entities[i] == "elements"):
                    self.flag_utils.SetFlag(self.flag, self.value, self.model_part.Elements)
                elif (self.entities[i] == "model_part"):
                    self.model_part.Set(self.flag, self.value)
