# Importing the Kratos Library
import KratosMultiphysics

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignScalarVariableToEntitiesProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignScalarVariableToEntitiesProcess(KratosMultiphysics.Process):
    """This process assigns a given value (scalar) to the entities belonging a certain submodelpart

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
            "help"            : "This process assigns a given value (scalar) to the entities belonging a certain submodelpart",
            "mesh_id"         : 0,
            "model_part_name" : "please_specify_model_part_name",
            "variable_name"   : "SPECIFY_VARIABLE_NAME",
            "interval"        : [0.0, 1e30],
            "value"           : 0.0,
            "local_axes"      : {},
            "entities"        : []
        }
        """
        )

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # Here i do a trick, since i want to allow "value" to be a string or a double value
        if settings.Has("value"):
            if settings["value"].IsString():
                default_settings["value"].SetString("0.0")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.value_is_numeric = False

        # Define entities
        if settings["entities"].size() == 0:
           raise Exception("This process requires a list of entities. The options are: nodes, conditions, elements and constraints")
        self.entities = []
        for i in range(settings["entities"].size()):
            self.entities.append(settings["entities"][i].GetString())

        # Set parameters of the processes
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name", settings["model_part_name"])
        params.AddValue("mesh_id", settings["mesh_id"])
        params.AddValue("value", settings["value"])
        params.AddValue("variable_name", settings["variable_name"])

        # Set processes
        self.aux_processes = []
        if settings["value"].IsNumber():
            self.value_is_numeric = True
            for i in range(len(self.entities)):
                if self.entities[i] == "nodes":
                    self.aux_processes.append( KratosMultiphysics.AssignScalarVariableToNodesProcess(self.model_part, params))
                elif self.entities[i] == "conditions":
                    self.aux_processes.append( KratosMultiphysics.AssignScalarVariableToConditionsProcess(self.model_part, params))
                elif self.entities[i] == "elements":
                    self.aux_processes.append( KratosMultiphysics.AssignScalarVariableToElementsProcess(self.model_part, params))
                elif self.entities[i] == "constraints":
                    self.aux_processes.append( KratosMultiphysics.AssignScalarVariableToMasterSlaveConstraintsProcess(self.model_part, params))
        else:
            params.AddValue("local_axes", settings["local_axes"])
            for i in range(len(self.entities)):
                if self.entities[i] == "nodes":
                    self.aux_processes.append( KratosMultiphysics.AssignScalarFieldToNodesProcess(self.model_part, params))
                elif self.entities[i] == "conditions":
                    self.aux_processes.append( KratosMultiphysics.AssignScalarFieldToConditionsProcess(self.model_part, params))
                elif self.entities[i] == "elements":
                    self.aux_processes.append( KratosMultiphysics.AssignScalarFieldToElementsProcess(self.model_part, params))

        # construct a variable_utils object to speedup fixing
        self.step_is_active = False

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):
            self.step_is_active = True
            for process in self.aux_processes:
                process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.step_is_active = False
