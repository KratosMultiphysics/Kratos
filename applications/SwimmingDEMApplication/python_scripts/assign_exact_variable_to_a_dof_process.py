from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as Kratos

def Factory(settings, Model):
    if not isinstance(settings, Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignExactVariableToADOFProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignExactVariableToADOFProcess(Kratos.Process):
    def __init__(self, Model, settings):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        Kratos.Process.__init__(self)

        default_settings = Kratos.Parameters("""
        {
            "help"                 : "This process sets a variable a certain scalar value in a given direction, for all the nodes belonging to a submodelpart. Uses assign_scalar_variable_to_conditions_process for each component",
            "mesh_id"              : 0,
            "model_part_name"      : "please_specify_model_part_name",
            "variable_name"        : "SPECIFY_VARIABLE_NAME",
            "exact_variable_name"  : "SPECIFY_EXACT_VARIABLE_NAME",
            "interval"             : [0.0, 1e30]
        }
        """)

        # Detect "End" as a tag and replace it by a large number
        if settings.Has("interval"):
            if settings["interval"][1].IsString():
                if settings["interval"][1].GetString() == "End":
                    settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
                else:
                    raise Exception("The second value of interval can be \"End\" or a number, interval currently:" + settings["interval"].PrettyPrintJsonString())

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.aux_processes = []

        self.destination_variable_name = settings["variable_name"].GetString()
        self.destination_variable = Kratos.KratosGlobals.GetVariable(self.destination_variable_name)

        self.exact_variable_name = settings["exact_variable_name"].GetString()
        self.exact_variable = Kratos.KratosGlobals.GetVariable(self.exact_variable_name)

        self.variable_type = Kratos.KratosGlobals.GetVariableType(self.destination_variable_name)
        if self.variable_type == "Array":
            self.destination_variable_name_component = []
            self.exact_variable_name_component = []
            for string in ["_X", "_Y", "_Z"]:
                self.destination_variable_name_component.append(Kratos.KratosGlobals.GetVariable((self.destination_variable_name + string)))
                self.exact_variable_name_component.append(Kratos.KratosGlobals.GetVariable((self.exact_variable_name + string)))

    def ExecuteBeforeSolutionLoop(self):
        """
        This method is executed in before initialize the solution loop
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.ExecuteInitializeSolutionStep()


    def ExecuteInitializeSolutionStep(self):
        """
        This method is executed in order to initialize the current step
        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        if self.variable_type == "Array":
            Kratos.VariableUtils().CopyVectorVar(self.exact_variable, self.destination_variable, self.model_part.Nodes)
            for node in self.model_part.Nodes:
                node.Fix(self.destination_variable_name_component[0])
                node.Fix(self.destination_variable_name_component[1])
                node.Fix(self.destination_variable_name_component[2])
        else:
            Kratos.VariableUtils().CopyScalarVar(self.exact_variable, self.destination_variable, self.model_part.Nodes)
            for node in self.model_part.Nodes:
                node.Fix(self.destination_variable)

    def ExecuteFinalizeSolutionStep(self):
        """
        This method is executed in order to finalize the current step
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        for process in self.aux_processes:
            process.ExecuteFinalizeSolutionStep()