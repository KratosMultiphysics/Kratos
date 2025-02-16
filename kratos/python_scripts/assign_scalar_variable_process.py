# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.read_csv_table_utility import ReadCsvTableUtility


def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception(f"Expecting a Parameters object, but got {settings}.")
    return AssignScalarVariableProcess(Model, settings["Parameters"])

class AssignScalarVariableProcess(KratosMultiphysics.Process):
    """This process sets a given scalar value for a certain variable in all the nodes of a submodelpart

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.

    Possible specifications for the key 'value' from settings:
    - double: the constant value is applied. Example: 1.0
    - string: the string is parsed as a function. Example: "(x*sin(y))*exp(-t^2)"
    - parameters: a csv table can be specified. Example:
        {
            "name"       : "csv_table",
            "filename"   : "path/to/file.csv",
            "delimiter"  : ",",
            "skiprows"   : 0,
        }
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)
        settings.ValidateAndAssignDefaults(self.GetSchema())

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        if not isinstance(self.variable, KratosMultiphysics.DoubleVariable) and not isinstance(self.variable, KratosMultiphysics.VectorVariable):
            msg = "Error in AssignScalarToNodesProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a scalar or a component"
            raise Exception(msg)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.mesh = self.model_part.GetMesh(settings["mesh_id"].GetInt())
        self.is_fixed = settings["constrained"].GetBool()

        self.value_is_numeric = False
        self.value_is_function = False
        if settings["value"].IsNumber():
            self.value_is_numeric = True
            self.value = settings["value"].GetDouble()
        elif settings["value"].IsString():
            self.value_is_function = True
            self.function_string = settings["value"].GetString()
            self.aux_function = KratosMultiphysics.GenericFunctionUtility(self.function_string, settings["local_axes"])

            if self.aux_function.DependsOnSpace():
                self.cpp_apply_function_utility = KratosMultiphysics.ApplyFunctionToNodesUtility(self.mesh.Nodes, self.aux_function )
        else:
            self.table = ReadCsvTableUtility(settings["value"]).Read(self.model_part)

        # Construct a variable_utils object to speedup fixing
        self.variable_utils = KratosMultiphysics.VariableUtils()
        self.step_is_active = False

    def ExecuteBeforeSolutionLoop(self):
        """This method is executed in before initialize the solution step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        """This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):

            self.step_is_active = True

            if self.is_fixed:
                self.variable_utils.ApplyFixity(self.variable, self.is_fixed, self.mesh.Nodes)

            if self.value_is_numeric:
                self.variable_utils.SetVariable(self.variable, self.value, self.mesh.Nodes)
            elif self.value_is_function:
                if self.aux_function.DependsOnSpace() == False: #depends on time only
                    self.value = self.aux_function.CallFunction(0.0,0.0,0.0,current_time,0.0,0.0,0.0)
                    self.variable_utils.SetVariable(self.variable, self.value, self.mesh.Nodes)
                else: #most general case - space varying function (possibly also time varying)
                    self.cpp_apply_function_utility.ApplyFunction(self.variable, current_time)
            else:
                self.value = self.table.GetValue(current_time)
                self.variable_utils.SetVariable(self.variable, self.value, self.mesh.Nodes)

    def ExecuteFinalizeSolutionStep(self):
        """This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if self.step_is_active:
            # Here we free all of the nodes in the mesh
            if self.is_fixed:
                fixity_status  = False
                self.variable_utils.ApplyFixity(self.variable, fixity_status, self.mesh.Nodes)

        self.step_is_active = False

    @classmethod
    def GetSchema(cls) -> KratosMultiphysics.Schema:
        definition: KratosMultiphysics.Parameters = KratosMultiphysics.Parameters("""{
            "title" : "AssignScalarVariableToEntitiesProcess",
            "description" : "A process for setting scalar kratos variables on all nodes of the provided model part.",
            "type" : "object",
            "properties" : {
                "model_part_name" : {
                    "description" : "Name of the model part to assign variables in.",
                    "type" : "string",
                    "pattern" : ".+"
                },
                "mesh_id" : {
                    "type" : "integer",
                    "default" : 0
                },
                "variable_name" : {
                    "description" : "Name of the variable to set.",
                    "type" : "string",
                    "pattern" : "[A-Z0-9_]+"
                },
                "interval" : {
                    "description" : "Populated at runtime."
                },
                "constrained" : {
                    "description" : "Set a Dirichlet condition on the provided variable of the specified model part's nodes if this is true.",
                    "type" : "boolean",
                    "default" : true
                },
                "value" : {
                    "description" : "Value to set the specified variable to. The value may be provided as a literal number, a string defining a function of x, y, z, t, X, Y, Z, or an object specifying a CSV file.",
                    "anyOf": [
                        {"type" : "bool"},
                        {"type" : "number"},
                        {"type" : "string"},
                        {"type" : "object", "description" : "Populated at runtime."}
                    ],
                    "default" : {}
                },
                "local_axes" : {
                    "decription" : "Local coordinate system to compute the value in if it is a string.",
                    "type" : "object",
                    "properties" : {
                        "origin" : {
                            "description" : "Origin of the local coordinate system.",
                            "type" : "array",
                            "minItems" : 3,
                            "maxItems" : 3,
                            "items" : "number",
                            "default" : [0,0,0]
                        },
                        "axes" : {
                            "description" : "3D rotation matrix defining the local basis vectors from the global unit vectors.",
                            "type" : "array",
                            "items" : {
                                "type" : "array",
                                "minItems" : 3,
                                "maxItems" : 3
                            },
                            "default" : [[1,0,0],[0,1,0],[0,0,1]]
                        }
                    },
                    "default" : {}
                }
            },
            "additionalProperties" : false,
            "required" : ["model_part_name", "variable_name"]
        }""")

        # Define the schema when the value refers to a CSV file.
        definition["properties"]["value"]["anyOf"][3] = ReadCsvTableUtility.GetSchema().GetDefinition()

        # Define the interval's schema.
        definition["properties"]["interval"] = KratosMultiphysics.IntervalUtility(KratosMultiphysics.Parameters()).GetSchema().GetDefinition()["properties"]["interval"]

        return KratosMultiphysics.Schema(definition)
