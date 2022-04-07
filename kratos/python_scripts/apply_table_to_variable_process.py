import KratosMultiphysics as KM
import csv

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTableToVariableProcess(model, settings["Parameters"])


class ApplyTableToVariableProcess(KM.Process):
    r"""This class is used to apply a table to a scalar variable from entities in a model part.

    Example of table input type:
    |----------------------------------------------|
    |   "table_input_type"       : "csv",          |
    |   "table_input_parameters" : {               |
    |       "file_name"      : "path/to/file.csv"  |
    |       "delimiter"      : ",",                |
    |       "skiprows"       : 0                   |
    |   }                                          |
    |----------------------------------------------|
    |   "table_input_type"       : "model_part",   |
    |   "table_input_parameters" : {               |
    |       "table_id"       : 1                   |
    |   }                                          |
    |----------------------------------------------|

    Some possible delimiters for csv:
    |--------|-------------------------------------|
    | "\t"   | tab                                 |
    | " "    | spaces                              |
    | ","    | comma (default)                     |
    |        | etc.                                |
    |--------|-------------------------------------|
    """


    @staticmethod
    def GetDefaultParameters():
        default_parameters = KM.Parameters("""{
            "model_part_name"        : "",
            "container"              : "Nodes",
            "historical_variable"    : true,
            "constrained"            : false,
            "variable_name"          : "",
            "interval"               : [0, "End"],
            "table_input_type"       : "csv",
            "table_input_parameters" : {}
        }""")
        return default_parameters


    def __init__(self, model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the model contaning the model_parts
        settings -- Kratos parameters containing solver settings.
        """

        KM.Process.__init__(self)
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.settings = settings
        self.model_part = model.GetModelPart(self.settings["model_part_name"].GetString())
        self.variable = KM.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())
        self.interval = KM.IntervalUtility(self.settings)
        self.container = getattr(self.model_part, self.settings["container"].GetString())

        table_retriever_name = self._table_input_types[self.settings["table_input_type"].GetString()]
        table_retriever = getattr(self, table_retriever_name)
        self.table = table_retriever(self.settings["table_input_parameters"])


    def ExecuteInitializeSolutionStep(self):
        """Apply the value and set fixity if defined."""
        time = self.model_part.ProcessInfo[KM.TIME]
        if self.interval.IsInInterval(time):
            value = self.table.GetValue(time)
            if self.settings["historical_variable"].GetBool():
                KM.VariableUtils().SetVariable(self.variable, value, self.container)
                if self.settings["constrained"].GetBool():
                    KM.VariableUtils().ApplyFixity(self.variable, True, self.container)
            else:
                KM.VariableUtils().SetNonHistoricalVariable(self.variable, value, self.container)


    def ExecuteFinalizeSolutionStep(self):
        """Unset fixity if defined."""
        time = self.model_part.ProcessInfo[KM.TIME]
        if self.interval.IsInInterval(time):
            if self.settings["historical_variable"].GetBool():
                if self.settings["constrained"].GetBool():
                    KM.VariableUtils().ApplyFixity(self.variable, False, self.container)


    def Check(self):
        """Check the correctness of the input."""
        if not self.settings["model_part_name"].GetString():
            raise Exception(self.__class__.__name__, "The model part name must be a valid string.")

        if not self.settings["variable_name"].GetString():
            raise Exception(self.__class__.__name__, "The variable name must be a valid string.")

        if not self.settings["container"].GetString() in ["Nodes", "Elements", "Conditions"]:
            raise Exception(self.__class__.__name__, "Wrong container type.")

        if self.settings["historical_variable"].GetBool():
            if self.settings["container"].GetString() != "Nodes":
                raise Exception(self.__class__.__name__, "The historical database only is valid for nodes.")

        if not isinstance(self.variable, KM.DoubleVariable):
            raise Exception(self.__class__.__name__, "Only scalar variables are supported.")


    def _RetrieveTableFromCsv(self, settings):
        """Helper function to retrieve a table from a csv file."""
        default_settings =  KM.Parameters("""{
            "file_name" : "",
            "delimiter" : ",",
            "skiprows"  : 0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        file_name = settings["file_name"].GetString()
        delimiter = settings["delimiter"].GetString()
        skiprows = settings["skiprows"].GetInt()

        table = KM.PiecewiseLinearTable()
        with open(file_name, 'r') as table_file:
            data = csv.reader(table_file, delimiter=delimiter, skipinitialspace=True)
            for _ in range(skiprows):
                next(data,None)
            for row in data:
                if row:  # skip empty rows
                    table.AddRow(float(row[0]), float(row[1]))
        return table


    def _RetrieveTableFromModelPart(self, settings):
        """Helper function to retrieve a table from a model part."""
        default_settings = KM.Parameters("""{
            "table_id" : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        table_id = settings["table_id"].GetInt()
        return self.model_part.GetTable(table_id)


    _table_input_types = {
        "csv"        : "_RetrieveTableFromCsv",
        "model_part" : "_RetrieveTableFromModelPart",
    }
