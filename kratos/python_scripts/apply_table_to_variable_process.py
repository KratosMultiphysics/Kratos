import KratosMultiphysics as KM
import pandas as pd

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTableToVariableProcess(model, settings["Parameters"])


class ApplyTableToVariableProcess(KM.Process):
    r"""This class is used to apply a table to a scalar variable from entities in a model part.
    
    Example of table input types:
    |----------------------------------------------|
    |   "table_input_type"       : "csv",          |
    |   "table_input_parameters" : {               |
    |       "file_name"      : "path/to/file.csv"  |
    |       "delimiter"      : ",",                |
    |       "decimal"        : ".",                |
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
    | "\s+"  | multiple spaces                     |
    | ","    | comma                               |
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

        table_retriever = self._table_input_types[self.settings["table_input_type"].GetString()]
        self.table = table_retriever(self.settings["table_input_parameters"], self.model_part)


    def ExecuteInitializeSolutionStep(self):
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
        time = self.model_part.ProcessInfo[KM.TIME]
        if self.interval.IsInInterval(time):
            if self.settings["historical_variable"].GetBool():
                if self.settings["constrained"].GetBool():
                    KM.VariableUtils().ApplyFixity(self.variable, False, self.container)


    def Check(self):
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


    def _RetrieveTableFromCsv(settings, _):

        default_settings =  KM.Parameters("""{
            "file_name" : "",
            "delimiter" : ",",
            "decimal"   : ".",
            "skiprows"  : 0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        file_name = settings["file_name"].GetString()
        delimiter = settings["delimiter"].GetString()
        decimal = settings["decimal"].GetString()
        skiprows = settings["skiprows"].GetInt()
        data = pd.read_csv(file_name, delimiter=delimiter, decimal=decimal, skiprows=skiprows, header=None)

        table = KM.PiecewiseLinearTable()
        for row in data.values:
            table.AddRow(row[0], row[1])
        return table


    def _RetrieveTableFromModelPart(settings, model_part):
        """Helper class to retrieve a table from a model part."""

        default_settings = KM.Parameters("""{
            "table_id" : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        table_id = settings["table_id"].GetInt()
        return model_part.GetTable(table_id)


    _table_input_types = {
        "csv"        : _RetrieveTableFromCsv,
        "model_part" : _RetrieveTableFromModelPart,
    }
