import KratosMultiphysics as KM
import csv

class ReadCsvTableUtility:
    r"""This class is used to retrieve a table from the specified parameters.

    The following parameters can be specified:
    |---------------|---------------------------------------------------|
    | "name"        | The type of input (csv_table)                     |
    |---------------|---------------------------------------------------|
    | "file_name"   | The file name                                     |
    |---------------|---------------------------------------------------|
    | "delimiter"   | ","  comma                                        |
    |               | ";"  semicolon                                    |
    |               | "\t" tab                                          |
    |               | " "  spaces                                       |
    |               |      etc                                          |
    |---------------|---------------------------------------------------|
    | "skiprows"    | The number of rows to skip before reading data    |
    |---------------|---------------------------------------------------|
    | "table_id"    | If >-1 the input table will be stored in the      |
    |               | model part                                        |
    |---------------|---------------------------------------------------|
    | "na_replace"  | The value to apply when N/A is read               |
    |---------------|---------------------------------------------------|
    """

    def __init__(self, settings):
        """Constructor of the csv table reader: validate the parameters.

        Keyword arguments:
        self -- It signifies an instance of the class.
        settings -- Kratos parameters containing solver settings.
        """
        default_settings =  KM.Parameters("""{
            "name"       : "csv_table",
            "filename"  : "",
            "delimiter"  : ",",
            "skiprows"   : 0,
            "table_id"   : -1,
            "na_replace" : 0.0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.file_name = settings["file_name"].GetString()
        self.delimiter = settings["delimiter"].GetString()
        self.skiprows = settings["skiprows"].GetInt()
        self.table_id = settings["table_id"].GetInt()
        self.na_replace = settings["na_replace"].GetDouble()

    def Read(self, model_part = None):
        """Read a csv table.

        Keyword arguments:
        self -- It signifies an instance of the class.
        model_part -- ModelPart where to store or apply the table.
        """
        table = KM.PiecewiseLinearTable()
        with open(self.file_name, 'r') as table_file:
            data = csv.reader(table_file, delimiter=self.delimiter, skipinitialspace=True)
            for _ in range(self.skiprows):
                next(data)
            for row in data:
                if row:  # skip empty rows
                    if len(row) < 2:
                        raise Exception("Only 2-column tables are supported. However, a {}-column row is found.".format(len(row)))
                    elif len(row) > 2:
                        KM.Logger.PrintWarning("Only 2-column tables are supported. However, a {}-column row is found. Extra columns will be ignored.".format(len(row)))
                    table.AddRow(self._Float(row[0]), self._Float(row[1]))
        if self.table_id > -1:
            if model_part:
                model_part.AddTable(self.table_id, table)
            else:
                err_msg = "Asking to save table with id {} but no model part is provided.".format(self.table_id)
                raise Exception(err_msg)
        return table

    def _Float(self, value):
        try:
            return float(value)
        except ValueError:
            return self.na_replace

