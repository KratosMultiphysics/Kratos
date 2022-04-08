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
            "file_name"  : "",
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

    def Read(self, model_part):
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
                        raise Exception(f"Only 2-column tables are supported. However, a {len(row)}-column row is found.")
                    elif len(row) > 2:
                        KM.Logger.PrintWarning(f"Only 2-column tables are supported. However, a {len(row)}-column row is found. Extra columns will be ignored.")
                    table.AddRow(self._Float(row[0]), self._Float(row[1]))
        if self.table_id > -1:
            model_part.AddTable(table)
        return table

    def _Float(self, value):
        try:
            return float(value)
        except ValueError:
            return self.na_replace

