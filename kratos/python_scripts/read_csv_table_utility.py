import KratosMultiphysics as KM
import csv

class ReadCsvTableUtility:
    r"""This class is used to retrieve a table from the specified parameters.

    The following parameters can be specified:
    |--------------------|------------------------------------------------------|
    | "name"             | The type of input (csv_table)                        |
    |--------------------|------------------------------------------------------|
    | "filename"         | The file name                                        |
    |--------------------|------------------------------------------------------|
    | "delimiter"        | ","  comma                                           |
    |                    | ";"  semicolon                                       |
    |                    | "\t" tab                                             |
    |                    | " "  spaces                                          |
    |                    |      etc                                             |
    |---------------- ---|------------------------------------------------------|
    | "skiprows"         | The number of rows to skip before reading data       |
    |--------------- ----|------------------------------------------------------|
    | "first_column_id"  | The index of the first column to read (zero-based)   |
    |----------------- --|------------------------------------------------------|
    | "second_column_id" | The index of the second column to read (zero-based)  |
    |--------------------|------------------------------------------------------|
    | "table_id"         | If >-1 the input table will be stored in the         |
    |                    | model part                                           |
    |--------------------|------------------------------------------------------|
    | "na_replace"       | The value to apply when N/A is read                  |
    |--------------------|------------------------------------------------------|
    """

    def __init__(self, settings):
        """Constructor of the csv table reader: validate the parameters.

        Keyword arguments:
        self -- It signifies an instance of the class.
        settings -- Kratos parameters containing solver settings.
        """
        default_settings =  KM.Parameters("""{
            "name"             : "csv_table",
            "filename"         : "",
            "delimiter"        : ",",
            "skiprows"         : 0,
            "first_column_id"  : 0,
            "second_column_id" : 1,
            "table_id"         : -1,
            "na_replace"       : 0.0
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.filename = settings["filename"].GetString()
        self.delimiter = settings["delimiter"].GetString()
        self.skiprows = settings["skiprows"].GetInt()
        self.first_column_id = settings["first_column_id"].GetInt()
        self.second_column_id = settings["second_column_id"].GetInt()
        self.table_id = settings["table_id"].GetInt()
        self.na_replace = settings["na_replace"].GetDouble()

    def Read(self, model_part = None):
        """Read a csv table.

        Keyword arguments:
        self -- It signifies an instance of the class.
        model_part -- ModelPart where to store or apply the table.
        """
        table = KM.PiecewiseLinearTable()
        minimum_columns = max(self.first_column_id, self.second_column_id) + 1
        row_id = self.skiprows
        with open(self.filename, 'r') as table_file:
            data = csv.reader(table_file, delimiter=self.delimiter, skipinitialspace=True)
            for _ in range(self.skiprows):
                next(data)
            for row in data:
                row_id += 1
                if row:  # skip empty rows
                    if len(row) < minimum_columns:
                        msg = self.__class__.__name__ + ". {}\n".format(self.filename)
                        msg += "There is not enough data, a {}-column row is found at line {}.\n".format(len(row), row_id)
                        msg += "In order to read the columns {} and {}, the table must have at least {} columns.".format(self.first_column_id, self.second_column_id, minimum_columns)
                        raise Exception(msg)
                    table.AddRow(self._Float(row[self.first_column_id], row_id), self._Float(row[self.second_column_id], row_id))
        if self.table_id > -1:
            if model_part != None:
                model_part.AddTable(self.table_id, table)
            else:
                err_msg = "Asking to save table with id {} but no model part is provided.".format(self.table_id)
                raise Exception(err_msg)
        return table

    def _Float(self, value, row_id):
        try:
            return float(value)
        except ValueError:
            KM.Logger.PrintWarning(self.__class__.__name__, "{} replaced by {} at row {}".format(value, self.na_replace, row_id))
            return self.na_replace

