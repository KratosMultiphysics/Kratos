# Project imports
import KratosMultiphysics as KM

# System imports
import csv
from typing import Optional # Optional
import pathlib # pathlib.Path


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

    def __init__(self, settings: KM.Parameters):
        """Constructor of the csv table reader: validate the parameters.

        Keyword arguments:
        self -- It signifies an instance of the class.
        settings -- Kratos parameters containing solver settings.
        """
        settings.ValidateAndAssignDefaults(self.GetSchema())

        self.filename: pathlib.Path = pathlib.Path(settings["filename"].GetString())
        self.delimiter: str = settings["delimiter"].GetString()
        self.skiprows: int = settings["skiprows"].GetInt()
        self.first_column_id: int = settings["first_column_id"].GetInt()
        self.second_column_id: int = settings["second_column_id"].GetInt()
        self.table_id: int = settings["table_id"].GetInt()
        self.na_replace: float = settings["na_replace"].GetDouble()

    def Read(self,
             model_part: Optional[KM.ModelPart] = None):
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

    def _Float(self, value: float, row_id: int) -> float:
        try:
            return float(value)
        except ValueError:
            KM.Logger.PrintWarning(self.__class__.__name__, "{} replaced by {} at row {}".format(value, self.na_replace, row_id))
            return self.na_replace

    @classmethod
    def GetSchema(cls) -> KM.Schema:
        return KM.Schema(KM.Parameters("""{
            "title" : "ReadCsvTableUility",
            "description" : "A utility class for reading data from a CSV file.",
            "properties" : {
                "name" : {
                    "description" : "Name of the table to refer to internally.",
                    "type" : "string",
                    "default" : "csv_table"
                },
                "filename" : {
                    "description" : "Path pointing to the CSV file.",
                    "type" : "string"
                },
                "delimiter" : {
                    "type" : "string",
                    "pattern" : "[,;\\t|]",
                    "default" : ","
                },
                "skiprows" : {
                    "description" : "Number of rows to skip before beginning to read data.",
                    "type" : "integer",
                    "default" : 0
                },
                "first_column_id" : {
                    "description" : "Zero-based index of the first column to read.",
                    "type" : "integer",
                    "minimum" : 0,
                    "default" : 0
                },
                "second_column_id" : {
                    "description" : "Zero-based index of the second column to read.",
                    "type" : "integer",
                    "minimum" : 0,
                    "default" : 1
                },
                "table_id" : {
                    "description" : "Identifier of the table to store by in the provided model part. If -1, it will not be stored.",
                    "type" : "integer",
                    "minimum" : -1,
                    "default" : -1
                },
                "na_replace" : {
                    "description" : "Value to replace NaNs by.",
                    "type" : "number",
                    "default" : 0.0
                }
            },
            "additionalProperties" : false,
            "required" : ["filename"]
        }"""))

