import KratosMultiphysics as KM
import csv

class ReadCsvTableUtility:
    """This class reads retrieves a table from the specified parameters."""

    def __init__(self, settings):
        """Constructor of the csv table reader: validate the parameters."""
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
        """Read a csv table."""
        table = KM.PiecewiseLinearTable()
        with open(self.file_name, 'r') as table_file:
            data = csv.reader(table_file, delimiter=self.delimiter, skipinitialspace=True)
            for _ in range(self.skiprows):
                next(data)
            for row in data:
                if row:  # skip empty rows
                    table.AddRow(self._Float(row[0]), self._Float(row[1]))
        if self.table_id > -1:
            model_part.AddTable(table)
        return table

    def _Float(self, value):
        try:
            return float(value)
        except ValueError:
            return self.na_replace

