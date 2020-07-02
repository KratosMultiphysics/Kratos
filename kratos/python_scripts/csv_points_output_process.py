import csv, os

# Importing the Kratos Library
import KratosMultiphysics

# other imports
from KratosMultiphysics.multiple_points_output_process import MultiplePointsOutputProcess

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return CSVPointsOutputProcess(Model, settings["Parameters"])

class CSVPointsOutputProcess(KratosMultiphysics.Process):
    """This process writes output for several points given in a .csv-file to an output file
    Internally it holds an object of type "MultiplePointsOutputProcess"
    Usage:
        - file path of an .csv-file for the import of sampling points need to be specified
        - the first line of the .csv-file should state the column names: "x,y,z"
        - the following rows should each be one point with its x,y and z value in the respective column of the .csv-file, separated by a comma 
    """
    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters('''{
            "help"              : "This process writes output for several points given in a .csv-file to a file. Internally it holds an object of type MultiplePointsOutputProcess",
            "model_part_name"   : "",
            "entity_type"       : "element",
            "csv_file_path"         : "",
            "output_variables"  : [],
            "historical_value"  : true,
            "search_tolerance"  : 1e-6,
            "print_format"      : "",
            "output_file_settings": {}
        }''')

        params.ValidateAndAssignDefaults(default_settings)

        # retrieving file path and name
        csv_file_path = params["csv_file_path"].GetString()
        if not csv_file_path:
            raise Exception('A file path has has to be provided!')

        # initialize array for sampling points
        points = []
        p_count = 0

        # open file and get points
        with open(csv_file_path, newline=None) as file:
            reader = csv.DictReader(file, delimiter=',')
            for p_count, row in enumerate(reader):
                points.append([float(row['x']), float(row['y']), float(row['z'])])
        
        # initialize position matrix for sampling points
        p_count += 1
        positions = KratosMultiphysics.Matrix(p_count, 3)
        
        # add position of sampling point
        for c in range(0, p_count):
            positions[c, 0] = points[c][0]
            positions[c, 1] = points[c][1]
            positions[c, 2] = points[c][2]

        params.RemoveValue("csv_file_path")
        params.RemoveValue("file_type")

        params.AddEmptyValue("positions")
        params["positions"].SetMatrix(positions)

        self.multiple_points_output_process = MultiplePointsOutputProcess(model, params)

    def ExecuteInitialize(self):
        self.multiple_points_output_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.multiple_points_output_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.multiple_points_output_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.multiple_points_output_process.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        self.multiple_points_output_process.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        self.multiple_points_output_process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        self.multiple_points_output_process.ExecuteFinalize()
