# Importing the Kratos Library
import KratosMultiphysics

# other imports
from KratosMultiphysics.multiple_points_output_process import MultiplePointsOutputProcess
import csv, os

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    params = settings["Parameters"]

    default_settings = KratosMultiphysics.Parameters('''{
        "help"              : "This process writes output for several points given in a .csv-file ('csv_file_path') to a file. Internally it holds an object of type MultiplePointsOutputProcess. Usage: The first line of the .csv-file should state the column names: 'x,y,z'. The following rows should each be one point with its x,y and z value in the respective column of the .csv-file, separated by a comma.",
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

    # adapt process parameters to MultiplePointsOutputProcess
    params.RemoveValue("csv_file_path")
    params.RemoveValue("file_type")
    params.AddEmptyValue("positions")
    params["positions"].SetMatrix(positions)

    # return and instance of MultiplePointsOutputProcess with all sampling positions from .csv-file
    return MultiplePointsOutputProcess(Model, params)
