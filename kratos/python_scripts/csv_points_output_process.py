# Importing the Kratos Library
import KratosMultiphysics

# other imports
from KratosMultiphysics.multiple_points_output_process import MultiplePointsOutputProcess
import csv

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    params = settings["Parameters"]

    default_settings = KratosMultiphysics.Parameters('''{
        "help"                 : "This process writes output for several points given in a .csv-file ('csv_file_path') to a file. Internally it holds an object of type MultiplePointsOutputProcess. Usage: The first line of the .csv-file should state the column names: 'x,y,z'. The following rows should each be one point with its x,y and z value in the respective column of the .csv-file, separated by a comma.",
        "model_part_name"      : "",
        "entity_type"          : "element",
        "interval"             : [0.0, 1e30],
        "csv_file_path"        : "",
        "output_variables"     : [],
        "historical_value"     : true,
        "search_tolerance"     : 1e-6,
        "print_format"         : "",
        "output_file_settings" : {}
    }''')

    params.ValidateAndAssignDefaults(default_settings)

    # retrieving file path and name
    csv_file_path = params["csv_file_path"].GetString()
    if not csv_file_path:
        raise Exception('A file path has has to be provided!')

    # initialize array for sampling points
    points = []

    # open file and get points
    with open(csv_file_path, mode='r', newline=None) as csv_file:
        reader = csv.DictReader(csv_file, delimiter=',')
        points = [(float(row['x']), float(row['y']), float(row['z'])) for row in reader]

    # initialize position matrix for sampling points
    positions = KratosMultiphysics.Matrix(len(points), 3)
    
    # add position of sampling point
    for i, coords in enumerate(points):
        positions[i, 0] = coords[0]
        positions[i, 1] = coords[1]
        positions[i, 2] = coords[2]

    # adapt process parameters to MultiplePointsOutputProcess
    params.RemoveValue("csv_file_path")
    params.AddEmptyValue("positions")
    params["positions"].SetMatrix(positions)

    # return and instance of MultiplePointsOutputProcess with all sampling positions from .csv-file
    return MultiplePointsOutputProcess(Model, params)
