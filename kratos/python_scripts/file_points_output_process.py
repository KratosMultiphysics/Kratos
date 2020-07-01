from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import csv, os

# Importing the Kratos Library
import KratosMultiphysics

# other imports
from KratosMultiphysics.multiple_points_output_process import MultiplePointsOutputProcess

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return FilePointsOutputProcess(Model, settings["Parameters"])

class FilePointsOutputProcess(KratosMultiphysics.Process):
    """This process writes output for several points given in a file to an output file
    Internally it holds an object of type "MultiplePointsOutputProcess"
    Usage:
        - file path of an .csv-file (!) for the import of sampling points need to be specified
        - the first line of the .csv-file should state the column names: "x,y,z"
        - the following rows should each be one point with its x,y and z value in the respective column of the .csv-file, separated by a comma 
        - the number of sampling points needs to be specified
    """
    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters('''{
            "help"              : "This process writes output for several points given in a file to a file. Internally it holds an object of type MultiplePointsOutputProcess",
            "model_part_name"   : "",
            "entity_type"       : "element",
            "file_path"         : "",
            "sampling_points"   : 3,
            "output_variables"  : [],
            "historical_value"  : true,
            "search_tolerance"  : 1e-6,
            "print_format"      : "",
            "output_file_settings": {}
        }''')

        params.ValidateAndAssignDefaults(default_settings)

        # retrieving file path and name
        file_path = params["file_path"].GetString()
        if len(file_path) == 0:
            raise Exception('A file path has has to be provided!')

        # get and check number of sampling points
        number_of_sampling_points = params["sampling_points"].GetInt()
        if number_of_sampling_points <= 0:
            raise Exception('The number of sampling points has to be larger than 0!')

        # check the entity type
        entity_type = params["entity_type"].GetString()
        if not entity_type == "element" or entity_type == "condition":
            raise Exception('"entity_type" can be either "element" or "condition"!')

        # initialize position matrix for sampling points
        positions = KratosMultiphysics.Matrix(number_of_sampling_points, 3)

        # check and read file
        if os.path.exists(file_path) and os.access(file_path, os.R_OK):

            # open file
            try:
                with open(file_path, newline=None) as file:
                    reader = csv.DictReader(file, delimiter=',')
                    for p_count, row in enumerate(reader):

                        # check if number of sampling points is already reached
                        if p_count == number_of_sampling_points:
                            KratosMultiphysics.Logger.PrintWarning("FilePointsOutputProcess: Number of points in file higher than number of samples!")
                            break

                        # add position of sampling point
                        positions[p_count,0] = float(row['x'])
                        positions[p_count,1] = float(row['y'])
                        positions[p_count,2] = float(row['z'])

            except:
                raise Exception('The .csv-file can not be read (correctly)!')
        else:
            raise Exception('The file does not exist or can not be accessed!')

        # check if enough points were found
        p_count += 1
        if p_count < number_of_sampling_points:
            KratosMultiphysics.Logger.PrintWarning("FilePointsOutputProcess: Found " + str(p_count) + " points in file. Number of samples given was " + str(number_of_sampling_points) + "!")

        params.RemoveValue("file_path")
        params.RemoveValue("file_name")
        params.RemoveValue("file_type")
        params.RemoveValue("sampling_points")

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
