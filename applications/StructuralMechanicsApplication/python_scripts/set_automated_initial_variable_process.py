# Importing the Kratos Library
from genericpath import isfile
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics import Logger
from KratosMultiphysics.read_csv_table_utility import ReadCsvTableUtility

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")

    default_settings = KM.Parameters(
        """{
            "help"                     : "This automates the application of initial strain/stress variables using csv tables",
            "model_part_name"          : "please_specify_model_part_name",
            "variable_name"            : "SPECIFY_VARIABLE_NAME",
            "hole_generatrix_axis"     : [0.0,0.0,1.0],
            "hole_generatrix_point"    : [0.0,0.0,0.0],
            "hole_radius_offset"       : 0.0,
            "last_layer"               : false,
            "initial_variable_table"     : {
                        "name"             : "csv_table",
                        "filename"         : "sample.csv",
                        "delimiter"        : ",",
                        "skiprows"         : 1,
                        "first_column_id"  : 0,
                        "second_column_id" : 1,
                        "table_id"         : -1,
                        "na_replace"       : 0.0
                    }
        }""")
    process_settings = settings["Parameters"]
    process_settings.ValidateAndAssignDefaults(default_settings)
    computing_model_part = Model[process_settings["model_part_name"].GetString()]

    if process_settings["initial_variable_table"]["filename"].GetString().find("/") != -1:
        filepath = process_settings["initial_variable_table"]["filename"].GetString().split("/")[0] + "/"
        filename = process_settings["initial_variable_table"]["filename"].GetString().split("/")[1]
    else:
        filepath = ""
        filename = process_settings["initial_variable_table"]["filename"].GetString()

    layer_string = filename.split("_")[0]
    layer_number_string = layer_string[-1]
    stress_string = filename.split("_")[1].split(".")[0]
    variable_string = process_settings["variable_name"].GetString().split("_")[0] + " " + process_settings["variable_name"].GetString().split("_")[1].split("_")[0]

    i = 0

    if not isfile(filepath + layer_string[:-1] + str(i + 1) + "_" + stress_string + ".csv"):
        ErrorMsg = "Table " + "\"" + layer_string[:-1] + str(i + 1) + "_" + stress_string + ".csv\" not found"
        raise RuntimeError(ErrorMsg)

    while isfile(filepath + layer_string + "_" + stress_string[:-1] + str(i + 1) + ".csv") and i < 6:
        
        table_id = int(layer_number_string + str(i))

        process_settings["initial_variable_table"]["filename"].SetString(filepath + layer_string + "_" + stress_string[:-1] + str(i+1)+".csv")
        process_settings["initial_variable_table"]["table_id"].SetInt(table_id)

        ReadCsvTableUtility(process_settings["initial_variable_table"]).Read(computing_model_part)
        
        if not isfile(filepath + layer_string + "_" + stress_string[:-1] + str(i + 2) + ".csv") or i == 5:
            Logger.PrintInfo("SetAutomatedInitialVariableProcess:: ", variable_string.capitalize() + " tables of " + layer_string + " were successfully imported")
     
            if not isfile(filepath + layer_string[:-1] + str(int(layer_number_string) + 1) + "_" + stress_string + ".csv") and not process_settings["last_layer"].GetBool():
                ErrorMsg = "Table " + "\"" + layer_string[:-1] + str(int(layer_number_string) + 1) + "_" + stress_string + ".csv\" not found"
                raise RuntimeError(ErrorMsg)

        i += 1
    
    return SMA.SetAutomatedInitialVariableProcess(computing_model_part, process_settings)