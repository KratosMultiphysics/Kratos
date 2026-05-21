# Importing the Kratos Library
from pathlib import Path
import KratosMultiphysics as KM
import KratosMultiphysics.ConstitutiveLawsApplication as CLA
from KratosMultiphysics import Logger
from KratosMultiphysics.read_csv_table_utility import ReadCsvTableUtility

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")

    """
    This process reads tables in csv format containing damage values as function of the spatial radial coordinate from the edge of a given hole.
    Tables are defined by layers ortogonally to the radial direction (generatrix hole direction). Each layer of elements must be allocated in a modelpart.
    Damage values provided by the tables (linear picewise tables) are interpolated according to each element centroid.
    Interpolated damage values are set to the variable DAMAGE.
    Each layer requires one single table.

    """

    default_settings = KM.Parameters(
        """{
            "help"                     : "This automates the application of initial damage values using csv tables",
            "model_part_name"          : "please_specify_model_part_name",
            "hole_generatrix_axis"     : [0.0,0.0,1.0],
            "hole_generatrix_point"    : [0.0,0.0,0.0],
            "hole_radius_offset"       : 0.0,
            "initial_damage_table"     : {
                        "name"             : "csv_table",
                        "filename"         : "sample.csv",
                        "delimiter"        : ",",
                        "skiprows"         : 1,
                        "first_column_id"  : 0,
                        "second_column_id" : 1,
                        "na_replace"       : 0.0
                    }
        }""")
    process_settings = settings["Parameters"]
    process_settings.ValidateAndAssignDefaults(default_settings)
    computing_model_part = Model[process_settings["model_part_name"].GetString()]

    default_table_id = KM.Parameters("""{
    "table_id": 0
    }""")
    process_settings["initial_damage_table"].AddValue("table_id", default_table_id)

    file_path = Path(process_settings["initial_damage_table"]["filename"].GetString())
    layer_name = file_path.name.split("_")[0]

    table_id = int("".join(layer_number for layer_number in layer_name if layer_number.isdigit()))
    process_settings["initial_damage_table"]["table_id"].SetInt(table_id)
    ReadCsvTableUtility(process_settings["initial_damage_table"]).Read(computing_model_part)

    Logger.PrintInfo("SetAutomatedInitialVariableProcess:: ", f"Table of {layer_name} was successfully imported")

    process_settings.RemoveValue("help")
    process_settings.RemoveValue("model_part_name")
    process_settings.RemoveValue("initial_damage_table")

    process_settings.AddEmptyValue("table_id").SetInt(table_id)

    return CLA.SetAutomatedInitialDamageProcess(computing_model_part, process_settings)