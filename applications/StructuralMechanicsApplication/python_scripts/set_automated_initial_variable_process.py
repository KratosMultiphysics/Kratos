# Importing the Kratos Library
from pathlib import Path
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics import Logger
from KratosMultiphysics.read_csv_table_utility import ReadCsvTableUtility

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    
    """
    This process reads tables in csv format containing the strain/stress components as function of the spatial radial coordinate from the edge of a given hole.
    Tables are defined by layers ortogonally to the radial direction (generatrix hole direction). Each layer of elements must be allocated in a modelpart.
    Stress components provided by the tables (linear picewise tables) are interpolated according to each element centroid. 
    Interpolated components are transferred to the initial variable (INITIAL_STRAIN_VECTOR/INITIAL_STRESS_VECTOR).
    A minimum of one table per layer must be provided to the process. 
    Missing stress component are considered to be zero in the INITIAL_STRAIN_VECTOR/INITIAL_STRESS_VECTOR variables. 
    Tables corresponding to out of range strain/stress components are neglect (e.g. Layer1_StresComp7, Layer2_StresComp8).
    """

    default_settings = KM.Parameters(
        """{
            "help"                     : "This automates the application of initial strain/stress variables using csv tables",
            "model_part_name"          : "please_specify_model_part_name",
            "variable_name"            : "SPECIFY_VARIABLE_NAME",
            "hole_generatrix_axis"     : [0.0,0.0,1.0],
            "hole_generatrix_point"    : [0.0,0.0,0.0],
            "hole_radius_offset"       : 0.0,
            "initial_variable_table"     : {
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
    process_settings["initial_variable_table"].AddValue("table_id", default_table_id)

    file_path = Path(process_settings["initial_variable_table"]["filename"].GetString())
    layer_name = file_path.name.split("_")[0]

    layer_list = [file for file in file_path.parent.iterdir() if file.name.split("_")[0] == layer_name]

    if not layer_list:
        ErrorMsg = f"Tables of {layer_name} not found"
        raise RuntimeError(ErrorMsg)
    else:
        component_list = []
        table_id_list = []
        out_of_range_component_list = []
        for i in range (0, len(layer_list)):
            component_number = int(layer_list[i].stem.split("_")[1][-1])
            if  0 < component_number < 7:
                component_list.append(component_number)
                table_id = int("".join(layer_number for layer_number in layer_name if layer_number.isdigit()) + str(component_number - 1))
                table_id_list.append(table_id)
                process_settings["initial_variable_table"]["filename"].SetString(layer_list[i].as_posix())
                process_settings["initial_variable_table"]["table_id"].SetInt(table_id)
                ReadCsvTableUtility(process_settings["initial_variable_table"]).Read(computing_model_part)
            else:
                out_of_range_component_list.append(component_number)

    raw_variable_name = process_settings["variable_name"].GetString()
    variable_name = raw_variable_name.split("_")[0] + " " + raw_variable_name.split("_")[1].split("_")[0]

    if len(out_of_range_component_list) != 0:
        if len(out_of_range_component_list) == 1:
            Logger.PrintWarning("SetAutomatedInitialVariableProcess:: ", f"{variable_name.capitalize()} component {out_of_range_component_list[0]} of {layer_name} is out of range. The corresponding table will be neglected")
        else:
            out_of_range_component_name = ", ".join(str(out_of_range_component) for (out_of_range_component) in out_of_range_component_list)
            Logger.PrintWarning("SetAutomatedInitialVariableProcess:: ", f"{variable_name.capitalize()} components {out_of_range_component_name} of {layer_name} are out of range. Correspoding tables will be negleted")   

    if len(component_list) < 6:
        missing_component_list = list(set(range(1,7)).difference(component_list))
        if len(missing_component_list) == 1:
            Logger.PrintWarning("SetAutomatedInitialVariableProcess:: ", f"Table correspoding to {variable_name.lower()} component {str(missing_component_list[0])} of {layer_name} not found. A zero entry will be added to the {raw_variable_name} variable")
        else:
            missing_component_name = ", ".join(str(missing_component) for (missing_component) in missing_component_list)
            Logger.PrintWarning("SetAutomatedInitialVariableProcess:: ", f"Tables correspoding to {variable_name.lower()} components {missing_component_name} of {layer_name} not found. Zero entries will be added to the {raw_variable_name} variable")
    else:
        Logger.PrintInfo("SetAutomatedInitialVariableProcess:: ", f"{variable_name.capitalize()} tables of {layer_name} were successfully imported")
    
    default_table_id_vector = KM.Parameters("""{
    "table_id_vector": [10,11,12,13,14,15]
    }""")
    process_settings.AddEmptyValue("table_id_vector").SetVector(table_id_list)

    process_settings.RemoveValue("help")
    process_settings.RemoveValue("model_part_name")
    process_settings.RemoveValue("initial_variable_table")

    return SMA.SetAutomatedInitialVariableProcess(computing_model_part, process_settings)