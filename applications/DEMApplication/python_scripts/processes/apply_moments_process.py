import KratosMultiphysics

# Import applications
import KratosMultiphysics.DEMApplication as DEM

# Other imports

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    process_settings = settings["Parameters"]

    folder_settings = KratosMultiphysics.Parameters("""{
            "help"                 : "This process applies loads over the particles in a certain submodelpart, for a certain time interval",
            "mesh_id"              : 0,
            "model_part_name"      : "please_specify_model_part_name",
            "variable_name"        : "object",
            "constrained"          : [true,true,true],
                "value"            : [10.0, "3*t", "x+y"],
                "table"            : [0, 0, 0],
            "interval"             : [0.0, 1e30]
        }""" )

    # Detect "End" as a tag and replace it by a large number
    if process_settings.Has("interval"):
        if process_settings["interval"][1].IsString():
            if process_settings["interval"][1].GetString() == "End":
                process_settings["interval"][1].SetDouble(1e30) # = default_settings["interval"][1]
            else:
                raise Exception("The second value of interval can be \"End\" or a number, interval currently:"+settings["interval"].PrettyPrintJsonString())


    process_settings.AddMissingParameters(folder_settings)

    if process_settings.Has("model_part_name"):
        computing_model_part = Model[process_settings["model_part_name"].GetString()]
    else: # using default name
        computing_model_part = Model["DEM"]

    process_settings.RemoveValue("help")

    return DEM.ApplyMomentsProcess(computing_model_part, process_settings)
