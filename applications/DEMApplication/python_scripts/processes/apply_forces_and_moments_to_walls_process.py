import KratosMultiphysics

# Import applications
import KratosMultiphysics.DEMApplication as DEM

# Other imports

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    process_settings = settings["Parameters"]

    folder_settings = KratosMultiphysics.Parameters("""{
            "help"                 : "This process applies loads over the rigid walls in a certain submodelpart, for a certain time interval",
            "model_part_name"      : "please_specify_model_part_name",
            "force_settings" : {
                "value"            : [10.0, "3*t", "x+y"],
                "table"            : [0, 0, 0]
            },
            "moment_settings" : {
                "value"            : [10.0, "3*t", "x+y"],
                "table"            : [0, 0, 0]
            },
            "interval"             : [0.0, 1e30]
        }""" )

    process_settings.AddMissingParameters(folder_settings)

    if process_settings.Has("model_part_name"):
        computing_model_part = Model[process_settings["model_part_name"].GetString()]
    else: # using default name
        computing_model_part = Model["DEM_FEM_boundary"]

    process_settings.RemoveValue("help")

    return DEM.ApplyForcesAndMomentsToWallsProcess(computing_model_part, process_settings)
