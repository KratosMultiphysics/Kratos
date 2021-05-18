import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    # Backwards compatibility helper
    if settings["Parameters"].Has("parallel_type"):
        settings["Parameters"].RemoveValue("parallel_type")
        warn_msg = "\"parallel_type\" is no longer needed. Removing from input settings."
        KratosMultiphysics.Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)
    if settings["Parameters"].Has("output_configuration"):
        settings["Parameters"].RemoveValue("output_configuration")
        warn_msg = "\"output_configuration\" is no longer needed as the visualization mesh is no longer print in the \"ApplyEmbeddedSkinVisualizationProcess\". Removing from input settings.\n"
        warn_msg += "Add your preferred output process in the \"output_process_list\" of the simulation settings (namely ProjectParameters.json)."
        KratosMultiphysics.Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)

    return KratosFluid.EmbeddedSkinVisualizationProcess(Model, settings["Parameters"])