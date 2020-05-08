import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyEmbeddedSkinVisualizationProcess(Model, settings["Parameters"])

class ApplyEmbeddedSkinVisualizationProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        # Backwards compatibility helper
        if settings.Has("parallel_type"):
            settings.RemoveValue("parallel_type")
            warn_msg = "\"parallel_type\" is no longer needed. Removing from input settings."
            KratosMultiphysics.Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)
        if settings.Has("output_configuration"):
            settings.RemoveValue("output_configuration")
            warn_msg = "\"output_configuration\" is no longer needed as the visualization mesh is no longer print in the \"ApplyEmbeddedSkinVisualizationProcess\". Removing from input settings.\n"
            warn_msg += "Add your preferred output process in the \"output_process_list\" of the simulation settings (namely ProjectParameters.json)."
            KratosMultiphysics.Logger.PrintWarning('\n\x1b[1;31mDEPRECATION-WARNING\x1b[0m', warn_msg)

        # Create the embedded skin visualization proces
        self.EmbeddedSkinVisualizationProcess = KratosFluid.EmbeddedSkinVisualizationProcess(
            Model,
            settings)

    def ExecuteInitialize(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteFinalize()
