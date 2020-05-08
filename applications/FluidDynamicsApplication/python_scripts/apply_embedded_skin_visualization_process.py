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

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                     : "origin_model_part",
            "visualization_model_part_name"       : "origin_model_part_visualization",
            "shape_functions"                     : "standard",
            "reform_model_part_at_each_time_step" : false,
            "visualization_variables"             : ["VELOCITY","PRESSURE"]
        } """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        # # Get the origin model part
        # self.origin_model_part = Model[settings["model_part_name"].GetString()]

        # # Set up the visualization model part
        # visualization_buffer_size = 1
        # self.visualization_model_part = Model.CreateModelPart(settings["visualization_model_part_name"].GetString(), visualization_buffer_size)
        # self.visualization_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.origin_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # # Add the visualization variables to the visualization model part.
        # for i_var in range(0, settings["visualization_variables"].size()):
        #     variable_name = settings["visualization_variables"][i_var].GetString()
        #     self.visualization_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(variable_name))

        # Create the embedded skin visualization proces
        self.EmbeddedSkinVisualizationProcess = KratosFluid.EmbeddedSkinVisualizationProcess(
            Model,
            settings)

    def ExecuteInitialize(self):
        # self.gid_output.ExecuteInitialize()
        self.EmbeddedSkinVisualizationProcess.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteBeforeSolutionLoop()
        # self.gid_output.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteInitializeSolutionStep()
        # self.gid_output.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteFinalizeSolutionStep()
        # self.gid_output.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteBeforeOutputStep()
        # if (self.gid_output.IsOutputStep()):
        #     self.gid_output.PrintOutput()

    def ExecuteAfterOutputStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteAfterOutputStep()
        # self.gid_output.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteFinalize()
        # self.gid_output.ExecuteFinalize()
