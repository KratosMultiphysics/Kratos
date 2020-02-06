import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyEmbeddedSkinVisualizationProcess(Model, settings["Parameters"])

class ApplyEmbeddedSkinVisualizationProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "parallel_type"                       : "OpenMP",
            "model_part_name"                     : "origin_model_part",
            "visualization_model_part_name"       : "origin_model_part_visualization",
            "shape_functions"                     : "standard",
            "reform_model_part_at_each_time_step" : false,
            "visualization_variables"             : ["VELOCITY","PRESSURE"],
            "output_configuration"                : {
                "result_file_configuration" : {
                    "gidpost_flags" : {
                        "GiDPostMode"           : "GiD_PostBinary",
                        "WriteDeformedMeshFlag" : "WriteDeformed",
                        "WriteConditionsFlag"   : "WriteConditions",
                        "MultiFileFlag"         : "SingleFile"
                    },
                    "file_label"          : "time",
                    "output_control_type" : "time",
                    "output_frequency"    : 0.1,
                    "body_output"         : true,
                    "node_output"         : false,
                    "skin_output"         : false,
                    "nodal_results"       : []
                },
                "point_data_configuration"  : []
            }
        } """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        # Get the origin model part
        self.origin_model_part = Model[settings["model_part_name"].GetString()]

        # Set up the visualization model part
        visualization_buffer_size = 1
        self.visualization_model_part = Model.CreateModelPart(settings["visualization_model_part_name"].GetString(), visualization_buffer_size)
        self.visualization_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.origin_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # Check that the nodal results array is empty
        if (settings["output_configuration"]["result_file_configuration"]["nodal_results"].size() != 0):
            error_msg = "The nodal_results field in output_configuration is not empty.\n Add the variables in the visualization_variables field instead."
            raise Exception(error_msg)

        # Add the visualization model part variables to the visualization model part.
        # Add them to the nodal_results GiD output process list as well.
        for i_var in range(0, settings["visualization_variables"].size()):
            variable_name = settings["visualization_variables"][i_var].GetString()
            settings["output_configuration"]["result_file_configuration"]["nodal_results"].Append(variable_name)
            self.visualization_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(variable_name))

        # Set an auxilar Kratos Parameters object to build the skin visualization process
        aux_params = KratosMultiphysics.Parameters("""{}""")
        aux_params.AddValue("shape_functions",settings["shape_functions"])
        aux_params.AddValue("visualization_variables",settings["visualization_variables"])
        aux_params.AddValue("reform_model_part_at_each_time_step",settings["reform_model_part_at_each_time_step"])

        self.EmbeddedSkinVisualizationProcess = KratosFluid.EmbeddedSkinVisualizationProcess(
            self.origin_model_part,
            self.visualization_model_part,
            aux_params)

        # Set the output variables and build the GiD output process
        if (settings["parallel_type"].GetString() == "OpenMP"):
            from KratosMultiphysics.gid_output_process import GiDOutputProcess
            self.gid_output = GiDOutputProcess(self.visualization_model_part, settings["visualization_model_part_name"].GetString(), settings["output_configuration"])
        elif (settings["parallel_type"].GetString() == "MPI"):
            from KratosMultiphysics.mpi.distributed_gid_output_process import DistributedGiDOutputProcess
            self.gid_output = DistributedGiDOutputProcess(self.visualization_model_part, settings["visualization_model_part_name"].GetString(), settings["output_configuration"])

    def ExecuteInitialize(self):
        self.gid_output.ExecuteInitialize()
        self.EmbeddedSkinVisualizationProcess.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteBeforeSolutionLoop()
        self.gid_output.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteInitializeSolutionStep()
        self.gid_output.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteFinalizeSolutionStep()
        self.gid_output.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteBeforeOutputStep()
        if (self.gid_output.IsOutputStep()):
            self.gid_output.PrintOutput()

    def ExecuteAfterOutputStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteAfterOutputStep()
        self.gid_output.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteFinalize()
        self.gid_output.ExecuteFinalize()
