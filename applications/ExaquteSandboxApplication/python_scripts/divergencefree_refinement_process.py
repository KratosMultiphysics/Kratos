import KratosMultiphysics
import KratosMultiphysics.ExaquteSandboxApplication as ExaquteSandboxApplication
import KratosMultiphysics.MeshingApplication as MeshingApplication

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DivergenceFreeRefinementProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class DivergenceFreeRefinementProcess(KratosMultiphysics.Process):
    def __init__(self,model,settings):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"                       : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "integration_end_point_control_value"   : 1e6,
                "integration_start_point_control_value" : 0.0,
                "average_settings": {
                    "model_part_name"                               : "PLEASE_SPECIFY_MODEL_PART_NAME",
                    "variables_list"                                : ["DIVERGENCE"],
                    "averaged_variables_list"                       : ["AVERAGED_DIVERGENCE"],
                    "time_averaging_container"                      : "ElementalNonHistorical",
                    "time_averaging_method"                         : "RootMeanSquare",
                    "integration_start_point_control_variable_name" : "TIME",
                    "integration_start_point_control_value"         : 0.2
                },
                "refinement_settings" : {
                    "refinement_strategy"     : "global_tolerance_strategy",
                    "reference_variable_name" : "AVERAGED_DIVERGENCE",
                    "minimal_size"            : 0.1,
                    "maximal_size"            : 10.0,
                    "mean_distribution_strategy":
                    {
                        "target_refinement_coefficient"       : 0.9,
                        "refinement_bound"                    : 2.0,
                        "reference_norm_name"                 : "VELOCITY_H1_SEMINORM"
                    },
                    "maximum_strategy":
                    {
                        "target_refinement_coefficient"       : 0.1,
                        "refinement_coefficient"              : 2.0
                    },
                    "global_tolerance_strategy":
                    {
                        "global_tolerance"                 : 0.1
                    }
                }
            }  """ )
        settings.ValidateAndAssignDefaults(default_parameters)
        self.integration_end_time = settings["integration_end_point_control_value"].GetDouble()
        self.model = model
        self.model_part_name = settings["model_part_name"].GetString()
        self.domain_size = self.model.GetModelPart(self.model_part_name).ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        # Divergence process
        # ready

        # Averaging process
        self.averaging_process_parameters = settings["average_settings"]

        # Metric construction process
        self.metric_divergencefree_parameters = settings["refinement_settings"]

    def Check(self):
        # Divergence process

        # Averaging process
        KratosMultiphysics.TimeAveragingProcess(self.model,self.averaging_process_parameters).Check()
        # Metric construction process

    def ExecuteInitialize(self):
        ExaquteSandboxApplication.DivergenceProcess(self.model.GetModelPart(self.model_part_name)).ExecuteInitialize()
        KratosMultiphysics.TimeAveragingProcess(self.model,self.averaging_process_parameters).ExecuteInitialize()
        if (self.domain_size == 2):
            ExaquteSandboxApplication.MetricDivergenceFreeProcess2D(\
                self.model.GetModelPart(self.model_part_name),self.metric_divergencefree_parameters).ExecuteInitialize()
        elif (self.domain_size == 3):
            ExaquteSandboxApplication.MetricDivergenceFreeProcess3D(\
                self.model.GetModelPart(self.model_part_name),self.metric_divergencefree_parameters).ExecuteInitialize()

    def ExecuteFinalizeSolutionStep(self):
        current_time = self.model.GetModelPart(self.model_part_name).ProcessInfo[KratosMultiphysics.TIME]
        if (current_time <= self.integration_end_time):
            ExaquteSandboxApplication.DivergenceProcess(self.model.GetModelPart(self.model_part_name)).ExecuteBeforeOutputStep()
            KratosMultiphysics.TimeAveragingProcess(self.model,self.averaging_process_parameters).ExecuteFinalizeSolutionStep()

    def ExecuteFinalize(self):
        # Calculate nodal h
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.model.GetModelPart(self.model_part_name))
        find_nodal_h.Execute()
        # Calculate divergence metric of the variable
        if (self.domain_size == 2):
            divergencefree_metric = ExaquteSandboxApplication.MetricDivergenceFreeProcess2D(\
                self.model.GetModelPart(self.model_part_name),self.metric_divergencefree_parameters)
        elif (self.domain_size == 3):
            divergencefree_metric = ExaquteSandboxApplication.MetricDivergenceFreeProcess3D(\
                self.model.GetModelPart(self.model_part_name),self.metric_divergencefree_parameters)
        divergencefree_metric.Execute()
        # Execute remeshing process
        remesh_parameters = KratosMultiphysics.Parameters()
        if (self.domain_size == 2):
            MmgProcess = MeshingApplication.MmgProcess2D(self.model.GetModelPart(self.model_part_name),remesh_parameters)
        elif (self.domain_size == 3):
            MmgProcess = MeshingApplication.MmgProcess3D(self.model.GetModelPart(self.model_part_name),remesh_parameters)
        MmgProcess.Execute()