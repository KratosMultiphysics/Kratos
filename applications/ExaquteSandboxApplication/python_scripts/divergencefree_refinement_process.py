import KratosMultiphysics
import KratosMultiphysics.ExaquteSandboxApplication as KratosExaqute
import KratosMultiphysics.MeshingApplication as KratosMeshing

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DivergenceFreeRefinementProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class DivergenceFreeRefinementProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name":"PLEASE_SPECIFY_MODEL_PART_NAME",
                "time_average_length" : 10.0,
                "time_average_start_coefficient" : 0.2,
                "refinement_settings" : {
                    "refinement_strategy" : "global_tolerance_strategy",
                    "reference_variable_name" : "DIVERGENCE_WEIGHTED",
                    "minimal_size"                     : 0.1,
                    "maximal_size"                     : 10.0,
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

        self.fluid_model_part = model[settings["model_part_name"].GetString()]
        self.domain_size = self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        # Weighted divergence process
        self.time_average_length = settings["time_average_length"].GetDouble()
        self.weighted_divergencefree_parameters = KratosMultiphysics.Parameters()
        self.weighted_divergencefree_parameters.AddValue("time_coefficient",settings["time_average_start_coefficient"])

        # Metric construction process
        self.metric_divergencefree_parameters = settings["refinement_settings"]

    def ExecuteFinalizeSolutionStep(self):
        current_time = self.fluid_model_part.ProcessInfo[KratosMultiphysics.TIME]
        if (current_time <= self.time_average_length):
            KratosExaqute.WeightedDivergenceCalculationProcess(self.fluid_model_part,self.weighted_divergencefree_parameters).Execute()

    def ExecuteFinalize(self):
        # Calculate nodal h
        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.fluid_model_part)
        find_nodal_h.Execute()

        # Calculate divergence metric of the variable
        if (self.domain_size == 2):
            divergencefree_metric = KratosExaqute.MetricDivergenceFreeProcess2D(self.fluid_model_part,self.metric_divergencefree_parameters)
        elif (self.domain_size == 3):
            divergencefree_metric = KratosExaqute.MetricDivergenceFreeProcess3D(self.fluid_model_part,self.metric_divergencefree_parameters)
        divergencefree_metric.Execute()

        # Execute remeshing process
        remesh_parameters = KratosMultiphysics.Parameters()
        if (self.domain_size == 2):
            MmgProcess = KratosMeshing.MmgProcess2D(self.fluid_model_part,remesh_parameters)
        elif (self.domain_size == 3):
            MmgProcess = KratosMeshing.MmgProcess3D(self.fluid_model_part,remesh_parameters)
        MmgProcess.Execute()
