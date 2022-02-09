import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import KratosMultiphysics.MeshingApplication as KratosMeshing


def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyPotentialFlowHessianRemeshingProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyPotentialFlowHessianRemeshingProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"   :"",
                "metric_parameters" : {
                    "minimal_size"                        : 0.001,
                    "maximal_size"                        : 10.0,
                    "enforce_current"                     : false,
                    "hessian_strategy_parameters":
                    {
                        "non_historical_metric_variable"  : true,
                        "estimate_interpolation_error"    : false,
                        "interpolation_error"             : 1e-2
                    },
                    "anisotropy_remeshing"                : false
                },
                "mmg_parameters"  :{
                    "discretization_type"              : "STANDARD",
                    "save_external_files"              : false,
                    "initialize_entities"              : false,
                    "echo_level"                       : 0
                }
            }  """ )
        settings.ValidateAndAssignDefaults(default_parameters)

        self.main_model_part = Model[settings["model_part_name"].GetString()]
        self.domain_size = self.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)
        if not self.domain_size == 2 and not self.domain_size == 3:
            raise(Exception("Domain size is different than 2 and different than 3. Please set the domain size correclty"))

        self.metric_parameters = settings["metric_parameters"]
        self.mmg_parameters = settings["mmg_parameters"]

        if self.metric_parameters["hessian_strategy_parameters"].Has("non_historical_metric_variable"):
            if not self.metric_parameters["hessian_strategy_parameters"]["non_historical_metric_variable"].GetBool():
               raise(Exception("Potential Flow remeshing process uses non historical velocity variable!"))
        else:
            self.metric_parameters["hessian_strategy_parameters"].AddEmptyValue("non_historical_metric_variable")
            self.metric_parameters["hessian_strategy_parameters"]["non_historical_metric_variable"].SetBool(True)

    def ExecuteFinalize(self):

        self.__ComputeNodalVelocity()

        self.__ComputeHessianMetric()

        self.__RemoveSubModelParts()

        self.__ExecuteRefinement()

    def __ComputeNodalVelocity(self):

        nodal_velocity_process = CPFApp.ComputeNodalValueProcess(self.main_model_part, ["VELOCITY"])
        nodal_velocity_process.Execute()

    def __ComputeHessianMetric(self):

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.main_model_part)
        find_nodal_h.Execute()

        metric_x = KratosMeshing.ComputeHessianSolMetricProcess(self.main_model_part, KratosMultiphysics.VELOCITY_X, self.metric_parameters)
        metric_x.Execute()
        metric_y = KratosMeshing.ComputeHessianSolMetricProcess(self.main_model_part, KratosMultiphysics.VELOCITY_Y, self.metric_parameters)
        metric_y.Execute()

        if (self.domain_size == 3):
            metric_z = KratosMeshing.ComputeHessianSolMetricProcess(self.main_model_part, KratosMultiphysics.VELOCITY_Z, self.metric_parameters)
            metric_z.Execute()

    def __RemoveSubModelParts(self):

        self.main_model_part.RemoveSubModelPart('wake_sub_model_part')
        self.main_model_part.RemoveSubModelPart('trailing_edge_sub_model_part')
        self.main_model_part.RemoveSubModelPart('fluid_computational_model_part')

    def __ExecuteRefinement(self):

        if (self.domain_size == 2):
            MmgProcess = KratosMeshing.MmgProcess2D(self.main_model_part, self.mmg_parameters)
            MmgProcess.Execute()
        elif (self.domain_size == 3):
            MmgProcess = KratosMeshing.MmgProcess3D(self.main_model_part, self.mmg_parameters)
            MmgProcess.Execute()



