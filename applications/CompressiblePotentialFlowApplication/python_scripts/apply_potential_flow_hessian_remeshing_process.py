import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import KratosMultiphysics.MeshingApplication as KratosMeshing
from KratosMultiphysics.gid_output_process import GiDOutputProcess
import time

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
                "metric_variables_list"               : ["VELOCITY"],
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
                    "preserve_flags"              : false,
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
        self.metric_variables_list = settings["metric_variables_list"].GetStringArray()

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

        nodal_velocity_process = CPFApp.ComputeNodalValueProcess(self.main_model_part, ["VELOCITY","PRESSURE_COEFFICIENT"])
        nodal_velocity_process.Execute()

    def __ComputeHessianMetric(self):

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.main_model_part)
        find_nodal_h.Execute()

        # KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D,self.main_model_part.Nodes)

        for node in self.main_model_part.Nodes:
            final_size = node.GetValue(KratosMultiphysics.NODAL_H)
            metric_tensor_2d = KratosMultiphysics.Vector(3, 0.0)
            metric_tensor_2d[0] =1/final_size/final_size
            metric_tensor_2d[1] =1/final_size/final_size
            node.SetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D, metric_tensor_2d)

        if "VELOCITY_POTENTIAL" in self.metric_variables_list:
            print("COMPUTING VELOCITY_POTENTIAL METRIC")
            copy_parameters = self.metric_parameters.Clone()
            metric_x = KratosMultiphysics.CompressiblePotentialFlowApplication.ComputePotentialHessianSolMetricProcess(self.main_model_part, copy_parameters)
            metric_x.Execute()

        if "PRESSURE_COEFFICIENT" in self.metric_variables_list:
            print("COMPUTING PRESSURE_COEFFICIENT METRIC")
            copy_parameters = self.metric_parameters.Clone()
            metric_x = KratosMeshing.ComputeHessianSolMetricProcess(self.main_model_part, KratosMultiphysics.PRESSURE_COEFFICIENT, copy_parameters)
            metric_x.Execute()

        if "VELOCITY" in self.metric_variables_list:
            print("COMPUTING VELOCITY METRIC")
            copy_parameters = self.metric_parameters.Clone()
            metric_x = KratosMeshing.ComputeHessianSolMetricProcess(self.main_model_part, KratosMultiphysics.VELOCITY_X, copy_parameters)
            metric_x.Execute()
            metric_y = KratosMeshing.ComputeHessianSolMetricProcess(self.main_model_part, KratosMultiphysics.VELOCITY_Y, copy_parameters)
            metric_y.Execute()

        if "ADJOINT_VELOCITY" in self.metric_variables_list:
            print("COMPUTING ADJOINT_VELOCITY METRIC")
            copy_parameters = self.metric_parameters.Clone()
            metric_x = KratosMeshing.ComputeHessianSolMetricProcess(self.main_model_part, CPFApp.ADJOINT_VELOCITY_X, copy_parameters)
            metric_x.Execute()
            metric_y = KratosMeshing.ComputeHessianSolMetricProcess(self.main_model_part, CPFApp.ADJOINT_VELOCITY_Y, copy_parameters)
            metric_y.Execute()

        if (self.domain_size == 3):
            metric_z = KratosMeshing.ComputeHessianSolMetricProcess(self.main_model_part, KratosMultiphysics.VELOCITY_Z, self.metric_parameters)
            metric_z.Execute()

    def __RemoveSubModelParts(self):

        self.main_model_part.RemoveSubModelPart('wake_sub_model_part')
        self.main_model_part.RemoveSubModelPart('trailing_edge_sub_model_part')
        self.main_model_part.RemoveSubModelPart('fluid_computational_model_part')
        self.main_model_part.RemoveSubModelPart('lower_surface_sub_model_part')
        self.main_model_part.RemoveSubModelPart('upper_surface_sub_model_part')

    def __ExecuteRefinement(self):
        ini_time = time.time()
        if (self.domain_size == 2):
            MmgProcess = KratosMeshing.MmgProcess2D(self.main_model_part, self.mmg_parameters)
            MmgProcess.Execute()
        elif (self.domain_size == 3):
            MmgProcess = KratosMeshing.MmgProcess3D(self.main_model_part, self.mmg_parameters)
            MmgProcess.Execute()
        print("REMESHED NUMBER OF NODES", self.main_model_part.NumberOfNodes())
        print("REMESHED NUMBER OF ELEMENTS", self.main_model_part.NumberOfElements())
        print("REMESHED TIME", time.time() - ini_time)



