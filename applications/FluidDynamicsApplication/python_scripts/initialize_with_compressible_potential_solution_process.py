import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
from KratosMultiphysics.CompressiblePotentialFlowApplication import ComputeNodalValueProcess

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return InitializeWithCompressiblePotentialSolutionProcess(model, settings["Parameters"])


class InitializeWithCompressiblePotentialSolutionProcess(KratosMultiphysics.Process):
    """Initializes the values by solving a steady-state problem with the compressible potential flow analysis stage"""
    def __init__(self, model, settings):
        super().__init__()
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model = model
        self.settings = settings

        self.analysis_parameters = self._GenerateAnalysisparameters(self.settings)
        self.analysis = None

        self.freestream_properties = self._GenerateFreeStreamProperties(settings["properties"])


    def GetDefaultParameters(self):
        return KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "",
            "volume_model_part_name" : "",
            "skin_parts" : [],
            "boundary_conditions_process_list" : [],
            "properties" : {
                "c_v" : 722.14,
                "gamma" : 1.4,
                "free_stream_density" : 1.19659,
                "free_stream_momentum" : 0,
                "free_stream_energy" : 2.53313e5
            }
        }
        """)

    def ExecuteInitialize(self):
        # Creating model part
        potential_mpart = self.model.CreateModelPart("potential_analysis_model_part")

        # Starting analysis
        self.analysis = PotentialFlowAnalysis(self.model, self.analysis_parameters)

        # Filling model part
        original_model_part = self.model[self.settings["model_part_name"].GetString()]
        KratosMultiphysics.MergeVariableListsUtility().Merge(original_model_part, potential_mpart)

        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        modeler.GenerateModelPart(original_model_part, potential_mpart, "Element2D3N", "LineCondition2D2N")

        # potential_mpart.ProcessInfo = pre_merger_pinfo

        # Running analysis
        self.analysis.Run()

        # Calculate the velocity and pressure nodal projections
        computing_model_part = self.analysis._GetSolver().GetComputingModelPart()
        nodal_variables = ["VELOCITY","PRESSURE"]
        nodal_value_process = ComputeNodalValueProcess(computing_model_part, nodal_variables)
        nodal_value_process.Execute()

        for node in potential_mpart.Nodes:
            local_props = self._ComputeLocalProperties(node)

            # Density initial condition
            for variable,value in local_props.items():
                node.SetSolutionStepValue(variable, 0, value)
                node.SetSolutionStepValue(variable, 1, value)
            

    @classmethod
    def _GenerateFreeStreamProperties(cls, settings):
        c_v = settings["c_v"].GetDouble()
        gamma = settings["gamma"].GetDouble()
        rho = settings["free_stream_density"].GetDouble()
        momentum = settings["free_stream_momentum"].GetDouble()
        total_energy = settings["free_stream_energy"].GetDouble()

        vel = momentum / rho
        temperature = (total_energy / rho - 0.5 * vel**2) / c_v
        sound_speed = (gamma*(gamma-1.0)*c_v*temperature)**0.5
        mach_number = vel / sound_speed

        freestream_properties = {}
        freestream_properties["c_v"] = c_v
        freestream_properties["gamma"] = gamma
        freestream_properties["rho"] = rho
        freestream_properties["temperature"] = temperature
        freestream_properties["c"] = sound_speed
        freestream_properties["M"] = mach_number

        return freestream_properties

    def _ComputeLocalProperties(self, node):
        vel = node.GetValue(KratosMultiphysics.VELOCITY)
        vel2 = vel[0] * vel[0] + vel[1] * vel[1]
        vel_norm = vel2**0.5
        mach = vel_norm / self.freestream_properties["c"]
        num = 1.0 + 0.5 * (self.freestream_properties["gamma"] - 1.0) * self.freestream_properties["M"]**2
        det = 1.0 + 0.5 * (self.freestream_properties["gamma"] - 1.0) * mach**2
        rho = self.freestream_properties["rho"] * (num / det)**(1.0 / (self.freestream_properties["gamma"] - 1.0))
        total_energy = rho * (self.freestream_properties["c_v"] * self.freestream_properties["temperature"] + 0.5 * vel2)

        local_properties = {}
        local_properties[KratosMultiphysics.DENSITY] = rho
        local_properties[KratosMultiphysics.VELOCITY] = vel
        local_properties[KratosMultiphysics.MOMENTUM] = vel*rho
        local_properties[KratosMultiphysics.TOTAL_ENERGY] = total_energy
        local_properties[KratosMultiphysics.PRESSURE] = vel = node.GetValue(KratosMultiphysics.PRESSURE)

        return local_properties

    @classmethod
    def _GenerateAnalysisparameters(cls, settings):
        defaults = KratosMultiphysics.Parameters("""
        {
            "problem_data"     : {
                "problem_name"  : "potential_analysis_model_part",
                "parallel_type" : "OpenMP",
                "echo_level"    : 0,
                "start_time"    : 0.0,
                "end_time"      : 1
            },
            "output_processes" : {
            },
            "solver_settings"  : {
                "model_part_name"          : "potential_analysis_model_part",
                "domain_size"              : 2,
                "solver_type"              : "potential_flow",
                "model_import_settings"    : {
                    "input_type"     : "use_input_model_part"
                },
                "formulation": {
                    "element_type" : "compressible"
                },
                "maximum_iterations"       : 10,
                "echo_level"               : 0,
                "volume_model_part_name"   : "{replaceme}",
                "skin_parts"               : [],
                "no_skin_parts"            : [],
                "reform_dofs_at_each_step" : false
            },
            "processes"        : {
                "boundary_conditions_process_list" : [],
                "auxiliar_process_list"            : []
            }
        }
        """)

        defaults["solver_settings"]["volume_model_part_name"] = settings["volume_model_part_name"]
        defaults["solver_settings"]["skin_parts"] = settings["skin_parts"]
        defaults["processes"]["boundary_conditions_process_list"] = settings["boundary_conditions_process_list"]

        return defaults