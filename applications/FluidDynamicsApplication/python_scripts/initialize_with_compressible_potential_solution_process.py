import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
from KratosMultiphysics.CompressiblePotentialFlowApplication import ComputeNodalValueProcess

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return InitializeWithCompressiblePotentialSolutionProcess(model, settings["Parameters"])


class InitializeWithCompressiblePotentialSolutionProcess(KratosMultiphysics.Process):

    """
    Initializes the values by solving a steady-state problem with the compressible potential flow analysis stage.

    Parameters
    ----------
     - model_part_name: The name of the modelpart to clone and run the analysis
            on. Should be the main_model_part of the parent analysis.
     - volume_model_part_name: The modelpart containing the fluid.
     - skin_parts: Inlets and outlets
     - boundary_conditions_process_list: Processes ran for the potential
            simulation.
     - element_type: Type of element. See PotentialFlowSolver._GetStrategyType
            to see the options
     - properties: These are used to convert VELOCITY into MOMENTUM, DENSITY
            and TOTAL_ENERGY

    Furthermore, if any of the listed items is missing in the parameters of a
    FarField process, they will be autofilled based on these properties.
    - free_stream_density
    - sound speed
    - heat capacity ratio
    - mach_infinity

    """

    def __init__(self, model, settings):
        """Stores the model and the validated parameters."""
        super().__init__()
        settings.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model = model
        self.settings = settings
        self.freestream_properties = self._GenerateFreeStreamProperties(self.settings["properties"])


    def ExecuteInitialize(self):
        """
        Automatically sets certain values of FarFieldProcess according to the freestream properties.
        """

        auto_filled_parameters = KratosMultiphysics.Parameters("""
        {
            "free_stream_density" : 0.0,
            "heat_capacity_ratio" : 0.0,
            "speed_of_sound" : 0.0,
            "mach_infinity" : 0.0
        }
        """)

        auto_filled_parameters["free_stream_density"].SetDouble(self.freestream_properties["rho"])
        auto_filled_parameters["heat_capacity_ratio"].SetDouble(self.freestream_properties["gamma"])
        auto_filled_parameters["speed_of_sound"].SetDouble(self.freestream_properties["c"])
        auto_filled_parameters["mach_infinity"].SetDouble(self.freestream_properties["M"])

        for process_parameters in self.settings["boundary_conditions_process_list"]:
            if process_parameters.Has("process_name") and \
                    process_parameters["process_name"].GetString() == "FarFieldProcess":
                process_parameters["Parameters"].AddMissingParameters(auto_filled_parameters)


    def ExecuteBeforeSolutionLoop(self):
        # Creating model part
        potential_mpart = self.model.CreateModelPart(self._GetPotentialModelPartName())
        original_model_part = self.model[self.settings["model_part_name"].GetString()]

        # Starting analysis
        time = original_model_part.ProcessInfo[KratosMultiphysics.TIME]
        analysis_parameters = self._GenerateAnalysisparameters(self.settings, time)
        analysis = PotentialFlowAnalysis(self.model, analysis_parameters)

        # Filling model part
        KratosMultiphysics.MergeVariableListsUtility().Merge(original_model_part, potential_mpart)
        potential_mpart.SetBufferSize(3)

        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        modeler.GenerateModelPart(original_model_part, potential_mpart, "Element2D3N", "LineCondition2D2N")

        # Removing outer solver interfearence
        if potential_mpart.HasSubModelPart("fluid_computational_model_part"):
            potential_mpart.RemoveSubModelPart("fluid_computational_model_part")

        # Running analysis
        analysis.Run()

        # Calculate the velocity and pressure nodal projections
        nodal_variables = ["VELOCITY","PRESSURE"]
        nodal_value_process = ComputeNodalValueProcess(potential_mpart, nodal_variables)
        nodal_value_process.Execute()

        for node in potential_mpart.Nodes:
            local_props = self._ComputeLocalProperties(node)

            # Setting initial conditions
            for variable,value in local_props.items():
                node.SetSolutionStepValue(variable, value)

        potential_mpart.CloneSolutionStep()


    def Clear(self):
        self.model.DeleteModelPart(self._GetPotentialModelPartName())


    def GetDefaultParameters(self):
        return KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "",
            "volume_model_part_name" : "",
            "skin_parts" : [],
            "boundary_conditions_process_list" : [],
            "element_type" : "compressible",
            "properties" : {
                "c_v" : 722.14,
                "gamma" : 1.4,
                "free_stream_density" : 1.19659,
                "free_stream_momentum" : 0,
                "free_stream_energy" : 2.53313e5
            }
        }
        """)


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
        energy = rho * (self.freestream_properties["c_v"] * self.freestream_properties["temperature"] + 0.5 * vel2)

        local_properties = {}
        local_properties[KratosMultiphysics.VELOCITY] = vel
        local_properties[KratosMultiphysics.PRESSURE] = node.GetValue(KratosMultiphysics.PRESSURE)
        local_properties[KratosMultiphysics.DENSITY] = rho
        local_properties[KratosMultiphysics.MOMENTUM] = vel*rho
        local_properties[KratosMultiphysics.TOTAL_ENERGY] = energy

        return local_properties


    @classmethod
    def _GetPotentialModelPartName(cls):
        """Giving it a unique name to avoid colisions"""
        return "initial_potential_process_model_part"


    @classmethod
    def _GenerateAnalysisparameters(cls, settings, simulation_time):
        defaults = KratosMultiphysics.Parameters("""
        {
            "problem_data"     : {
                "problem_name"  : "initialize_with_compressible_potential_solution_process_internal_analyisis",
                "parallel_type" : "OpenMP",
                "echo_level"    : 0,
                "start_time"    : 0.0,
                "end_time"      : 1.0
            },
            "output_processes" : {
            },
            "solver_settings"  : {
                "model_part_name"          : "{replaceme}",
                "domain_size"              : 2,
                "solver_type"              : "potential_flow",
                "model_import_settings"    : {
                    "input_type"     : "use_input_model_part"
                },
                "formulation": {
                    "element_type" : "{replaceme}"
                },
                "maximum_iterations"       : 10,
                "echo_level"               : 0,
                "volume_model_part_name"   : "{replaceme}",
                "skin_parts"               : ["{replaceme}"],
                "no_skin_parts"            : [],
                "reform_dofs_at_each_step" : false
            },
            "processes"        : {
                "boundary_conditions_process_list" : [],
                "auxiliar_process_list"            : []
            }
        }
        """)

        defaults["solver_settings"]["model_part_name"].SetString(cls._GetPotentialModelPartName())
        defaults["solver_settings"]["volume_model_part_name"] = settings["volume_model_part_name"]
        defaults["solver_settings"]["formulation"]["element_type"] = settings["element_type"]
        defaults["solver_settings"]["skin_parts"] = settings["skin_parts"]
        defaults["processes"]["boundary_conditions_process_list"] = settings["boundary_conditions_process_list"]

        # Ensuring time ends at the begging of simulation
        defaults["problem_data"]["end_time"].SetDouble(simulation_time)
        defaults["problem_data"]["start_time"].SetDouble(simulation_time - 1.0)

        return defaults
