import KratosMultiphysics
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
from KratosMultiphysics.CompressiblePotentialFlowApplication import ComputeNodalValueProcess
from bin.RelWithDebInfo.KratosMultiphysics import analysis_stage

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return InitializeWithCompressiblePotentialSolutionProcess(model, settings["Parameters"])


class InitializeWithCompressiblePotentialSolutionProcess(KratosMultiphysics.Process):
    """Initializes the values by solving a steady-state problem with the compressible potential flow analysis stage"""
    def __init__(self, model, settings):
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.model = model
        self.settings = settings

        self.analysis_parameters = self._GenerateAnalysisparameters(self.settings)
        self.analysis = None


    def GetDefaultParameters(self):
        return KratosMultiphysics.Parameters("""
        {
            "model_part_name" : ""
            "volume_model_part_name" : "",
            "skin_parts" : [],
            "materials_filename" : "FluidMaterials.json",
            "boundary_conditions_process_list" : [],
            "free_stream_density" : 1.19659,
            "free_stream_momentum" : 0,
            "free_stream_energy" : 2.53313e5
        }
        """)

    def ExecuteInitialize(self):
        self._GeneratePotentialModelPart()

        self.analysis = PotentialFlowAnalysis(self.model, self.analysis_parameters)
        self.analysis.Run()

        # Calculate the velocity and pressure nodal projections
        computing_model_part = self.analysis._GetSolver().GetComputingModelPart()
        nodal_variables = ["VELOCITY","PRESSURE"]
        nodal_value_process = ComputeNodalValueProcess(computing_model_part, nodal_variables)
        nodal_value_process.Execute()

        # Transfer the potential flow values as initial condition for the compressible problem
        c_v = computing_model_part.Properties.GetValue("SPECIFIC_HEAT")
        gamma = computing_model_part.Properties.GetValue("HEAT_CAPACITY_RATIO")

        inlet_properties = self._ComputeInletProperties(c_v, gamma)

        for node in computing_model_part.Nodes():
            rho, momentum, total_energy = self._ComputeLocalProperties(node, gamma, c_v, inlet_properties)

            # Density initial condition
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0, rho)
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 1, rho)
            # Momentum initial condition
            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM, 0, momentum)
            node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM, 1, momentum)
            # Total energy initial condition
            node.SetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY, 0, total_energy)
            node.SetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY, 1, total_energy)

    def _ComputeInletProperties(self, c_v, gamma):
        inlet_rho = self.settings["free_stream_rho"].GetDouble()
        inlet_momentum = self.settings["free_stream_momentum"].GetDouble()
        inlet_energy = self.settings["free_stream_energy"].GetDouble()

        inlet_v = inlet_momentum / inlet_rho
        inlet_temp = (inlet_energy / inlet_rho - 0.5 * inlet_v**2) / c_v
        inlet_c = (gamma*(gamma-1.0)*c_v*inlet_temp)**0.5
        inlet_mach = inlet_v / inlet_c

        inlet_properties = {}
        inlet_properties["rho"] = inlet_rho
        inlet_properties["temperature"] = inlet_temp
        inlet_properties["c"] = inlet_c
        inlet_properties["M"] = inlet_mach

        return inlet_properties

    def _ComputeLocalProperties(self, node, gamma, c_v, inlet_properties):
        vel = node.GetValue(KratosMultiphysics.VELOCITY)
        vel2 = vel[0] * vel[0] + vel[1] * vel[1]
        vel_norm = vel2**0.5
        mach = vel_norm / inlet_properties["c"]
        num = 1.0 + 0.5 * (gamma - 1.0) * inlet_properties["M"]**2
        det = 1.0 + 0.5 * (gamma - 1.0) * mach**2
        rho = inlet_properties["rho"] * (num / det)**(1.0 / (gamma - 1.0))
        total_energy = rho * (c_v * inlet_properties["temperature"] + 0.5 * vel2)

        return rho, vel_norm*rho, total_energy

    def _GeneratePotentialModelPart(self):
        original_model_part = self.model[self.settings["model_part_name"].GetString()]
        new_model_part = self.model.CreateModelPart("model_part_potential_analysis")

        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        modeler.GenerateModelPart(original_model_part, new_model_part, "Element2D3N", "LineCondition2D2N")

    @classmethod
    def _GenerateAnalysisparameters(cls, settings):
        defaults = KratosMultiphysics.Parameters("""
        {
            "problem_data"     : {
                "problem_name"  : "internal_potential_solver",
                "parallel_type" : "OpenMP",
                "echo_level"    : 0,
                "start_time"    : 0.0,
                "end_time"      : 1
            },
            "output_processes" : {
            },
            "solver_settings"  : {
                "model_part_name"          : "model_part_potential_analysis",
                "domain_size"              : 2,
                "solver_type"              : "potential_flow",
                "model_import_settings"    : {
                    "input_type"     : "use_input_model_part"
                },
                "material_import_settings" : {
                    "materials_filename" : "{replace_me}"
                },
                "formulation": {
                    "element_type" : "compressible"
                },
                "maximum_iterations"       : 10,
                "echo_level"               : 0,
                "volume_model_part_name"   : "{replace_me}",
                "skin_parts"               : [{replace_me}],
                "no_skin_parts"            : [],
                "reform_dofs_at_each_step" : false
            },
            "processes"        : {
                "boundary_conditions_process_list" : [{replace_me}],
                "auxiliar_process_list"            : []
            }
        }
        """)
        defaults["solver_settings"]["model_part_name"] = settings["model_part_name"]
        defaults["solver_settings"]["volume_model_part_name"] = settings["volume_model_part_name"]
        defaults["solver_settings"]["skin_parts"] = settings["skin_parts"]
        defaults["solver_settings"]["materials_filename"] = settings["materials_filename"]
        defaults["boundary_conditions_process_list"] = settings["boundary_conditions_process_list"]

        return defaults
