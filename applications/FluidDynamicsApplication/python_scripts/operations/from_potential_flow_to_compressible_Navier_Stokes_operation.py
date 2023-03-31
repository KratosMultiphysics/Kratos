import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as KratosCompressible

class FromPotentialFlowToCompressibleNavierStokesOperation(KratosMultiphysics.Operation):

    def __init__(self, model, settings):
        super().__init__()

        # Cannot validate as settings may differ among input types
        settings.AddMissingParameters(self.__GetDefaultSettings())

        # Declare required member variables
        self.model = model
        self.settings = settings

    def Execute(self):

        origin_model_part_name = self.settings["origin_model_part"].GetString()
        destination_model_part_name = self.settings["destination_model_part"].GetString()
        variable_names = self.settings["variable_names_to_pass"].GetStringArray()

        # saving the modelparts
        origin_model_part = self.model.GetModelPart(origin_model_part_name)
        destination_model_part = self.model.GetModelPart(destination_model_part_name)

        # computing nodal values
        nodal_value_process = KratosCompressible.ComputeNodalValueProcess(origin_model_part, variable_names)
        nodal_value_process.Execute()

        # Transfer the potential flow values as initial condition for the compressible problem
        tem_0 = 273
        gamma = 1.4
        c_v = 722.14
        c = (gamma*(gamma-1.0)*c_v*tem_0)**0.5
        free_stream_rho = 1.0
        free_stream_mach = 0.8

        if len(origin_model_part.Nodes) != len(destination_model_part.Nodes):
            raise RuntimeError("CopyValueOperation: Model parts have a different number of nodes!")

        for potential_node, ns_node in zip(origin_model_part.Nodes, destination_model_part.Nodes):

            pot_v = potential_node.GetValue(KratosMultiphysics.VELOCITY)

            pot_v_norm_2 = pot_v[0] * pot_v[0] + pot_v[1] * pot_v[1]
            pot_v_norm = pot_v_norm_2**0.5
            mach = pot_v_norm / c
            num = 1.0 + 0.5 * (gamma - 1.0) * free_stream_mach**2
            den = 1.0 + 0.5 * (gamma - 1.0) * mach**2
            pot_rho = free_stream_rho * (num / den)**(1.0 / (gamma - 1.0))
            aux_tot_ener = pot_rho * (c_v * tem_0 + 0.5 * pot_v_norm_2)

            # Density initial condition
            ns_node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0, pot_rho)
            # Momentum initial condition
            ns_node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM, 0, pot_rho * pot_v)
            # Total energy initial condition
            ns_node.SetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY, 0, aux_tot_ener) 
            

    def __GetDefaultSettings(self):
        default_settings = KratosMultiphysics.Parameters('''{
                        "origin_model_part"      : "",
                        "destination_model_part" : "",
                        "variable_names_to_pass" : [""],
                        "echo_level"             : 0
        }''')
        return default_settings
    
    @staticmethod
    def Create(settings, model):
        return FromPotentialFlowToCompressibleNavierStokesOperation(settings, model)
