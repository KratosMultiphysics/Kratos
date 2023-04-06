import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.CompressiblePotentialFlowApplication as KratosCompressible

class PotentialToCompressibleNavierStokesOperation(KratosMultiphysics.Operation):

    def __init__(self, model, settings):
        super().__init__()

        # Cannot validate as settings may differ among input types
        settings.AddMissingParameters(self.__GetDefaultParameters())

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
        reference_temperature = self.settings["reference_temperature"].GetDouble()
        gamma = origin_model_part.ProcessInfo.GetValue(KratosFluid.HEAT_CAPACITY_RATIO)
        sound_velocity = origin_model_part.ProcessInfo.GetValue(KratosMultiphysics.SOUND_VELOCITY)
        free_stream_rho = origin_model_part.ProcessInfo.GetValue(KratosCompressible.FREE_STREAM_DENSITY)
        free_stream_mach = origin_model_part.ProcessInfo.GetValue(KratosCompressible.FREE_STREAM_MACH)
        specific_heat = (sound_velocity**2 / (gamma * reference_temperature)) / (gamma - 1.0)

        if origin_model_part.NumberOfNodes() != destination_model_part.NumberOfNodes():
            raise RuntimeError("Origin and destination modelparts have a different number of nodes.")

        for potential_node, ns_node in zip(origin_model_part.Nodes, destination_model_part.Nodes):

            velocity = potential_node.GetValue(KratosMultiphysics.VELOCITY)
            velocity_norm_2 = velocity[0] * velocity[0] + velocity[1] * velocity[1]
            velocity_norm = velocity_norm_2**0.5
            mach = velocity_norm / sound_velocity
            num = 1.0 + 0.5 * (gamma - 1.0) * free_stream_mach**2
            den = 1.0 + 0.5 * (gamma - 1.0) * mach**2
            density = free_stream_rho * (num / den)**(1.0 / (gamma - 1.0))
            energy = specific_heat * reference_temperature + 0.5 * velocity_norm_2

            # Density initial condition
            ns_node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0, density)
            # Momentum initial condition
            ns_node.SetSolutionStepValue(KratosMultiphysics.MOMENTUM, 0, density * velocity)
            # Total energy initial condition
            ns_node.SetSolutionStepValue(KratosMultiphysics.TOTAL_ENERGY, 0, density * energy) 
            
    @classmethod
    def __GetDefaultParameters(self):
        return KratosMultiphysics.Parameters('''{
                        "origin_model_part"      : "",
                        "destination_model_part" : "",
                        "variable_names_to_pass" : [""],
                        "reference_temperature"  : 273,
                        "echo_level"             : 0
        }''')
    
    @staticmethod
    def Create(settings, model):
        return PotentialToCompressibleNavierStokesOperation(settings, model)
