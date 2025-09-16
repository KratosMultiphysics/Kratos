# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

def Create(*args):
    return ComputeResultantsOperation(*args)

class ComputeResultantsOperation(CoSimulationCouplingOperation):
    """This operation computes the resultant forces and moments
    on ModelPart by integrating the reactions of all its nodes.
    """
    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.reference_point = self.settings["reference_point"].GetVector()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)


    def Execute(self):

        # Calculate the total reaction
        model_part = self.interface_data.GetModelPart()
        reaction_force, reaction_moment = KM.ForceAndTorqueUtils.ComputeEquivalentForceAndTorque(
            model_part,
            self.reference_point,
            KM.REACTION,
            KM.REACTION_MOMENT
        )

        # Sign is flipped to go from reaction to action (force)
        force = -1*reaction_force
        moment= -1*reaction_moment
        
        # Save the forces in the model part
        model_part[KMC.RESULTANT_FORCE] = force
        model_part[KMC.RESULTANT_MOMENT] = moment


    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED",
            "reference_point" : [0.0, 0.0, 0.0]
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults



