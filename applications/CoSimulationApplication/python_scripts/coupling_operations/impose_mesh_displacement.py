# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC
import KratosMultiphysics.MeshMovingApplication as KMM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# Other imports
import numpy as np

def Create(*args):
    return ImposeMeshDisplacementOperation(*args)

class ImposeMeshDisplacementOperation(CoSimulationCouplingOperation):
    """This operation takes the global displacements of a model part
    and applies the corresponding displacements to each node in the mesh.
    """
    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.reference_point = self.settings["reference_point"].GetVector()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)
        

    def Execute(self):
        
        # Get global displacements
        model_part = self.interface_data.GetModelPart()
        displacement = model_part[KMC.GLOBAL_DISPLACEMENT]
        rotation = model_part[KMC.GLOBAL_ROTATION]
        angle = np.linalg.norm(rotation)

        # If the angle is 0, the axis doesn't really matter
        if angle != 0:
            axis = rotation / angle
        else:
            axis = [1,0,0]

        # Apply displacement to mesh
        KMM.MoveModelPart(
            model_part,
            axis,                 # rotation axis
            angle,                # rotation angle
            self.reference_point, # one point of the rotation axis
            displacement)         # translation
        

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED",
            "reference_point" : [0.0, 0.0, 0.0]
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
