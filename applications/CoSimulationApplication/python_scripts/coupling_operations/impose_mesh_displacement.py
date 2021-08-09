# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# Other imports
import numpy as np

def Create(*args):
    return ComputeResultantsOperation(*args)

class ComputeResultantsOperation(CoSimulationCouplingOperation):
    """This operation computes the Normals (NORMAL) on a given ModelPart
    """
    def __init__(self, settings, solver_wrappers, process_info):
        super().__init__(settings, process_info)
        solver_name = self.settings["solver"].GetString()
        data_name = self.settings["data_name"].GetString()
        self.reference_point = self.settings["reference_point"].GetVector()
        self.interface_data = solver_wrappers[solver_name].GetInterfaceData(data_name)

    def Initialize(self):
        pass

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def InitializeCouplingIteration(self):
        pass

    def FinalizeCouplingIteration(self):
        pass

    def Execute(self):
        
        model_part = self.interface_data.GetModelPart()

        displacement = model_part[KMC.GLOBAL_DISPLACEMENT]
        rotation = model_part[KMC.GLOBAL_ROTATION]
        angle = np.linalg.norm(rotation)

        if angle != 0:
            axis = rotation / angle

            KM.MeshMovingApplication.MoveModelPart(
                model_part,
                axis,                 # rotation axis
                angle,                # rotation angle
                self.reference_point, # one point of the rotation axis
                displacement)         # translation

        '''
        for node in model_part.GetCommunicator().LocalMesh().Nodes:

            dx = displacement[0]
            dy = displacement[1]
            dz = displacement[2]

            if list(rotation) != [0,0,0]:
                
                # TODO: Some things can go out of the loop
                #print()
                rotation = np.array(rotation)
                #print("Rotation vector: "+str(rotation))
                unitary_axis = rotation / np.linalg.norm(rotation)
                #print("Unitary axis: "+str(unitary_axis))
                angle = np.linalg.norm(rotation)
                #print("Angle: "+str(angle))
                node_coord = np.array([node.X0, node.Y0, node.Z0])
                #print("Node coord (original): "+str(node_coord))

                node_coord_rotated = node_coord*np.cos(angle) + np.cross(unitary_axis, node_coord)*np.sin(angle) + unitary_axis*np.dot(unitary_axis,node_coord)*(1-np.cos(angle))
                #print("Node coord (rotated): "+str(node_coord_rotated))



            node.SetSolutionStepValue(KM.MESH_DISPLACEMENT_X, dx)
            node.SetSolutionStepValue(KM.MESH_DISPLACEMENT_Y, dy)
            node.SetSolutionStepValue(KM.MESH_DISPLACEMENT_Z, dx)
        '''

    def PrintInfo(self):
        pass

    def Check(self):
        pass

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"    : "UNSPECIFIED",
            "data_name" : "UNSPECIFIED",
            "reference_point" : [0.0, 0.0, 0.0]
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
