import KratosMultiphysics
import KratosMultiphysics.MeshMovingApplication as MeshMovingApplication


def Factory(params, Model):
    if(type(params) != KratosMultiphysics.Parameters):
        raise Exception(
            'expected input shall be a Parameters object, encapsulating a json string')
    return MeshMovingApplication.ImposeMeshMotionProcess(Model, params["Parameters"])