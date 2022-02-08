# Import Python modules
import json
import numpy

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.RomApplication as KratosROM

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string.")
    return ApplyHRomVisualizationMeshProjectionProcess(model, settings["Parameters"])

class ApplyHRomVisualizationMeshProjectionProcess(KratosMultiphysics.Process):
    """Auxiliary process to interface the cpp implementation of HRomVisualizationMeshProjectionProcess."""

    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        self.hrom_visualization_mesh_projection_process = KratosROM.HRomVisualizationMeshProjectionProcess(model, settings)

    def ExecuteInitialize(self):
        self.hrom_visualization_mesh_projection_process.ExecuteInitialize()

    def ExecuteBeforeOutputStep(self):
        self.hrom_visualization_mesh_projection_process.ExecuteBeforeOutputStep()