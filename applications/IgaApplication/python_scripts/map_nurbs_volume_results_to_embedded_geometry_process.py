# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IGA
from KratosMultiphysics import kratos_utilities

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MapNurbsVolumeResultsToEmbeddedGeometryProcess(model, settings["Parameters"])

class MapNurbsVolumeResultsToEmbeddedGeometryProcess(KratosMultiphysics.Process):
    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)
        self.params = params
        self.process = IGA.MapNurbsVolumeResultsToEmbeddedGeometryProcess(model, params)

    def ExecuteBeforeOutputStep(self):
        nodal_variable_list = kratos_utilities.GenerateVariableListFromInput(self.params["nodal_results"])

        for nodal_variable in nodal_variable_list:
            if( nodal_variable == KratosMultiphysics.DISPLACEMENT):
                self.process.MapNodalValues(nodal_variable)
            else:
                warn_msg  = "Mapping of Variable '" + nodal_variable.Name() + "' is not available."
                KratosMultiphysics.Logger.PrintWarning("::[MapNurbsVolumeResultsToEmbeddedGeometryProcess]::", warn_msg)

