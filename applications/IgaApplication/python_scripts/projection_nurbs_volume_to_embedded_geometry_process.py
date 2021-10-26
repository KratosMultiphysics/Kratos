# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IGA
from KratosMultiphysics import kratos_utilities

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ProjectionNurbsVolumeToEmbeddedGeometryProcess(model, settings["Parameters"])

class ProjectionNurbsVolumeToEmbeddedGeometryProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)
        self.params = params
        self.process = IGA.ProjectionNurbsVolumeToEmbeddedGeometryProcess(model, params)

    def ExecuteBeforeOutputStep(self):
        nodal_variable_list = kratos_utilities.GenerateVariableListFromInput(self.params["nodal_results"])
        #gauss_point_variables = kratos_utilities.GenerateVariableListFromInput(self.params["gauss_point_results"])
        print("0")
        for nodal_variable in nodal_variable_list:
            self.process.MapNodalValues(nodal_variable)
        print("1")