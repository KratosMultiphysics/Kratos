from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as KratosIga


def Factory(settings, model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return IgaOutputIntegrationDomainProcess(model, settings["Parameters"])

class IgaOutputIntegrationDomainProcess(KratosMultiphysics.Process):

    def __init__(self, model, params):
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""{
            "output_file_name"           : "",
            "model_part_name"            : "",
            "output_geometry_elements"   : true,
            "output_geometry_conditions" : false,
            "output_brep_elements"       : false,
            "output_brep_conditions"     : false
        }""")
        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[self.params["model_part_name"].GetString()]

        self.output_file_name = self.params["output_file_name"].GetString()
        if not self.output_file_name.endswith("_integrationdomain.json"):
            self.output_file_name += "_integrationdomain.json"

        self.output_geometry_elements = self.params["output_geometry_elements"].GetBool()
        self.output_geometry_conditions = self.params["output_geometry_conditions"].GetBool()
        self.output_brep_elements = self.params["output_brep_elements"].GetBool()
        self.output_brep_conditions = self.params["output_brep_conditions"].GetBool()

    def ExecuteBeforeSolutionLoop(self):
        with open(self.output_file_name, 'w') as output_file:
            output_geometry_integration_points =  "\"geometry_integration_points\":["
            output_brep_integration_points =  "\"boundary_integration_points\":[ "
            for element in self.model_part.Elements:
                if not (element.Has(KratosIga.TANGENTS)):
                    if (self.output_geometry_elements):
                        local_coords = element.GetValue(KratosIga.LOCAL_COORDINATES)
                        output_geometry_integration_points += "[" + str(element.Id) + "," + str(element.GetValue(KratosIga.BREP_ID)) + ",[" + str(local_coords[0]) + "," + str(local_coords[1]) + "]],"
                else:
                    if (self.output_brep_elements):
                        local_coords = element.GetValue(KratosIga.LOCAL_COORDINATES)
                        output_brep_integration_points += "[" + str(element.Id) + "," + str(element.GetValue(KratosIga.BREP_ID)) + ",[" + str(local_coords[0]) + "," + str(local_coords[1]) + "]],"

            for condition in self.model_part.Conditions:
                if not (condition.Has(KratosIga.TANGENTS)):
                    if (self.output_geometry_conditions):
                        local_coords = condition.GetValue(KratosIga.LOCAL_COORDINATES)
                        output_geometry_integration_points += "[" + str(condition.Id) + "," + str(condition.GetValue(KratosIga.BREP_ID)) + ",[" + str(local_coords[0]) + "," + str(local_coords[1]) + "]],"
                else:
                    if (self.output_brep_conditions):
                        local_coords = element.GetValue(KratosIga.LOCAL_COORDINATES)
                        output_brep_integration_points += "[" + str(condition.Id) + "," + str(condition.GetValue(KratosIga.BREP_ID)) + ",[" + str(local_coords[0]) + "," + str(local_coords[1]) + "]],"
            #Cut off last ","
            output_geometry_integration_points = output_geometry_integration_points[:-1]
            output_geometry_integration_points += "],"
            output_brep_integration_points = output_brep_integration_points[:-1]
            output_brep_integration_points += "]"

            output_file.write("{" + output_geometry_integration_points + output_brep_integration_points + "}")
