from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.NurbsBrepApplication import *

KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, model_part):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return NurbsBrepProcess(model_part, settings)

class NurbsBrepProcess(KratosMultiphysics.Process):

    def __init__(self,model_part,params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "process_name": "NurbsBrepProcess",
            "cad_geometry_file_name": "geometry.json",
            "integration_domain_model_part_name": "IntegrationDomain",
            "parameters": {
                "integration_domain_parameter": {
                    "integration_scheme": "triangulation",
                    "accuracy": 1E-07,
                    "shape_function_derivatives_order": 2,
                    "number_projection_iterations": 20,
                    "polygon_discretization": 10,
                    "integration_domains": {
                    "faces": true,
                    "faces_embedded": false,
                    "faces_reversed": false,
                    "edges": true,
                    "coupling_edges": true
                  }
              },
              "geometry_refinement": [
                {
                  "selection": "ALL",
                  "parameters": {
                    "knot_insertions_u": [ ],
                    "knot_insertions_v": [ ],
                    "multiply_knots_u": 0,
                    "multiply_knots_v": 0,
                    "max_element_size_u": 0,
                    "max_element_size_v": 0,
                    "min_order_p": 0,
                    "min_order_q": 0,
                    "order_elevation_p": 0,
                    "order_elevation_q": 0
                  }
                }
              ]
            },
            "compute_surface_area": true,
            "write_integration_domain_global": true,
            "write_integration_domain_json_local": true,
            "integration_domain_file_name": "problem_name_integrationdomain.json"
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model_part

    def ExecuteInitialize(self):
        self.model_part_integration_domain = KratosMultiphysics.ModelPart(self.params["integration_domain_model_part_name"].GetString())

        with open(self.params["cad_geometry_file_name"].GetString(),'r') as cad_geometry_file:
            cad_geometry = KratosMultiphysics.Parameters( cad_geometry_file.read())
        self.this_modeler = KratosMultiphysics.NurbsBrepApplication.NurbsBrepModeler(self.model_part)
        self.geometry_reader = KratosMultiphysics.NurbsBrepApplication.BrepModelGeometryReader(cad_geometry)
        self.this_modeler.LoadGeometry(self.geometry_reader)

        self.this_modeler.ApplyGeometryRefinement(self.params["parameters"]["geometry_refinement"])

        self.this_modeler.CreateIntegrationDomain(self.params["parameters"]["integration_domain_parameter"], self.model_part_integration_domain)

        if(self.params["compute_surface_area"].GetBool()):
            self.this_modeler.ComputeArea(self.model_part_integration_domain)

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        if(self.params["write_integration_domain_global"].GetBool()):
            self.geometry_reader.WriteGaussPoints(self.model_part_integration_domain)
        if (self.params["write_integration_domain_json_local"].GetBool()):
            self.geometry_reader.WriteGaussPointsJson(self.model_part_integration_domain, self.params["integration_domain_file_name"].GetString())



