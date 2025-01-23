# Kratos imports
import KratosMultiphysics as KM
import KratosMultiphysics.ConvectionDiffusionApplication as CDA

# External imports 
from shapely.geometry import Polygon
import numpy as np

def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DefineActiveElementsProcess(Model, settings["Parameters"])


class DefineActiveElementsProcess(KM.Process):
    def __init__(self, model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
            {
                "background_model_part_name"                 : "please_specify_model_part_name",
                "embedded_body_model_part_name"              : "please_specify_model_part_name",
                "keep_external_domain"                       : false
            }
            """
            )
        
        self.settings = settings

        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model = model

        self.keep_external_domain = self.settings["keep_external_domain"].GetBool()

        # Retrieve from the input file the model part names 
        self.background_model_part_name = self.settings["background_model_part_name"].GetString()
        self.embedded_body_model_part_name = self.settings["embedded_body_model_part_name"].GetString()
        
        # Get the model parts 
        self.background_model_part = self.model.GetModelPart(self.background_model_part_name)
        self.embedded_body_model_part = self.model.GetModelPart(self.embedded_body_model_part_name)


    def Check(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        # Create sub model parts containing the elements intersected by the boundary and the active elements (to be integrated)
        intersected_elements_sub_model_part = self.background_model_part.CreateSubModelPart("intersected_elements")
        active_elements_sub_model_part = self.background_model_part.CreateSubModelPart("active_elements")
        inactive_elements_sub_model_part = self.background_model_part.CreateSubModelPart("inactive_elements")

        # Create a sub model part containing the discretization of the embedded body or boundary
        embedded_body_boundary_model_part = self.model.CreateModelPart("embedded_body_boundary")
        for condition in self.embedded_body_model_part.Conditions:
            embedded_body_boundary_model_part.AddCondition(condition)

        # Construct a polygon with the embedded body representation
        points_x_coordinate = []
        points_y_coordinate = []
        for condition in self.embedded_body_model_part.Conditions:
            condition_nodes = condition.GetNodes()
            points_x_coordinate.append(condition_nodes[1].X)
            points_y_coordinate.append(condition_nodes[1].Y)
        points = list(zip(points_x_coordinate, points_y_coordinate))
        boundary_polygon = Polygon(points)


        # Iterate over the elements in the background mesh and define if they are active or no
        for element in self.background_model_part.Elements:
            nodes_x_coordinate = []
            nodes_y_coordinate = []

            # Recover the nodes coordinates for an element 
            for node in element.GetNodes():
                nodes_x_coordinate.append(node.X)
                nodes_y_coordinate.append(node.Y)
                if (len(nodes_x_coordinate) >= 4):
                    break

            # Create a polygon for the element
            points_list = list(zip(nodes_x_coordinate, nodes_y_coordinate))     
            element_polygon = Polygon(points_list)
            
            # Check for intersection between the embedded body and an element
            intersection = element_polygon.intersection(boundary_polygon)

            if (intersection.area > 0.0 and (np.abs(element_polygon.area - intersection.area)/element_polygon.area) < 1.0 and (np.abs(element_polygon.area - intersection.area)/element_polygon.area) > 0.0):
                intersected_elements_sub_model_part.AddElement(element)
                # active_elements_sub_model_part.AddElement(element)
            elif (intersection.area > 0.0 and (np.abs(element_polygon.area - intersection.area)/element_polygon.area) <= 1.0 and self.keep_external_domain == False):
                active_elements_sub_model_part.AddElement(element)
            elif ((np.abs(element_polygon.area - intersection.area)/element_polygon.area) > 0.95 and self.keep_external_domain == True):
                active_elements_sub_model_part.AddElement(element)
            else:
                inactive_elements_sub_model_part.AddElement(element)

        # Desactivate the elements in inactive_elements_sub_model_part. This elements are desactivated because they shouldn't be assembled in the LHS and RHS
        for element in inactive_elements_sub_model_part.Elements:
            element.Set(KM.ACTIVE, False)
        
