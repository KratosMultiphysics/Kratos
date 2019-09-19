"""
    # edit:     19 September 2019
"""

import KratosMultiphysics
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.MeshingApplication as KratosMesh

from geo_processor import GeoProcessor

class GeoModel( GeoProcessor ):

    def __init__( self ):
        super(GeoMesher, self).__init__()

        self.HasModelPart = False
        self.HasExtrusionHeight = False

    # def DefineElementType( self, type_of_element ):

    #     if ( not self.HasModelPart ):
    #         KratosMultiphysics.Logger.PrintWarning("GeoModel", "No model part has been set for this function")
    #         return

    def GenerateCfdModelPart(self):

        current_model = self.ModelPart.GetModel()

        if current_model.HasModelPart("CFD_model_part"):
            current_model.DeleteModelPart("CFD_model_part")

        self.CfdModelPart = current_model.CreateModelPart("CFD_model_part")
        self.CfdModelPart.CreateSubModelPart("Parts_Fluid")
        self.CfdModelPart.CreateSubModelPart("Inlet")
        self.CfdModelPart.CreateSubModelPart("Outlet")
        self.CfdModelPart.CreateSubModelPart("Slip")
        self.CfdModelPart.CreateSubModelPart("NoSlip")

        self.HasModelPart = True
    

    def FillPartsFluid(self, model_part, elem_sub_model):
        # we fill element into Parts_Fluid sub model part
        # we perform a check if the nodes and the elements are already into the main model part
        # and, if not, we create it and we added it into sub model part

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Model part has to be set, first.")
            return

        if not self.CfdModelPart.HasSubModelPart("Parts_Fluid"):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"Parts_Fluid\" created right now!")
            self.CfdModelPart.CreateSubModelPart("Parts_Fluid")

        nodes_main = self.CfdModelPart.GetNodes()       # array with nodes into main model part
        elems_main = self.CfdModelPart.GetElements()    # array with elements into main model part

        parts_fluid = self.CfdModelPart.GetSubModePart("Parts_Fluid")

        for elem in self.CfdModelPart.Elements:
            nodes = elem.GetNodes()
            for node in nodes:
                if not node in nodes_main:
                    # if the node it is not into main model part, we create it
                    self.CfdModelPart.CreateNewNode(node.Id, node.X, node.Y, node.Z)
                # we add the node into sub model part
                parts_fluid.AddNode(node, 0)
            
            if not elem in elems_main:
                # if the element it is not into main model part, we create it
                self.CfdModelPart.CreateNewElement("Element3D4N", elem.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id, nodes[3].Id], self.CfdModelPart.GetProperties()[0])
            # we add the element into sub model part
            parts_fluid.AddElement(elem, 0)




    def FillInlet(self, model_part, elem_sub_model):
        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Model part has to be set, first.")
            return

        if not self.CfdModelPart.HasSubModelPart("Inlet"):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"Inlet\" created right now!")
            self.CfdModelPart.CreateSubModelPart("Inlet")

    def FillOutlet(self, model_part, elem_sub_model):
        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Model part has to be set, first.")
            return

        if not self.CfdModelPart.HasSubModelPart("Outlet"):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"Outlet\" created right now!")
            self.CfdModelPart.CreateSubModelPart("Outlet")

    def FillSlip(self, model_part, elem_sub_model):
        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Model part has to be set, first.")
            return

        if not self.CfdModelPart.HasSubModelPart("Slip"):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"Slip\" created right now!")
            self.CfdModelPart.CreateSubModelPart("Slip")

    def FillNoslip(self, model_part, elem_sub_model):
        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Model part has to be set, first.")
            return

        if not self.CfdModelPart.HasSubModelPart("NoSlip"):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"NoSlip\" created right now!")
            self.CfdModelPart.CreateSubModelPart("NoSlip")