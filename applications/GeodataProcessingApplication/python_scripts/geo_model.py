"""
    # edit:     19 September 2019 -> added FillPartsFluid, FillInlet, FillOutlet, FillSlip, FillNoslip functions
    # edit:     20 September 2019
"""

import KratosMultiphysics
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.MeshingApplication as KratosMesh

from geo_processor import GeoProcessor

class GeoModel( GeoProcessor ):

    def __init__( self ):
        super(GeoModel, self).__init__()

        self.HasModelPart = False
        self.HasCfdModelPart = False


    def GenerateCfdModelPart(self):

        # current_model = KratosMultiphysics.Model()
        current_model = self.ModelPart.GetModel()
        
        self.CfdModelPart = current_model.CreateModelPart("CFD_model_part")
        
        # self.CfdModelPart.CreateSubModelPart("Parts_Fluid")
        # self.CfdModelPart.CreateSubModelPart("Inlet")
        # self.CfdModelPart.CreateSubModelPart("Outlet")
        # self.CfdModelPart.CreateSubModelPart("Slip")
        # self.CfdModelPart.CreateSubModelPart("NoSlip")

        # we set the DENSITY and DYNAMIC_VISCOSITY values
        prop = self.CfdModelPart.GetProperties()[0]
        prop.SetValue(KratosMultiphysics.DENSITY, 1)
        prop.SetValue(KratosMultiphysics.DYNAMIC_VISCOSITY, 0.002)

        self.HasCfdModelPart = True
    

    def GetGeoCfdModelPart(self):
        
        if self.HasCfdModelPart:
            return self.CfdModelPart
        else:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "No CFD model part can be returned")
    

    def FillPartsFluid(self, elem_sub_model_name):
        # we fill element into Parts_Fluid sub model part
        # we perform a check if the nodes and the elements are already into the main model part
        # and, if not, we create it and we added it into sub model part

        if not self.HasCfdModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "CFD model part has to be set, first.")
            return
        
        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "A model part has to be set, first.")
            return

        if not self.CfdModelPart.HasSubModelPart("Parts_Fluid"):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"Parts_Fluid\" created right now!")
            self.CfdModelPart.CreateSubModelPart("Parts_Fluid")

        nodes_main = self.CfdModelPart.GetNodes()       # array with nodes into main model part
        elems_main = self.CfdModelPart.GetElements()    # array with elements into main model part

        fluid_sub_model = self.CfdModelPart.GetSubModelPart("Parts_Fluid")

        elem_sub_model = self.ModelPart.GetSubModelPart(elem_sub_model_name)
        list_nodes = []
        list_elems = []
        for elem in elem_sub_model.Elements:
            nodes = elem.GetNodes()
            for node in nodes:
                if not node in nodes_main:
                    # if the node it is not into main model part, we create it
                    self.CfdModelPart.CreateNewNode(node.Id, node.X, node.Y, node.Z)
                # we add the node id into list
                list_nodes.append(node.Id)
            
            if not elem in elems_main:
                # if the element it is not into main model part, we create it
                self.CfdModelPart.CreateNewElement("Element3D4N", elem.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id, nodes[3].Id], self.CfdModelPart.GetProperties()[0])
            # we add the element id into a list
            list_elems.append(elem.Id)
        
        list_nodes = list(set(list_nodes))      # to avoid duplicate
        list_elems = list(set(list_elems))      # to avoid duplicate
        fluid_sub_model.AddNodes(list_nodes)
        fluid_sub_model.AddElements(list_elems)


    def FillInlet(self, cond_sub_model_name):
        # we fill conditions into Inlet sub model part
        # we perform a check if the nodes and the condiions are already into the main model part
        # and, if not, we create it and we added it into sub model part

        if not self.HasCfdModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "CFD model part has to be set, first.")
            return

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "A model part has to be set, first.")
            return

        if not self.CfdModelPart.HasSubModelPart("Inlet"):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"Inlet\" created right now!")
            self.CfdModelPart.CreateSubModelPart("Inlet")
        
        nodes_main = self.CfdModelPart.GetNodes()       # array with nodes into main model part
        conds_main = self.CfdModelPart.GetConditions()  # array with conditions into main model part

        inlet_sub_model = self.CfdModelPart.GetSubModelPart("Inlet")

        cond_sub_model = self.ModelPart.GetSubModelPart(cond_sub_model_name)
        list_nodes = []
        list_conds = []
        for cond in cond_sub_model.Conditions:
            nodes = cond.GetNodes()
            for node in nodes:
                if not node in nodes_main:
                    # if the node it is not into main model part, we create it
                    self.CfdModelPart.CreateNewNode(node.Id, node.X, node.Y, node.Z)
                # we add the node id into list
                list_nodes.append(node.Id)
            
            if not cond in conds_main:
                # if the condition it is not into main model part, we create it
                self.CfdModelPart.CreateNewCondition("SurfaceCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], self.CfdModelPart.GetProperties()[0])
            # we add the condition id into a list
            list_conds.append(cond.Id)
        
        list_nodes = list(set(list_nodes))      # to avoid duplicate
        list_conds = list(set(list_conds))      # to avoid duplicate
        inlet_sub_model.AddNodes(list_nodes)
        inlet_sub_model.AddConditions(list_conds)


    def FillOutlet(self, cond_sub_model_name):
        # we fill conditions into Outlet sub model part
        # we perform a check if the nodes and the condiions are already into the main model part
        # and, if not, we create it and we added it into sub model part
        
        if not self.HasCfdModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "CFD model part has to be set, first.")
            return

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "A model part has to be set, first.")
            return

        if not self.CfdModelPart.HasSubModelPart("Outlet"):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"Outlet\" created right now!")
            self.CfdModelPart.CreateSubModelPart("Outlet")
        
        nodes_main = self.CfdModelPart.GetNodes()       # array with nodes into main model part
        conds_main = self.CfdModelPart.GetConditions()  # array with conditions into main model part

        outlet_sub_model = self.CfdModelPart.GetSubModelPart("Outlet")

        cond_sub_model = self.ModelPart.GetSubModelPart(cond_sub_model_name)
        list_nodes = []
        list_conds = []
        for cond in cond_sub_model.Conditions:
            nodes = cond.GetNodes()
            for node in nodes:
                if not node in nodes_main:
                    # if the node it is not into main model part, we create it
                    self.CfdModelPart.CreateNewNode(node.Id, node.X, node.Y, node.Z)
                # we add the node id into list
                list_nodes.append(node.Id)
            
            if not cond in conds_main:
                # if the condition it is not into main model part, we create it
                self.CfdModelPart.CreateNewCondition("SurfaceCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], self.CfdModelPart.GetProperties()[0])
            # we add the condition id into a list
            list_conds.append(cond.Id)
        
        list_nodes = list(set(list_nodes))      # to avoid duplicate
        list_conds = list(set(list_conds))      # to avoid duplicate
        outlet_sub_model.AddNodes(list_nodes)
        outlet_sub_model.AddConditions(list_conds)


    def FillSlip(self, cond_sub_model_name):
        # we fill conditions into Slip sub model part
        # we perform a check if the nodes and the condiions are already into the main model part
        # and, if not, we create it and we added it into sub model part
        
        if not self.HasCfdModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "CFD model part has to be set, first.")
            return

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "A model part has to be set, first.")
            return

        if not self.CfdModelPart.HasSubModelPart("Slip"):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"Slip\" created right now!")
            self.CfdModelPart.CreateSubModelPart("Slip")
        
        nodes_main = self.CfdModelPart.GetNodes()       # array with nodes into main model part
        conds_main = self.CfdModelPart.GetConditions()  # array with conditions into main model part

        slip_sub_model = self.CfdModelPart.GetSubModelPart("Slip")

        cond_sub_model = self.ModelPart.GetSubModelPart(cond_sub_model_name)
        list_nodes = []
        list_conds = []
        for cond in cond_sub_model.Conditions:
            nodes = cond.GetNodes()
            for node in nodes:
                if not node in nodes_main:
                    # if the node it is not into main model part, we create it
                    self.CfdModelPart.CreateNewNode(node.Id, node.X, node.Y, node.Z)
                # we add the node id into list
                list_nodes.append(node.Id)
            
            if not cond in conds_main:
                # if the condiion it is not into main model part, we create it
                self.CfdModelPart.CreateNewCondition("SurfaceCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], self.CfdModelPart.GetProperties()[0])
            # we add the condition id into a list
            list_conds.append(cond.Id)
        
        list_nodes = list(set(list_nodes))      # to avoid duplicate
        list_conds = list(set(list_conds))      # to avoid duplicate
        slip_sub_model.AddNodes(list_nodes)
        slip_sub_model.AddConditions(list_conds)


    def FillNoslip(self, cond_sub_model_name):
        # we fill conditions into NoSlip sub model part
        # we perform a check if the nodes and the condiions are already into the main model part
        # and, if not, we create it and we added it into sub model part
        
        if not self.HasCfdModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "CFD model part has to be set, first.")
            return

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "A model part has to be set, first.")
            return

        if not self.CfdModelPart.HasSubModelPart("NoSlip"):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "Sub Model Part \"NoSlip\" created right now!")
            self.CfdModelPart.CreateSubModelPart("NoSlip")
        
        nodes_main = self.CfdModelPart.GetNodes()       # array with nodes into main model part
        conds_main = self.CfdModelPart.GetConditions()  # array with conditions into main model part

        noslip_sub_model = self.CfdModelPart.GetSubModelPart("NoSlip")

        cond_sub_model = self.ModelPart.GetSubModelPart(cond_sub_model_name)
        list_nodes = []
        list_conds = []
        for cond in cond_sub_model.Conditions:
            nodes = cond.GetNodes()
            for node in nodes:
                if not node in nodes_main:
                    # if the node it is not into main model part, we create it
                    self.CfdModelPart.CreateNewNode(node.Id, node.X, node.Y, node.Z)
                # we add the node id into list
                list_nodes.append(node.Id)
            
            if not cond in conds_main:
                # if the condition it is not into main model part, we create it
                self.CfdModelPart.CreateNewCondition("SurfaceCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], self.CfdModelPart.GetProperties()[0])
            # we add the condition id into a list
            list_conds.append(cond.Id)
        
        list_nodes = list(set(list_nodes))      # to avoid duplicate
        list_conds = list(set(list_conds))      # to avoid duplicate
        noslip_sub_model.AddNodes(list_nodes)
        noslip_sub_model.AddConditions(list_conds)
