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

        KratosGeo.FillCfdModelpartUtilities(self.CfdModelPart).FillPartsFluid(self.ModelPart, elem_sub_model_name)



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
        
        KratosGeo.FillCfdModelpartUtilities(self.CfdModelPart).FillInlet(self.ModelPart, cond_sub_model_name)


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
        
        KratosGeo.FillCfdModelpartUtilities(self.CfdModelPart).FillOutlet(self.ModelPart, cond_sub_model_name)


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
        
        KratosGeo.FillCfdModelpartUtilities(self.CfdModelPart).FillSlip(self.ModelPart, cond_sub_model_name)


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
        
        KratosGeo.FillCfdModelpartUtilities(self.CfdModelPart).FillNoslip(self.ModelPart, cond_sub_model_name)
