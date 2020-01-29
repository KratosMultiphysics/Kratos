import KratosMultiphysics as KratosMultiphysics
import KratosMultiphysics.CableNetApplication as CableNetApplication

from KratosMultiphysics import Logger

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return EdgeCableElementProcess(Model, settings["Parameters"])



class EdgeCableElementProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        default_settings = KratosMultiphysics.Parameters("""
        {
            "edge_sub_model_part_name"  : "Structure.example_part",
            "element_type"              : "cable",
            "node_id_order"             : [1,2,3],
            "element_id"                : 1,
            "property_id"               : 1
        }
        """)
        default_settings.ValidateAndAssignDefaults(settings)

        # The computing model part
        self.edge_model_part = Model[settings["edge_sub_model_part_name"].GetString()]
        self.edge_cable_element_process = CableNetApplication.EdgeCableElementProcess(self.edge_model_part, settings)


    def ExecuteInitialize(self):
        self.edge_cable_element_process.ExecuteInitialize()
        Logger.PrintInfo("Initialized","EdgeCableElementProcess")

