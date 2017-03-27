 
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library 
import KratosMultiphysics 
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication

KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ContactProcess(Model, settings["Parameters"])

class ContactProcess(KratosMultiphysics.Process):
  
    def __init__(self,model_part,params):
        
        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"             : "",
            "origin_model_part_name"      : "",
            "destination_model_part_name" : "",
            "contact_type"                : "MortarMethod",
            "search_factor"               : 1.5,
            "active_check_factor"         : 0.01,
            "max_number_results"          : 1000,
            "augmentation_normal"         : 0.0e0,
            "augmentation_tangent"        : 0.0e0,
            "simplify_geometry"           : false,
            "type_search"                 : "InRadius",
            "integration_order"           : 5
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
   
        self.main_model_part = model_part[self.params["model_part_name"].GetString()]

        self.dimension = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        self.o_model_part = model_part[self.params["origin_model_part_name"].GetString()]
        self.d_model_part = model_part[self.params["destination_model_part_name"].GetString()]
        
        self.o_interface = self.o_model_part.GetSubModelPart("Interface")
        self.d_interface = self.d_model_part.GetSubModelPart("Interface")
        
        self.search_factor            = self.params["search_factor"].GetDouble() 
        self.active_check_factor      = self.params["active_check_factor"].GetDouble() 
        self.max_number_results       = self.params["max_number_results"].GetInt() 
        self.augmentation_normal      = self.params["augmentation_normal"].GetDouble()
        self.augmentation_tangent     = self.params["augmentation_tangent"].GetDouble()
        self.simplify_geometry        = self.params["simplify_geometry"].GetBool()
        self.integration_order        = self.params["integration_order"].GetInt() 
        if self.params["type_search"].GetString() == "InRadius":
             self.type_search = 0
        
    def ExecuteInitialize(self):
        
        # Appending the conditions created to the computing_model_part
        computing_model_part = self.main_model_part.GetSubModelPart("computing_domain")
        computing_model_part.CreateSubModelPart("Contact")
        interface_computing_model_part = computing_model_part.GetSubModelPart("Contact")
            
        for node in self.o_interface.Nodes:
            interface_computing_model_part.AddNode(node, 0)  
            node.Set(KratosMultiphysics.INTERFACE,True)
        del(node)
        
        for node in self.d_interface.Nodes:
            interface_computing_model_part.AddNode(node, 0)
            node.Set(KratosMultiphysics.INTERFACE,True)
        del(node)
        
        self.Preprocess  = KratosMultiphysics.ContactStructuralMechanicsApplication.InterfacePreprocessCondition(self.main_model_part)
        
        if self.params["contact_type"].GetString() == "MortarMethod":
            condition_name = "MortarContact"
        
        #print("MODEL PART BEFORE CREATING INTERFACE")
        #print(self.main_model_part) 
        
        final_string = ""
        
        # It should create the conditions automatically
        if (self.dimension == 2):
            self.Preprocess.GenerateInterfacePart2D(self.o_model_part, self.o_interface, condition_name, final_string, self.simplify_geometry) 
            self.Preprocess.GenerateInterfacePart2D(self.d_model_part, self.d_interface, condition_name, final_string, self.simplify_geometry) 
        else:
            self.Preprocess.GenerateInterfacePart3D(self.o_model_part, self.o_interface, condition_name, final_string, self.simplify_geometry) 
            self.Preprocess.GenerateInterfacePart3D(self.d_model_part, self.d_interface, condition_name, final_string, self.simplify_geometry) 

        #print("MODEL PART AFTER CREATING INTERFACE")
        #print(self.main_model_part)
        
        for cond in self.o_interface.Conditions:
            interface_computing_model_part.AddCondition(cond)    
        del(cond)
        
        for cond in self.d_interface.Conditions:
            interface_computing_model_part.AddCondition(cond)  
        del(cond)

        self.contact_search = KratosMultiphysics.ContactStructuralMechanicsApplication.TreeContactSearch(self.o_interface, self.d_interface, self.max_number_results)
        
        if self.params["contact_type"].GetString() == "MortarMethod":
            ZeroVector = KratosMultiphysics.Vector(3) 
            ZeroVector[0] = 0.0
            ZeroVector[1] = 0.0
            ZeroVector[2] = 0.0
            
            # Initilialize weighted variables and LM
            for node in self.d_interface.Nodes:
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_GAP, 0.0)
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_SLIP, 0.0)
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_FRICTION, 0.0)
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.AUXILIAR_ACTIVE, False)
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.AUXILIAR_SLIP, False)
                node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)
                node.SetValue(KratosMultiphysics.NORMAL, ZeroVector)
                node.SetValue(KratosMultiphysics.TANGENT_XI, ZeroVector)
                node.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector)
                #node.Set(KratosMultiphysics.SLAVE, True)
            del node
            
            # Setting the master conditions 
            for cond in self.o_interface.Conditions:
                cond.SetValue(KratosMultiphysics.NORMAL, ZeroVector) 
                cond.SetValue(KratosMultiphysics.TANGENT_XI, ZeroVector) 
                cond.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector) 
                #cond.Set(KratosMultiphysics.MASTER, True) # TODO: This is not supposed o be necessary
            del cond
            
            # Setting the slave conditions 
            for cond in self.d_interface.Conditions:
                cond.SetValue(KratosMultiphysics.NORMAL, ZeroVector) 
                cond.SetValue(KratosMultiphysics.TANGENT_XI, ZeroVector) 
                cond.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector) 
            del cond
            
            self.contact_search.CreatePointListMortar()
            self.contact_search.InitializeMortarConditions(self.active_check_factor, self.augmentation_normal, self.augmentation_tangent, self.integration_order)
        
    def ExecuteBeforeSolutionLoop(self):
        self.contact_search.TotalClearMortarConditions();
    
    def ExecuteInitializeSolutionStep(self):
        #for cond in self.d_interface.Conditions:
            #print(cond.Is(KratosMultiphysics.ACTIVE))
        
        if self.params["contact_type"].GetString() == "MortarMethod":    
            self.contact_search.UpdateMortarConditions(self.search_factor, self.type_search)
            #self.contact_search.CheckMortarConditions()
            
        #for cond in self.d_interface.Conditions:
            #print(cond.Is(KratosMultiphysics.ACTIVE))
        
    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        if self.params["contact_type"].GetString() == "MortarMethod":
            self.contact_search.UpdatePointListMortar()
            self.contact_search.PartialClearMortarConditions()
            
        #for cond in self.d_interface.Conditions:
            #print(cond.Is(KratosMultiphysics.ACTIVE))
            
    def ExecuteFinalize(self):
        pass
