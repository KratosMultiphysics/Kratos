 
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library 
import KratosMultiphysics 
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication

KratosMultiphysics.CheckForPreviousImport()

def CalculateLastIdCondition(model_part):
    cond_id = 0
    for condition in model_part.Conditions:
        cond_id += 1

    return cond_id

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MeshTyingProcess(Model, settings["Parameters"])

class MeshTyingProcess(KratosMultiphysics.Process):
  
    def __init__(self,model_part,params):
        
        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name"             : "",
            "origin_model_part_name"      : "",
            "destination_model_part_name" : "",
            "type_variable"               : "Components",
            "geometry_element"            : "Quadrilateral",
            "search_factor"               : 1.5,
            "active_check_factor"         : 0.01,
            "max_number_results"          : 1000,
            "type_search"                 : "InRadius",
            "integration_order"           : 1
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
   
        self.main_model_part = model_part[self.params["model_part_name"].GetString()]

        self.o_model_part = model_part[self.params["origin_model_part_name"].GetString()]
        self.d_model_part = model_part[self.params["destination_model_part_name"].GetString()]
        
        self.o_interface = self.o_model_part.GetSubModelPart("Interface")
        self.d_interface = self.d_model_part.GetSubModelPart("Interface")
        
        self.search_factor            = self.params["search_factor"].GetDouble() 
        self.active_check_factor      = self.params["active_check_factor"].GetDouble() 
        self.max_number_results       = self.params["max_number_results"].GetInt() 

        self.type_variable            = self.params["type_variable"].GetString() 
        self.geometry_element         = self.params["geometry_element"].GetString() 
        self.integration_order        = self.params["integration_order"].GetInt() 
        if self.params["type_search"].GetString() == "InRadius":
             self.type_search = 0
        
    def ExecuteInitialize(self):
        
        # Appending the conditions created to the computing_model_part
        computing_model_part = self.main_model_part.GetSubModelPart("computing_domain")
        computing_model_part.CreateSubModelPart("MeshTying")
        interface_computing_model_part = computing_model_part.GetSubModelPart("MeshTying")
            
        for node in self.o_interface.Nodes:
            interface_computing_model_part.AddNode(node, 0)  
            node.Set(KratosMultiphysics.INTERFACE, True)
        del(node)
        
        for node in self.d_interface.Nodes:
            interface_computing_model_part.AddNode(node, 0)
            node.Set(KratosMultiphysics.INTERFACE, True)
        del(node)
        
        self.Preprocess = KratosMultiphysics.ContactStructuralMechanicsApplication.InterfacePreprocessCondition()
        
        condition_name = "MeshTyingMortar"
        
        #print("MODEL PART BEFORE CREATING INTERFACE")
        #print(self.main_model_part) 
        
        # It should create the conditions automatically
        initial_id = CalculateLastIdCondition(self.main_model_part)
        self.Preprocess.GenerateInterfacePart(self.o_model_part, self.o_interface, condition_name, initial_id, self.geometry_element + self.type_variable, False) 
        initial_id = CalculateLastIdCondition(self.main_model_part)
        self.Preprocess.GenerateInterfacePart(self.d_model_part, self.d_interface, condition_name, initial_id, self.geometry_element + self.type_variable, False) 

        #print("MODEL PART AFTER CREATING INTERFACE")
        #print(self.main_model_part)
        
        for cond in self.o_interface.Conditions:
            interface_computing_model_part.AddCondition(cond)    
        del(cond)
        
        for cond in self.d_interface.Conditions:
            interface_computing_model_part.AddCondition(cond)  
        del(cond)

        self.contact_search = KratosMultiphysics.ContactStructuralMechanicsApplication.TreeContactSearch(self.o_interface, self.d_interface, self.max_number_results)
        
        ZeroVector = KratosMultiphysics.Vector(3) 
        ZeroVector[0] = 0.0
        ZeroVector[1] = 0.0
        ZeroVector[2] = 0.0
        
        # Initilialize weighted variables and LM
        for node in self.d_interface.Nodes:
            node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)
            node.SetValue(KratosMultiphysics.NORMAL,      ZeroVector)
            node.SetValue(KratosMultiphysics.TANGENT_XI,  ZeroVector)
            node.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector)
            #node.Set(KratosMultiphysics.SLAVE, True)
        del node
        
        # Setting the master conditions 
        for cond in self.o_interface.Nodes:
            cond.SetValue(KratosMultiphysics.NORMAL,      ZeroVector) 
            cond.SetValue(KratosMultiphysics.TANGENT_XI,  ZeroVector) 
            cond.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector) 
            #cond.Set(KratosMultiphysics.MASTER, True) # TODO: This is not supposed o be necessary
        del cond
        
        # Setting the slave conditions 
        for cond in self.d_interface.Nodes:
            cond.SetValue(KratosMultiphysics.NORMAL,      ZeroVector) 
            cond.SetValue(KratosMultiphysics.TANGENT_XI,  ZeroVector) 
            cond.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector) 
        del cond
        
        self.contact_search.CreatePointListMortar()
        self.contact_search.InitializeMeshTyingMortarConditions(self.active_check_factor, self.integration_order)
        
    def ExecuteBeforeSolutionLoop(self):
        if (self.type_variable == "Scalar"):
            self.contact_search.TotalClearMeshTyingScalarMortarConditions()
        else:
            self.contact_search.TotalClearMeshTyingComponentsMortarConditions()
            
        self.contact_search.UpdateMortarConditions(self.search_factor, self.type_search)
    
    def ExecuteInitializeSolutionStep(self):
        pass
        #for cond in self.d_interface.Conditions:
            #print(cond.Is(KratosMultiphysics.ACTIVE))
        
        #self.contact_search.UpdateMortarConditions(self.search_factor, self.type_search)
        #self.contact_search.CheckMortarConditions()
            
        #for cond in self.d_interface.Conditions:
            #print(cond.Is(KratosMultiphysics.ACTIVE))
        
    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass
        #self.contact_search.UpdatePointListMortar()
        #if (self.type_variable == "Scalar"):
            #self.contact_search.PartialClearMeshTyingMortarScalarConditions()
        #else:
            #self.contact_search.PartialClearMeshTyingMortarComponentsConditions()
            
        #for cond in self.d_interface.Conditions:
            #print(cond.Is(KratosMultiphysics.ACTIVE))
            
    def ExecuteFinalize(self):
        pass
