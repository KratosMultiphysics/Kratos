 
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
    return MeshTyingProcess(Model, settings["Parameters"])

class MeshTyingProcess(KratosMultiphysics.Process):
  
    def __init__(self,model_part,params):
        
        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "mesh_id"                     : 0,
            "model_part_name"             : "Structure",
            "computing_model_part_name"   : "computing_domain",
            "mesh_tying_model_part"          : "Contact_Part",
            "assume_master_slave"         : "",
            "type_variable"               : "Components",
            "geometry_element"            : "Quadrilateral",
            "search_factor"               : 1.5,
            "active_check_factor"         : 0.01,
            "max_number_results"          : 1000,
            "type_search"                 : "InRadius",
            "integration_order"           : 2
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
   
        self.main_model_part = model_part[self.params["model_part_name"].GetString()]
        self.computing_model_part_name = self.params["computing_model_part_name"].GetString()

        self.dimension = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        self.mesh_tying_model_part = model_part[self.params["mesh_tying_model_part"].GetString()]
        
        if (self.params["assume_master_slave"].GetString() != ""):
            for node in model_part[self.params["assume_master_slave"].GetString()].Nodes:
                node.Set(KratosMultiphysics.SLAVE, True)
        
        self.search_factor            = self.params["search_factor"].GetDouble() 
        self.active_check_factor      = self.params["active_check_factor"].GetDouble() 
        self.max_number_results       = self.params["max_number_results"].GetInt() 

        self.type_variable            = self.params["type_variable"].GetString() 
        self.geometry_element         = self.params["geometry_element"].GetString() 
        self.integration_order        = self.params["integration_order"].GetInt() 
        self.type_search              = self.params["type_search"].GetString()
        
    def ExecuteInitialize(self):
        
        # Appending the conditions created to the self.main_model_part
        computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        if (computing_model_part.HasSubModelPart("Contact")):
            interface_model_part = computing_model_part.GetSubModelPart("Contact")
        else:
            interface_model_part = computing_model_part.CreateSubModelPart("Contact")
            
        for prop in computing_model_part.GetProperties():
            prop[KratosMultiphysics.ContactStructuralMechanicsApplication.INTEGRATION_ORDER_CONTACT] = self.integration_order
            
        for node in self.mesh_tying_model_part.Nodes:
            node.Set(KratosMultiphysics.INTERFACE, True)
        del(node)
        
        self.Preprocess = KratosMultiphysics.ContactStructuralMechanicsApplication.InterfacePreprocessCondition(self.main_model_part)
        
        condition_name = "MeshTyingMortar"
        
        #print("MODEL PART BEFORE CREATING INTERFACE")
        #print(self.main_model_part) 
        
        # It should create the conditions automatically
        if (self.dimension == 2):
            self.Preprocess.GenerateInterfacePart2D(computing_model_part, self.mesh_tying_model_part, condition_name, self.geometry_element + self.type_variable, False) 
        else:
            self.Preprocess.GenerateInterfacePart3D(computing_model_part, self.mesh_tying_model_part, condition_name, self.geometry_element + self.type_variable, False) 

        # When all conditions are simultaneously master and slave
        if (self.params["assume_master_slave"].GetString() == ""):
            for cond in self.mesh_tying_model_part.Conditions:
                cond.Set(KratosMultiphysics.SLAVE, True)

        #print("MODEL PART AFTER CREATING INTERFACE")
        #print(self.main_model_part)
        
        # We copy the conditions to the ContactSubModelPart
        for cond in self.mesh_tying_model_part.Conditions:
            interface_model_part.AddCondition(cond)    
        del(cond)
        for node in self.mesh_tying_model_part.Nodes:
            interface_model_part.AddNode(node, 0)    
        del(node)
        
        self.contact_search = KratosMultiphysics.ContactStructuralMechanicsApplication.TreeContactSearch(computing_model_part, self.max_number_results, self.active_check_factor, self.type_search)
        
        ZeroVector = KratosMultiphysics.Vector(3) 
        ZeroVector[0] = 0.0
        ZeroVector[1] = 0.0
        ZeroVector[2] = 0.0
        
        # Initilialize weighted variables and LM
        for node in self.mesh_tying_model_part.Nodes:
            if (self.type_variable == "Scalar"):
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_SCALAR_RESIDUAL, 0.0)
            else:
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_VECTOR_RESIDUAL, ZeroVector)
            node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)
            node.SetValue(KratosMultiphysics.NORMAL,      ZeroVector)
            node.SetValue(KratosMultiphysics.TANGENT_XI,  ZeroVector)
            node.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector)
        del node
        
        # Setting the conditions 
        for cond in self.mesh_tying_model_part.Conditions:
            cond.SetValue(KratosMultiphysics.NORMAL,      ZeroVector) 
            cond.SetValue(KratosMultiphysics.TANGENT_XI,  ZeroVector) 
            cond.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector) 
        del cond
        
        self.contact_search.CreatePointListMortar()
        self.contact_search.InitializeMortarConditions()
        
    def ExecuteBeforeSolutionLoop(self):
        if (self.type_variable == "Scalar"):
            self.contact_search.TotalClearScalarMortarConditions()
        else:
            self.contact_search.TotalClearComponentsMortarConditions()
            
        self.contact_search.UpdateMortarConditions(self.search_factor)
    
    def ExecuteInitializeSolutionStep(self):
        pass
        #for cond in self.d_interface.Conditions:
            #print(cond.Is(KratosMultiphysics.ACTIVE))
        
        #self.contact_search.UpdateMortarConditions(self.search_factor)
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
            #self.contact_search.PartialClearScalarMortarConditions()
        #else:
            #self.contact_search.PartialClearComponentsMortarConditions()
            
        #for cond in self.d_interface.Conditions:
            #print(cond.Is(KratosMultiphysics.ACTIVE))
            
    def ExecuteFinalize(self):
        pass
