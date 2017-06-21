 
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library 
import KratosMultiphysics 
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

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
            "bucket_size"                 : 4,
            "filter_candidates"           : true,
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

        self.type_variable            = self.params["type_variable"].GetString() 
        self.geometry_element         = self.params["geometry_element"].GetString() 
        self.integration_order        = self.params["integration_order"].GetInt() 
        
    def ExecuteInitialize(self):
        
        # Appending the conditions created to the self.main_model_part
        computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        if (computing_model_part.HasSubModelPart("Contact")):
            interface_model_part = computing_model_part.GetSubModelPart("Contact")
        else:
            interface_model_part = computing_model_part.CreateSubModelPart("Contact")
            
        for prop in computing_model_part.GetProperties():
            prop[ContactStructuralMechanicsApplication.INTEGRATION_ORDER_CONTACT] = self.integration_order
            
        for node in self.mesh_tying_model_part.Nodes:
            node.Set(KratosMultiphysics.INTERFACE, True)
        del(node)
        
        self.Preprocess = ContactStructuralMechanicsApplication.InterfacePreprocessCondition(self.main_model_part)
        
        condition_name = "MeshTyingMortar"
        
        #print("MODEL PART BEFORE CREATING INTERFACE")
        #print(self.main_model_part) 
        
        # It should create the conditions automatically
        interface_parameters = KratosMultiphysics.Parameters("""{"condition_name": "", "final_string": "", "simplify_geometry": false}""")
        interface_parameters["condition_name"].SetString(condition_name)
        interface_parameters["final_string"].SetString(self.geometry_element + self.type_variable)
        if (self.dimension == 2):
            self.Preprocess.GenerateInterfacePart2D(computing_model_part, self.mesh_tying_model_part, interface_parameters) 
        else:
            self.Preprocess.GenerateInterfacePart3D(computing_model_part, self.mesh_tying_model_part, interface_parameters) 

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
        
        # Creating the search
        search_parameters = KratosMultiphysics.Parameters("""{}""")
        search_parameters.AddValue("type_search",self.params["type_search"])
        search_parameters.AddValue("allocation_size",self.params["max_number_results"])
        search_parameters.AddValue("bucket_size",self.params["bucket_size"])
        search_parameters.AddValue("search_factor",self.params["search_factor"])
        search_parameters.AddValue("active_check_factor",self.params["active_check_factor"])
        search_parameters.AddValue("filter_candidates",self.params["filter_candidates"])
        self.contact_search = ContactStructuralMechanicsApplication.TreeContactSearch(computing_model_part, search_parameters)
        
        ZeroVector = KratosMultiphysics.Vector(3) 
        ZeroVector[0] = 0.0
        ZeroVector[1] = 0.0
        ZeroVector[2] = 0.0
        
        # Initilialize weighted variables and LM
        for node in self.mesh_tying_model_part.Nodes:
            if (self.type_variable == "Scalar"):
                node.SetSolutionStepValue(ContactStructuralMechanicsApplication.WEIGHTED_SCALAR_RESIDUAL, 0.0)
            else:
                node.SetSolutionStepValue(ContactStructuralMechanicsApplication.WEIGHTED_VECTOR_RESIDUAL, ZeroVector)
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
            
        self.contact_search.UpdateMortarConditions()
    
    def ExecuteInitializeSolutionStep(self):
        pass
        
    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass
            
    def ExecuteFinalize(self):
        pass
