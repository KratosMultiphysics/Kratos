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
    return ALMContactProcess(Model, settings["Parameters"])

# NOTE: To gain efficency we can do the distinction between MASTER and SLAVE (the automatic search is expensive). Maybe looking to the submodelparts the nodes belong
class ALMContactProcess(KratosMultiphysics.Process):
  
    def __init__(self,model_part,params):
        
        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "mesh_id"                     : 0,
            "model_part_name"             : "Structure",
            "computing_model_part_name"   : "computing_domain",
            "contact_model_part"          : "Contact_Part",
            "assume_master_slave"         : "",
            "contact_type"                : "Frictionless",
            "search_factor"               : 1.5,
            "active_check_factor"         : 0.01,
            "max_number_results"          : 1000,
            "normal_variation"            : false,
            "manual_ALM"                  : false,
            "manual_ALM_parameters"       :
            {
                "penalty"                   : 0.0,
                "scale_factor"              : 1.0
            },
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

        self.contact_model_part = model_part[self.params["contact_model_part"].GetString()]
        
        if (self.params["assume_master_slave"].GetString() != ""):
            for node in model_part[self.params["assume_master_slave"].GetString()].Nodes:
                node.Set(KratosMultiphysics.SLAVE, True)
        
        self.search_factor            = self.params["search_factor"].GetDouble() 
        self.active_check_factor      = self.params["active_check_factor"].GetDouble() 
        self.max_number_results       = self.params["max_number_results"].GetInt() 
        self.normal_variation         = self.params["normal_variation"].GetBool()
        self.integration_order        = self.params["integration_order"].GetInt() 
        self.type_search = self.params["type_search"].GetString()
        
    def ExecuteInitialize(self):
        
        # Appending the conditions created to the self.main_model_part
        computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        if (computing_model_part.HasSubModelPart("Contact")):
            interface_model_part = computing_model_part.GetSubModelPart("Contact")
        else:
            interface_model_part = computing_model_part.CreateSubModelPart("Contact")

        if (self.normal_variation == True):
            computing_model_part.Set(KratosMultiphysics.INTERACTION, True)

        if (self.params["manual_ALM"].GetBool() == False):
            # Computing the scale factors or the penalty parameters (10 * E_mean/h_mean)
            find_nodal_h = KratosMultiphysics.FindNodalHProcess(computing_model_part)
            find_nodal_h.Execute()

            import statistics as stat
            nodal_h_values = []
            for node in computing_model_part.Nodes:
                nodal_h_values.append(node.GetSolutionStepValue(KratosMultiphysics.NODAL_H))
                
            mean_h = stat.mean(nodal_h_values)
                
            elem_E_values = []
            for elem in computing_model_part.Elements:
                prop = elem.Properties
                elem_E_values.append(prop[KratosMultiphysics.YOUNG_MODULUS])
                
            mean_E = stat.mean(elem_E_values)
            
            # Penalty and scalar factor
            penalty = 10.0 * mean_E/mean_h
            scale_factor = 10.0 * mean_E/mean_h
        else:
            # Penalty and scalar factor
            penalty = self.params["manual_ALM_parameters"]["penalty"].GetDouble()
            scale_factor = self.params["manual_ALM_parameters"]["scale_factor"].GetDouble()
        
        for prop in computing_model_part.GetProperties():
            prop[KratosMultiphysics.ContactStructuralMechanicsApplication.PENALTY_FACTOR] = penalty
            prop[KratosMultiphysics.ContactStructuralMechanicsApplication.SCALE_FACTOR] = scale_factor
            prop[KratosMultiphysics.ContactStructuralMechanicsApplication.INTEGRATION_ORDER_CONTACT] = self.integration_order
            
        for node in self.contact_model_part.Nodes:
            node.Set(KratosMultiphysics.INTERFACE, True)
        del(node)
        
        self.Preprocess = KratosMultiphysics.ContactStructuralMechanicsApplication.InterfacePreprocessCondition(computing_model_part)
        
        if self.params["contact_type"].GetString() == "Frictionless":
            condition_name = "ALMFrictionlessMortarContact"
        
        #print("MODEL PART BEFORE CREATING INTERFACE")
        #print(computing_model_part) 
        
        # It should create the conditions automatically
        if (self.dimension == 2):
            self.Preprocess.GenerateInterfacePart2D(computing_model_part, self.contact_model_part, condition_name, "", False) 
        else:
            self.Preprocess.GenerateInterfacePart3D(computing_model_part, self.contact_model_part, condition_name, "", False) 

        # When all conditions are simultaneously master and slave
        if (self.params["assume_master_slave"].GetString() == ""):
            for cond in self.contact_model_part.Conditions:
                cond.Set(KratosMultiphysics.SLAVE, True)
            del(cond)
            
        #print("MODEL PART AFTER CREATING INTERFACE")
        #print(computing_model_part)

        # We copy the conditions to the ContactSubModelPart
        for cond in self.contact_model_part.Conditions:
            interface_model_part.AddCondition(cond)    
        del(cond)
        for node in self.contact_model_part.Nodes:
            interface_model_part.AddNode(node, 0)    
        del(node)

        self.contact_search = KratosMultiphysics.ContactStructuralMechanicsApplication.TreeContactSearch(computing_model_part, self.max_number_results, self.active_check_factor, self.type_search)
        
        if self.params["contact_type"].GetString() == "Frictionless":
            ZeroVector = KratosMultiphysics.Vector(3) 
            ZeroVector[0] = 0.0
            ZeroVector[1] = 0.0
            ZeroVector[2] = 0.0
            
            # Initilialize weighted variables and LM
            for node in self.contact_model_part.Nodes:
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_GAP, 0.0)
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.AUXILIAR_ACTIVE, False)
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.AUXILIAR_SLIP,   False)
                node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)
                node.SetValue(KratosMultiphysics.NORMAL,      ZeroVector)
                node.SetValue(KratosMultiphysics.TANGENT_XI,  ZeroVector)
                node.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector)
            del(node)
            
            # Setting the conditions 
            for cond in self.contact_model_part.Conditions:
                cond.SetValue(KratosMultiphysics.NORMAL,      ZeroVector) 
                cond.SetValue(KratosMultiphysics.TANGENT_XI,  ZeroVector) 
                cond.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector) 
            del(cond)
            
            self.contact_search.CreatePointListMortar()
            self.contact_search.InitializeMortarConditions()
        
    def ExecuteBeforeSolutionLoop(self):
        if self.params["contact_type"].GetString() == "Frictionless":  
            self.contact_search.TotalClearALMFrictionlessMortarConditions()
    
    def ExecuteInitializeSolutionStep(self):
        if self.params["contact_type"].GetString() == "Frictionless":    
            self.contact_search.UpdateMortarConditions(self.search_factor)
            #self.contact_search.CheckMortarConditions()
        
    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        if self.params["contact_type"].GetString() == "Frictionless":
            self.contact_search.UpdatePointListMortar()
            self.contact_search.PartialClearALMFrictionlessMortarConditions()
            
    def ExecuteFinalize(self):
        pass
