 
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library 
import KratosMultiphysics 
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication

KratosMultiphysics.CheckForPreviousImport()

# TODO: Finish me!!!!!

def CalculateLastIdCondition(model_part):
    cond_id = 0
    for condition in model_part.Conditions:
        cond_id += 1

    return cond_id

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
            "contact_type"                : "Frictionless",
            "search_factor"               : 1.5,
            "active_check_factor"         : 0.01,
            "max_number_results"          : 1000,
            "simplify_geometry"           : false,
            "normal_variation"            : false,
            "type_search"                 : "InRadius",
            "integration_type"            : "Collocation",
            "integration_order"           : 5
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
        self.simplify_geometry        = self.params["simplify_geometry"].GetBool()
        self.normal_variation         = self.params["normal_variation"].GetBool()
        self.integration_type         = self.params["integration_type"].GetString() 
        self.integration_order        = self.params["integration_order"].GetInt() 
        if self.params["type_search"].GetString() == "InRadius":
             self.type_search = 0
        
    def ExecuteInitialize(self):
        
        # Appending the conditions created to the computing_model_part
        computing_model_part = self.main_model_part.GetSubModelPart("computing_domain")
        computing_model_part.CreateSubModelPart("Contact")
        interface_computing_model_part = computing_model_part.GetSubModelPart("Contact")

        if (self.normal_variation == True):
            computing_model_part.Set(KratosMultiphysics.INTERACTION, True)

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
            
        penalty = 10.0 * mean_E/mean_h
        scale_factor = 10.0 * mean_E/mean_h
        
        for prop in computing_model_part.GetProperties():
            prop[KratosMultiphysics.ContactStructuralMechanicsApplication.PENALTY_FACTOR] = penalty
            prop[KratosMultiphysics.ContactStructuralMechanicsApplication.SCALE_FACTOR] = scale_factor
            
        for node in self.o_interface.Nodes:
            interface_computing_model_part.AddNode(node, 0)  
            node.Set(KratosMultiphysics.INTERFACE,True)
        del(node)
        
        for node in self.d_interface.Nodes:
            interface_computing_model_part.AddNode(node, 0)
            node.Set(KratosMultiphysics.INTERFACE,True)
        del(node)
        
        self.Preprocess  = KratosMultiphysics.ContactStructuralMechanicsApplication.InterfacePreprocessCondition()
        
        if self.params["contact_type"].GetString() == "Frictionless":
            condition_name = "AugmentedLagrangianMethodFrictionlessMortar"
        
        #print("MODEL PART BEFORE CREATING INTERFACE")
        #print(self.main_model_part) 
        
        # It should create the conditions automatically
        initial_id = CalculateLastIdCondition(self.main_model_part)
        self.Preprocess.GenerateInterfacePart(self.o_model_part, self.o_interface, condition_name, initial_id, "", self.simplify_geometry) 
        initial_id = CalculateLastIdCondition(self.main_model_part)
        self.Preprocess.GenerateInterfacePart(self.d_model_part, self.d_interface, condition_name, initial_id, "", self.simplify_geometry) 

        #print("MODEL PART AFTER CREATING INTERFACE")
        #print(self.main_model_part)
        
        for cond in self.o_interface.Conditions:
            interface_computing_model_part.AddCondition(cond)    
        del(cond)
        
        for cond in self.d_interface.Conditions:
            interface_computing_model_part.AddCondition(cond)  
        del(cond)

        self.contact_search = KratosMultiphysics.ContactStructuralMechanicsApplication.TreeContactSearch(self.o_interface, self.d_interface, self.max_number_results)
        
        if self.params["contact_type"].GetString() == "Frictionless":
            ZeroVector = KratosMultiphysics.Vector(3) 
            ZeroVector[0] = 0.0
            ZeroVector[1] = 0.0
            ZeroVector[2] = 0.0
            
            # Initilialize weighted variables and LM
            for node in self.d_interface.Nodes:
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_GAP, 0.0)
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.AUXILIAR_ACTIVE, False)
                node.SetValue(KratosMultiphysics.ContactStructuralMechanicsApplication.AUXILIAR_SLIP, False)
                node.SetValue(KratosMultiphysics.NODAL_AREA, 0.0)
                node.SetValue(KratosMultiphysics.NORMAL, ZeroVector)
                node.SetValue(KratosMultiphysics.TANGENT_XI, ZeroVector)
                node.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector)
                #node.Set(KratosMultiphysics.SLAVE, True)
            del node
            
            # Setting the master conditions 
            for cond in self.o_interface.Nodes:
                cond.SetValue(KratosMultiphysics.NORMAL, ZeroVector) 
                cond.SetValue(KratosMultiphysics.TANGENT_XI, ZeroVector) 
                cond.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector) 
                #cond.Set(KratosMultiphysics.MASTER, True) # TODO: This is not supposed o be necessary
            del cond
            
            # Setting the slave conditions 
            for cond in self.d_interface.Nodes:
                cond.SetValue(KratosMultiphysics.NORMAL, ZeroVector) 
                cond.SetValue(KratosMultiphysics.TANGENT_XI, ZeroVector) 
                cond.SetValue(KratosMultiphysics.TANGENT_ETA, ZeroVector) 
            del cond
            
            self.contact_search.CreatePointListMortar()
            self.contact_search.InitializeALMFrictionlessMortarConditions(self.active_check_factor, self.integration_order)
        
    def ExecuteBeforeSolutionLoop(self):
        if self.params["contact_type"].GetString() == "Frictionless":  
            self.contact_search.TotalClearALMFrictionlessMortarConditions();
    
    def ExecuteInitializeSolutionStep(self):
        #for cond in self.d_interface.Conditions:
            #print(cond.Is(KratosMultiphysics.ACTIVE))
        
        if self.params["contact_type"].GetString() == "Frictionless":    
            self.contact_search.UpdateMortarConditions(self.search_factor, self.type_search)
            #self.contact_search.CheckMortarConditions()
            
        #for cond in self.d_interface.Conditions:
            #print(cond.Is(KratosMultiphysics.ACTIVE))
        
    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        if self.params["contact_type"].GetString() == "Frictionless":
            self.contact_search.UpdatePointListMortar()
            self.contact_search.PartialClearALMFrictionlessMortarConditions()
            
        #for cond in self.d_interface.Conditions:
            #print(cond.Is(KratosMultiphysics.ACTIVE))
            
    def ExecuteFinalize(self):
        pass
