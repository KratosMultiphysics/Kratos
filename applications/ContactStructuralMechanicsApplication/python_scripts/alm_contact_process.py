from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library 
import KratosMultiphysics 
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ALMContactProcess(Model, settings["Parameters"])

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
            "frictional_law"              : "Coulomb",
            "search_factor"               : 2.0,
            "active_check_factor"         : 0.01,
            "max_number_results"          : 1000,
            "bucket_size"                 : 4,
            "normal_variation"            : false,
            "manual_ALM"                  : false,
            "penalty"                     : 0.0,
            "scale_factor"                : 1.0,
            "type_search"                 : "InRadius",
            "integration_order"           : 2,
            "debug_mode"                  : false
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
            model_part_slave = self.main_model_part.GetSubModelPart(self.params["assume_master_slave"].GetString())
            for node in model_part_slave.Nodes:
                node.Set(KratosMultiphysics.SLAVE, True)
            del(node)
        
        self.normal_variation         = self.params["normal_variation"].GetBool()
        self.integration_order        = self.params["integration_order"].GetInt() 
        self.frictional_law = self.params["frictional_law"].GetString()
        self.debug_mode = self.params["debug_mode"].GetBool()
        
        # Debug
        if (self.debug_mode == True):
            self.label = 0
            self.output_file = "POSTSEARCH"

            self.gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
            self.singlefile = KratosMultiphysics.MultiFileFlag.SingleFile
            self.deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteDeformed
            self.write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly
        
    def ExecuteInitialize(self):
        # Appending the conditions created to the self.main_model_part
        computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        if (computing_model_part.HasSubModelPart("Contact")):
            interface_model_part = computing_model_part.GetSubModelPart("Contact")
        else:
            interface_model_part = computing_model_part.CreateSubModelPart("Contact")

        # We consider frictional contact (We use the SLIP flag because was the easiest way)
        if self.params["contact_type"].GetString() == "Frictional":
            computing_model_part.Set(KratosMultiphysics.SLIP, True) 
        else:
            computing_model_part.Set(KratosMultiphysics.SLIP, False) 
            
        # We recompute the normal at each iteration (false by default)
        if (self.normal_variation == True):
            computing_model_part.Set(KratosMultiphysics.INTERACTION, True) # TODO: Pending to implement
        else:
            computing_model_part.Set(KratosMultiphysics.INTERACTION, False)
        
        # Copying the properties in the contact model part
        self.contact_model_part.SetProperties(computing_model_part.GetProperties())
        
        # Setting the integration order
        for prop in computing_model_part.GetProperties():
            prop[ContactStructuralMechanicsApplication.INTEGRATION_ORDER_CONTACT] = self.integration_order
            
        for node in self.contact_model_part.Nodes:
            node.Set(KratosMultiphysics.INTERFACE, True)
        del(node)
        
        self.Preprocess = ContactStructuralMechanicsApplication.InterfacePreprocessCondition(computing_model_part)
        
        if self.params["contact_type"].GetString() == "Frictionless":
            condition_name = "ALMFrictionlessMortarContact"
        elif self.params["contact_type"].GetString() == "Frictional":
            condition_name = "ALMFrictionalMortarContact"
        
        #print("MODEL PART BEFORE CREATING INTERFACE")
        #print(computing_model_part) 
        
        # It should create the conditions automatically
        interface_parameters = KratosMultiphysics.Parameters("""{"condition_name": "", "final_string": "", "simplify_geometry": false}""")
        interface_parameters["condition_name"].SetString(condition_name)
        if (self.dimension == 2):
            self.Preprocess.GenerateInterfacePart2D(computing_model_part, self.contact_model_part, interface_parameters) 
        else:
            self.Preprocess.GenerateInterfacePart3D(computing_model_part, self.contact_model_part, interface_parameters) 

        # When all conditions are simultaneously master and slave
        if (self.params["assume_master_slave"].GetString() == ""):
            for node in self.contact_model_part.Conditions:
                node.Set(KratosMultiphysics.SLAVE, True)
            del(node)
            for cond in self.contact_model_part.Conditions:
                cond.Set(KratosMultiphysics.SLAVE, True)
            del(cond)

        if (self.params["manual_ALM"].GetBool() == False):
            # Computing the scale factors or the penalty parameters (10 * E_mean/h_mean)
            self.find_nodal_h = KratosMultiphysics.FindNodalHProcess(computing_model_part)
            self.find_nodal_h.Execute()
            
            self.alm_var_process = ContactStructuralMechanicsApplication.ALMVariablesCalculationProcess(self.contact_model_part)
            self.alm_var_process.Execute()
        else:
            # Penalty and scalar factor
            penalty = self.params["penalty"].GetDouble()
            scale_factor = self.params["scale_factor"].GetDouble()
            
            # We set the values in the process info
            self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.PENALTY_FACTOR] = penalty
            self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.SCALE_FACTOR]   = scale_factor
            
        #print("MODEL PART AFTER CREATING INTERFACE")
        #print(computing_model_part)

        # We copy the conditions to the ContactSubModelPart
        for cond in self.contact_model_part.Conditions:
            interface_model_part.AddCondition(cond)    
        del(cond)
        for node in self.contact_model_part.Nodes:
            interface_model_part.AddNode(node, 0)    
        del(node)

        # Creating the search
        search_parameters = KratosMultiphysics.Parameters("""{}""")
        search_parameters.AddValue("type_search",self.params["type_search"])
        search_parameters.AddValue("allocation_size",self.params["max_number_results"])
        search_parameters.AddValue("bucket_size",self.params["bucket_size"])
        search_parameters.AddValue("search_factor",self.params["search_factor"])
        search_parameters.AddValue("active_check_factor",self.params["active_check_factor"])
        self.contact_search = ContactStructuralMechanicsApplication.TreeContactSearch(computing_model_part, search_parameters)
        
        # We initialize the conditions    
        self.alm_init_var = ContactStructuralMechanicsApplication.ALMFastInit(self.contact_model_part) 
        self.alm_init_var.Execute()
        
        self.contact_search.CreatePointListMortar()
        self.contact_search.InitializeMortarConditions()
        
    def ExecuteBeforeSolutionLoop(self):
        if self.params["contact_type"].GetString() == "Frictionless":  
            self.contact_search.TotalClearALMFrictionlessMortarConditions()
        else:
            self.contact_search.TotalClearComponentsMortarConditions()
    
    def ExecuteInitializeSolutionStep(self):
        self.contact_search.UpdateMortarConditions()
        #self.contact_search.CheckMortarConditions()
            
        # Debug
        if (self.debug_mode == True):
            self.label += 1
            
            gid_io = KratosMultiphysics.GidIO(self.output_file+"_STEP_"+str(self.label), self.gid_mode, self.singlefile, self.deformed_mesh_flag, self.write_conditions)
            
            gid_io.InitializeMesh(self.label)
            gid_io.WriteMesh(self.main_model_part.GetMesh())
            gid_io.FinalizeMesh()
            gid_io.InitializeResults(self.label, self.main_model_part.GetMesh())
            
            gid_io.WriteNodalFlags(KratosMultiphysics.ACTIVE, "ACTIVE", self.main_model_part.Nodes, self.label)
            gid_io.WriteNodalResults(KratosMultiphysics.DISPLACEMENT, self.main_model_part.Nodes, self.label, 0)
            
            gid_io.FinalizeResults()
            
            #raise NameError("DEBUG")
        
    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        if self.params["contact_type"].GetString() == "Frictionless":
            self.contact_search.PartialClearALMFrictionlessMortarConditions()
        else:
            self.contact_search.PartialClearComponentsMortarConditions()
            
    def ExecuteFinalize(self):
        pass
            
