from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library 
import KratosMultiphysics 
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ALMContactProcess(Model, settings["Parameters"])

import python_process
##all the processes python processes should be derived from "python_process"
class ALMContactProcess(python_process.PythonProcess):
  
    def __init__(self,model_part,params):
        
        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "mesh_id"                     : 0,
            "model_part_name"             : "Structure",
            "computing_model_part_name"   : "computing_domain",
            "contact_model_part"          : "Contact_Part",
            "axisymmetric"                : false,
            "assume_master_slave"         : "",
            "contact_type"                : "Frictionless",
            "frictional_law"              : "Coulomb",
            "search_factor"               : 2.0,
            "active_check_factor"         : 0.01,
            "max_number_results"          : 1000,
            "bucket_size"                 : 4,
            "normal_variation"            : "NO_DERIVATIVES_COMPUTATION",
            "manual_ALM"                  : false,
            "stiffness_factor"            : 1.0,
            "penalty_scale_factor"        : 1.0,
            "use_scale_factor"            : true,
            "penalty"                     : 0.0,
            "scale_factor"                : 1.0e0,
            "tangent_factor"              : 0.1,
            "type_search"                 : "InRadius",
            "check_gap"                   : "CheckMapping",
            "database_step_update"        : 1,
            "integration_order"           : 3,
            "max_gap_factor"              : 1.0e-3,
            "debug_mode"                  : false,
            "remeshing_with_contact_bc"   : false
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
   
        self.main_model_part = model_part[self.params["model_part_name"].GetString()]
        self.computing_model_part_name = self.params["computing_model_part_name"].GetString()

        self.dimension = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        self.contact_model_part = model_part[self.params["contact_model_part"].GetString()]
        
        self.axisymmetric  = self.params["axisymmetric"].GetBool()
        if (self.axisymmetric == True) and (self.dimension == 3):
            raise NameError("3D and axisymmetric makes no sense")
        if (self.params["normal_variation"].GetString() == "NO_DERIVATIVES_COMPUTATION"):
            self.normal_variation = 0
        elif (self.params["normal_variation"].GetString() == "ELEMENTAL_DERIVATIVES"):
            self.normal_variation = 1
        elif (self.params["normal_variation"].GetString() == "NODAL_ELEMENTAL_DERIVATIVES"):
            self.normal_variation = 2
        else:
            raise NameError("The options to normal derivatives are: NO_DERIVATIVES_COMPUTATION, ELEMENTAL_DERIVATIVES, NODAL_ELEMENTAL_DERIVATIVES")
        self.database_step_update = self.params["database_step_update"].GetInt() 
        self.database_step = 0
        self.frictional_law = self.params["frictional_law"].GetString()
        self.debug_mode = self.params["debug_mode"].GetBool()
        
        # Debug
        if (self.debug_mode == True):
            self.output_file = "POSTSEARCH"

            self.gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
            self.singlefile = KratosMultiphysics.MultiFileFlag.SingleFile
            self.deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
            self.write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly
        
    def ExecuteInitialize(self):
        # The computing model part
        computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        
        # We compute NODAL_H that can be used in the search and some values computation
        self.find_nodal_h = KratosMultiphysics.FindNodalHProcess(computing_model_part)
        self.find_nodal_h.Execute()
        
        # Assigning master and slave sides
        self._assign_slave_nodes()
        
        # Appending the conditions created to the self.main_model_part
        if (computing_model_part.HasSubModelPart("Contact")):
            preprocess = False
            interface_model_part = computing_model_part.GetSubModelPart("Contact")
        else:
            preprocess = True
            interface_model_part = computing_model_part.CreateSubModelPart("Contact")

        # We consider frictional contact (We use the SLIP flag because was the easiest way)
        if self.params["contact_type"].GetString() == "Frictional":
            computing_model_part.Set(KratosMultiphysics.SLIP, True) 
        else:
            computing_model_part.Set(KratosMultiphysics.SLIP, False) 
            
        # We recompute the normal at each iteration (false by default)
        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.CONSIDER_NORMAL_VARIATION] = self.normal_variation
        # We set the max gap factor for the gap adaptation
        max_gap_factor = self.params["max_gap_factor"].GetDouble()
        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.ADAPT_PENALTY] = (max_gap_factor > 0.0)
        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.MAX_GAP_FACTOR] = max_gap_factor
        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.ACTIVE_CHECK_FACTOR] = self.params["active_check_factor"].GetDouble()
        
        # We set the value that scales in the tangent direction the penalty and scale parameter
        if self.params["contact_type"].GetString() == "Frictional":
            self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.TANGENT_FACTOR] = self.params["tangent_factor"].GetDouble()
        
        # Copying the properties in the contact model part
        self.contact_model_part.SetProperties(computing_model_part.GetProperties())
        
        # Setting the integration order and active check factor
        for prop in computing_model_part.GetProperties():
            prop[ContactStructuralMechanicsApplication.INTEGRATION_ORDER_CONTACT] = self.params["integration_order"].GetInt()
            
        for node in self.contact_model_part.Nodes:
            node.Set(KratosMultiphysics.INTERFACE, True)
        del(node)
        
        #If the conditions doesn't exist we create them
        if (preprocess == True):
            self._interface_preprocess(computing_model_part)
        else:
            master_slave_process = ContactStructuralMechanicsApplication.MasterSlaveProcess(computing_model_part) 
            master_slave_process.Execute()
        
        # We initialize the contact values
        self._initialize_contact_values(computing_model_part)

        # When all conditions are simultaneously master and slave
        self._assign_slave_conditions()

        # We initialize the ALM parameters
        self._initialize_alm_parameters(computing_model_part)

        # We copy the conditions to the ContactSubModelPart
        if (preprocess == True):
            for cond in self.contact_model_part.Conditions:
                interface_model_part.AddCondition(cond)    
            del(cond)
            for node in self.contact_model_part.Nodes:
                interface_model_part.AddNode(node, 0)   
            del(node)

        # Creating the search
        self._create_main_search(computing_model_part)
        
        # We initialize the conditions    
        alm_init_var = ContactStructuralMechanicsApplication.ALMFastInit(self.contact_model_part) 
        alm_init_var.Execute()
        
        # We initialize the search utility
        self.contact_search.CreatePointListMortar()
        self.contact_search.InitializeMortarConditions()
        
    def ExecuteBeforeSolutionLoop(self):
        pass
    
    def ExecuteInitializeSolutionStep(self):
        self.database_step += 1
        self.global_step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        
        if (self.database_step >= self.database_step_update or self.global_step == 1):
            # We solve one linear step with a linear strategy if needed
            # Clear current pairs
            self._clear_sets(self.contact_search)
            # Update database
            self.contact_search.UpdateMortarConditions()
            #self.contact_search.CheckMortarConditions()
                
            # Debug
            if (self.debug_mode == True):
               self._debug_output(self.global_step, "")
        
    def ExecuteFinalizeSolutionStep(self):
        if (self.params["remeshing_with_contact_bc"].GetBool() == True):
            self._transfer_slave_to_master()

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        modified = self.main_model_part.Is(KratosMultiphysics.MODIFIED)
        if (modified == False and (self.database_step >= self.database_step_update or self.global_step == 1)):
            self._clear_sets(self.contact_search)
            self.database_step = 0
            
    def ExecuteFinalize(self):
        pass

    def _clear_sets(self, contact_search):
        if self.params["contact_type"].GetString() == "Frictionless":  
            contact_search.ClearALMFrictionlessMortarConditions()
        else:
            contact_search.ClearComponentsMortarConditions()

    def _assign_slave_conditions(self):
        if (self.params["assume_master_slave"].GetString() == ""):
            for cond in self.contact_model_part.Conditions:
                cond.Set(KratosMultiphysics.SLAVE, True)
            del(cond)
            
    def _assign_slave_nodes(self):
        if (self.params["assume_master_slave"].GetString() != ""):
            for node in self.contact_model_part.Nodes:
                node.Set(KratosMultiphysics.SLAVE, False)
                node.Set(KratosMultiphysics.MASTER, True)
            del(node)
            model_part_slave = self.main_model_part.GetSubModelPart(self.params["assume_master_slave"].GetString())
            for node in model_part_slave.Nodes:
                node.Set(KratosMultiphysics.SLAVE, True)
                node.Set(KratosMultiphysics.MASTER, False)
            del(node)
            
    def _interface_preprocess(self, computing_model_part):
        self.interface_preprocess = ContactStructuralMechanicsApplication.InterfacePreprocessCondition(computing_model_part)
                    
        #print("MODEL PART BEFORE CREATING INTERFACE")
        #print(computing_model_part) 
        
        # It should create the conditions automatically
        interface_parameters = KratosMultiphysics.Parameters("""{"simplify_geometry": false}""")
        if (self.dimension == 2):
            self.interface_preprocess.GenerateInterfacePart2D(computing_model_part, self.contact_model_part, interface_parameters) 
        else:
            self.interface_preprocess.GenerateInterfacePart3D(computing_model_part, self.contact_model_part, interface_parameters) 
            
        #print("MODEL PART AFTER CREATING INTERFACE")
        #print(computing_model_part)
            
    def _initialize_contact_values(self, computing_model_part):
        # We consider frictional contact (We use the SLIP flag because was the easiest way)
        if self.params["contact_type"].GetString() == "Frictional":
            computing_model_part.Set(KratosMultiphysics.SLIP, True) 
        else:
            computing_model_part.Set(KratosMultiphysics.SLIP, False) 
            
        # We recompute the normal at each iteration (false by default)
        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.CONSIDER_NORMAL_VARIATION] = self.normal_variation
        # We set the max gap factor for the gap adaptation
        max_gap_factor = self.params["max_gap_factor"].GetDouble()
        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.ADAPT_PENALTY] = (max_gap_factor > 0.0)
        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.MAX_GAP_FACTOR] = max_gap_factor
        
        # We set the value that scales in the tangent direction the penalty and scale parameter
        if self.params["contact_type"].GetString() == "Frictional":
            self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.TANGENT_FACTOR] = self.params["tangent_factor"].GetDouble()
        
        # Copying the properties in the contact model part
        self.contact_model_part.SetProperties(computing_model_part.GetProperties())
        
        # Setting the integration order and active check factor
        for prop in computing_model_part.GetProperties():
            prop[ContactStructuralMechanicsApplication.INTEGRATION_ORDER_CONTACT] = self.params["integration_order"].GetInt() 
            prop[ContactStructuralMechanicsApplication.ACTIVE_CHECK_FACTOR] = self.params["active_check_factor"].GetDouble()
            
    def _initialize_alm_parameters(self, computing_model_part):
        if (self.params["manual_ALM"].GetBool() == False):
            # Computing the scale factors or the penalty parameters (StiffnessFactor * E_mean/h_mean)
            alm_var_parameters = KratosMultiphysics.Parameters("""{}""")
            alm_var_parameters.AddValue("stiffness_factor",self.params["stiffness_factor"])
            alm_var_parameters.AddValue("penalty_scale_factor",self.params["penalty_scale_factor"])
            self.alm_var_process = ContactStructuralMechanicsApplication.ALMVariablesCalculationProcess(self.contact_model_part,KratosMultiphysics.NODAL_H, alm_var_parameters)
            self.alm_var_process.Execute()
            # We don't consider scale factor
            if (self.params["use_scale_factor"].GetBool() == False):
                self.main_model_part.ProcessInfo[KratosMultiphysics.SCALE_FACTOR] = 1.0
        else:
            # We set the values in the process info
            self.main_model_part.ProcessInfo[KratosMultiphysics.INITIAL_PENALTY] = self.params["penalty"].GetDouble()
            self.main_model_part.ProcessInfo[KratosMultiphysics.SCALE_FACTOR] = self.params["scale_factor"].GetDouble()
            
        # We print the parameters considered
        print("The parameters considered finally are: ")            
        print("SCALE_FACTOR: ","{:.2e}".format(self.main_model_part.ProcessInfo[KratosMultiphysics.SCALE_FACTOR]))
        print("INITIAL_PENALTY: ","{:.2e}".format(self.main_model_part.ProcessInfo[KratosMultiphysics.INITIAL_PENALTY]))
            
    def _create_main_search(self, computing_model_part):
        if self.params["contact_type"].GetString() == "Frictionless":
            if self.normal_variation == 2:
                if self.axisymmetric == True:
                    condition_name = "ALMNVFrictionlessAxisymMortarContact"
                else:
                    condition_name = "ALMNVFrictionlessMortarContact"
            else:
                if self.axisymmetric == True:
                    condition_name = "ALMFrictionlessAxisymMortarContact"
                else:
                    condition_name = "ALMFrictionlessMortarContact"
        elif self.params["contact_type"].GetString() == "Frictional":
            if self.normal_variation == 2:
                if self.axisymmetric == True:
                    condition_name = "ALMNVFrictionalAxisymMortarContact"
                else:
                    condition_name = "ALMNVFrictionalMortarContact"
            else:
                if self.axisymmetric == True:
                    condition_name = "ALMFrictionalAxisymMortarContact"
                else:
                    condition_name = "ALMFrictionalMortarContact"
        search_parameters = KratosMultiphysics.Parameters("""{"condition_name": "", "final_string": ""}""")
        search_parameters.AddValue("type_search",self.params["type_search"])
        search_parameters.AddValue("check_gap",self.params["check_gap"])
        search_parameters.AddValue("allocation_size",self.params["max_number_results"])
        search_parameters.AddValue("bucket_size",self.params["bucket_size"])
        search_parameters.AddValue("search_factor",self.params["search_factor"])
        search_parameters["condition_name"].SetString(condition_name)
        for cond in computing_model_part.Conditions:
            number_nodes = len(cond.GetNodes())
            break
        del(cond)
        if (self.dimension == 2):
            self.contact_search = ContactStructuralMechanicsApplication.TreeContactSearch2D2N(computing_model_part, search_parameters)
        else:
            if (number_nodes == 3):
                self.contact_search = ContactStructuralMechanicsApplication.TreeContactSearch3D3N(computing_model_part, search_parameters)
            else:
                self.contact_search = ContactStructuralMechanicsApplication.TreeContactSearch3D4N(computing_model_part, search_parameters)
    
    def _transfer_slave_to_master(self):
    
        for cond in self.contact_model_part.Conditions:
            break
        num_nodes = len(cond.GetNodes())
    
        # We use the search utility
        self._reset_search()
        self.contact_search.UpdateMortarConditions()
        #self.contact_search.CheckMortarConditions()
        
        map_parameters = KratosMultiphysics.Parameters("""
        {
            "echo_level"                       : 0,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2,
            "inverted_master_slave_pairing"    : true
        }
        """)
        
        computing_model_part = self.main_model_part.GetSubModelPart(self.computing_model_part_name)
        interface_model_part = computing_model_part.GetSubModelPart("Contact")
        if (self.dimension == 2): 
            #if self.params["contact_type"].GetString() == "Frictional":
                #mortar_mapping0 = KratosMultiphysics.SimpleMortarMapperProcess2D2NVectorHistorical(interface_model_part, KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER, map_parameters)
            #else:
                #mortar_mapping0 = KratosMultiphysics.SimpleMortarMapperProcess2D2NDoubleHistorical(interface_model_part, KratosMultiphysics.NORMAL_CONTACT_STRESS, map_parameters)
            mortar_mapping1 = KratosMultiphysics.SimpleMortarMapperProcess2D2NDoubleNonHistorical(interface_model_part, ContactStructuralMechanicsApplication.AUGMENTED_NORMAL_CONTACT_PRESSURE, map_parameters)
        else:
            if (num_nodes == 3): 
                #if self.params["contact_type"].GetString() == "Frictional":
                    #mortar_mapping0 = KratosMultiphysics.SimpleMortarMapperProcess3D3NVectorHistorical(interface_model_part, KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER, map_parameters)
                #else:
                    #mortar_mapping0 = KratosMultiphysics.SimpleMortarMapperProcess3D3NDoubleHistorical(interface_model_part, KratosMultiphysics.NORMAL_CONTACT_STRESS, map_parameters)
                mortar_mapping1 = KratosMultiphysics.SimpleMortarMapperProcess3D3NDoubleNonHistorical(interface_model_part, ContactStructuralMechanicsApplication.AUGMENTED_NORMAL_CONTACT_PRESSURE, map_parameters)
            else:
                #if self.params["contact_type"].GetString() == "Frictional":
                    #mortar_mapping0 = KratosMultiphysics.SimpleMortarMapperProcess3D4NVectorHistorical(interface_model_part, KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER, map_parameters)
                #else:
                    #mortar_mapping0 = KratosMultiphysics.SimpleMortarMapperProcess3D4NDoubleHistorical(interface_model_part, KratosMultiphysics.NORMAL_CONTACT_STRESS, map_parameters)
                mortar_mapping1 = KratosMultiphysics.SimpleMortarMapperProcess3D4NDoubleNonHistorical(interface_model_part, ContactStructuralMechanicsApplication.AUGMENTED_NORMAL_CONTACT_PRESSURE, map_parameters)
                    
        #mortar_mapping0.Execute()
        mortar_mapping1.Execute()
        
        # Transfering the AUGMENTED_NORMAL_CONTACT_PRESSURE to NORMAL_CONTACT_STRESS
        for node in interface_model_part.Nodes:
            aug_press = node.GetValue(ContactStructuralMechanicsApplication.AUGMENTED_NORMAL_CONTACT_PRESSURE)
            node.SetValue(KratosMultiphysics.NORMAL_CONTACT_STRESS, aug_press)
        del(node)
    
        self._reset_search()
    
    def _reset_search(self):
        self.contact_search.InvertSearch()
        self.contact_search.ResetContactOperators()
        self.contact_search.CreatePointListMortar()
        self.contact_search.InitializeMortarConditions()
        
    def _debug_output(self, label, name):

        gid_io = KratosMultiphysics.GidIO(self.output_file+name+"_STEP_"+str(label), self.gid_mode, self.singlefile, self.deformed_mesh_flag, self.write_conditions)
        
        gid_io.InitializeMesh(label)
        gid_io.WriteMesh(self.main_model_part.GetMesh())
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(label, self.main_model_part.GetMesh())
        
        gid_io.WriteNodalFlags(KratosMultiphysics.INTERFACE, "INTERFACE", self.main_model_part.Nodes, label)
        gid_io.WriteNodalFlags(KratosMultiphysics.ACTIVE, "ACTIVE", self.main_model_part.Nodes, label)
        gid_io.WriteNodalFlags(KratosMultiphysics.SLAVE, "SLAVE", self.main_model_part.Nodes, label)
        gid_io.WriteNodalResults(KratosMultiphysics.NORMAL, self.main_model_part.Nodes, label, 0)
        gid_io.WriteNodalResultsNonHistorical(ContactStructuralMechanicsApplication.AUGMENTED_NORMAL_CONTACT_PRESSURE, self.main_model_part.Nodes, label)
        gid_io.WriteNodalResultsNonHistorical(KratosMultiphysics.NODAL_AREA, self.main_model_part.Nodes, label)
        gid_io.WriteNodalResults(KratosMultiphysics.DISPLACEMENT, self.main_model_part.Nodes, label, 0)
        dynamic = False
        for node in self.main_model_part.Nodes:
            if (node.SolutionStepsDataHas(KratosMultiphysics.VELOCITY_X) == True):
                dynamic = True
            break
        del(node)
        if (dynamic == True):
            gid_io.WriteNodalResults(KratosMultiphysics.VELOCITY, self.main_model_part.Nodes, label, 0)
            gid_io.WriteNodalResults(KratosMultiphysics.ACCELERATION, self.main_model_part.Nodes, label, 0)
        gid_io.WriteNodalResults(KratosMultiphysics.NORMAL_CONTACT_STRESS, self.main_model_part.Nodes, label, 0)
        gid_io.WriteNodalResults(ContactStructuralMechanicsApplication.WEIGHTED_GAP, self.main_model_part.Nodes, label, 0)
        gid_io.WriteNodalResultsNonHistorical(ContactStructuralMechanicsApplication.NORMAL_GAP, self.main_model_part.Nodes, label)
        gid_io.WriteNodalResultsNonHistorical(ContactStructuralMechanicsApplication.AUXILIAR_COORDINATES, self.main_model_part.Nodes, label)
        gid_io.WriteNodalResultsNonHistorical(ContactStructuralMechanicsApplication.DELTA_COORDINATES, self.main_model_part.Nodes, label)
        
        gid_io.FinalizeResults()
        
        #raise NameError("DEBUG")
