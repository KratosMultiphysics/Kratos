from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
import KratosMultiphysics.ContactMechanicsApplication as KratosContact
KratosMultiphysics.CheckForPreviousImport()

import remesh_domains_process

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ContactDomainProcess(Model, settings["Parameters"])


class ContactDomainProcess(remesh_domains_process.RemeshDomainsProcess):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)
        
        self.model_part = Model[custom_settings["model_part_name"].GetString()]
    
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"       : "Solid Domain",
            "meshing_control_type"  : "step",
            "meshing_frequency"     : 1.0,
            "meshing_before_output" : true,
            "meshing_domains"       : []
        }
        """)
 
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.echo_level        = 1
        self.domain_size       = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        self.meshing_frequency = self.settings["meshing_frequency"].GetDouble()
        
        self.meshing_control_is_time = False
        meshing_control_type   = self.settings["meshing_control_type"].GetString()
        if(meshing_control_type == "time"):
            self.meshing_control_is_time = True
        elif(meshing_control_type == "step"):
            self.meshing_control_is_time = False

        #construct meshing domains
        self.meshing_domains = []
        domains_list = self.settings["meshing_domains"]
        self.number_of_domains = domains_list.size()
        for i in range(0,self.number_of_domains):
            item = domains_list[i]
            domain_module = __import__(item["python_file_name"].GetString())
            domain = domain_module.CreateMeshingDomain(self.model_part,item)
            self.meshing_domains.append(domain)

        # mesh modeler initial values
        self.remesh_domains_active = False
        self.neighbours_search_performed = False 

        self.step_count   = 1
        self.counter      = 1
        self.next_meshing = 0.0
        self.meshing_before_output = self.settings["meshing_before_output"].GetBool()
                       
    #
    def ExecuteInitialize(self):

        # Execute initialize base class
        if( self.model_part.ProcessInfo[MESHING_ACTIVE] == False ):
            remesh_domains_process.RemesDomainsProcess.ExecuteInitialize(self)

    ###

    #
    def RemeshDomains(self):

        if( self.echo_level > 0 ):
            print("::[Meshing_Process]:: CONTACT SEARCH", self.counter)
            
        meshing_options = KratosMultiphysics.Flags()
        self.model_meshing = KratosContact.ContactModelMeshing(self.model_part, meshing_options, self.echo_level)
        
        self.model_meshing.ExecuteInitialize()

        for domain in self.meshing_domains:
            domain.ExecuteMeshing();
 
        self.model_meshing.ExecuteFinalize()
        
        self.counter += 1 


        # schedule next meshing
        if(self.meshing_frequency >= 0.0):
            if(self.meshing_control_is_time):
                while(self.next_meshing <= time):
                    self.next_meshing += self.meshing_frequency
            else:
                while(self.next_meshing <= self.step_count):
                    self.next_meshing += self.meshing_frequency
                        



