from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.PfemApplication as KratosPfem
KratosMultiphysics.CheckForPreviousImport()

#from multiprocessing import Pool

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return RemeshDomainsProcess(Model, settings["Parameters"])


class RemeshDomainsProcess(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)
        
        self.main_model_part = Model[custom_settings["model_part_name"].GetString()]
    
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "echo_level"            : 1,
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
        self.echo_level        = self.settings["echo_level"].GetInt()

        self.dimension         = self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]
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
            domain_module = __import__(item["python_module"].GetString())
            domain = domain_module.CreateMeshingDomain(self.main_model_part,item)
            self.meshing_domains.append(domain)

        # mesh modeler initial values
        self.remesh_domains_active = False
        for domain in self.meshing_domains:
            if( domain.Active() ):
                self.remesh_domains_active = True

        self.step_count   = 1
        self.counter      = 1
        self.next_meshing = 0.0
        self.meshing_before_output = self.settings["meshing_before_output"].GetBool()
     
    #
    def ExecuteInitialize(self):

        # check restart
        self.restart = False
        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
            self.restart = True         
            self.step_count = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
            
            if self.meshing_control_is_time:
                self.next_meshing  = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] + self.meshing_frequency
            else:
                self.next_meshing = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + self.meshing_frequency
        else:
            #if commented in the first step meshing is applied
            self.next_meshing = self.meshing_frequency
            # it must be initialized if restart is called//no

        self.main_model_part.ProcessInfo.SetValue(KratosPfem.INITIALIZED_DOMAINS, False);
        
        # initialize all meshing domains 
        if( self.remesh_domains_active ):    

            print("::[Meshing_Process]:: Initialize Domains")
            import domain_utilities
            domain_utils = domain_utilities.DomainUtilities()            
            domain_utils.InitializeDomains(self.main_model_part,self.echo_level)

            for domain in self.meshing_domains:
                domain.SetEchoLevel(self.echo_level)
                domain.Initialize()
                #domain.Check()
                         

    ###

    #
    def ExecuteInitializeSolutionStep(self):

        self.step_count += 1

    #
    def ExecuteBeforeOutputStep(self):
        
        if(self.remesh_domains_active):
            if( self.meshing_before_output ):
                self.main_model_part.ProcessInfo[KratosPfem.MESHING_STEP_PERFORMED] = False
                if(self.IsMeshingStep()):
                    self.RemeshDomains()
 
    #
    def ExecuteAfterOutputStep(self):

        if(self.remesh_domains_active):
            if( not self.meshing_before_output ):
                self.main_model_part.ProcessInfo[KratosPfem.MESHING_STEP_PERFORMED] = False
                if(self.IsMeshingStep()):
                    self.RemeshDomains()
                    
    ###

    #
    def ExecuteMeshing(domain):
        domain.ExecuteMeshing()

    #
    def RemeshDomains(self):

        print("")
        print("::[Meshing_Process]:: MESHING DOMAIN...( call:", self.counter,")")
            
        meshing_options = KratosMultiphysics.Flags()
        self.model_meshing = KratosPfem.ModelMeshing(self.main_model_part, meshing_options, self.echo_level)
        
        self.model_meshing.ExecuteInitialize()

        #serial
        for domain in self.meshing_domains:
            domain.ExecuteMeshing()
        
        
        #parallel (not working pickling instances not enabled)
        #domains_number = len(self.meshing_domains)
        #if(domains_number>8):
        #    domains_number = 8
        
        #pool = Pool(domains_number)
        #pool.map(self.ExecuteMeshing,self.meshing_domains)
        #pool.close()
        #pool.joint()        
        #
        
        self.model_meshing.ExecuteFinalize()
        
        self.counter += 1 
        
        self.main_model_part.ProcessInfo[KratosPfem.MESHING_STEP_PERFORMED] = True
        
        # schedule next meshing
        if(self.meshing_frequency > 0.0): # note: if == 0 always active
            if(self.meshing_control_is_time):
                time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
                while(self.next_meshing <= time):
                    self.next_meshing += self.meshing_frequency
            else:
                while(self.next_meshing <= self.step_count):
                    self.next_meshing += self.meshing_frequency

    #
    def GetMeshingStep(self):
        return self.counter

    #
    def IsMeshingStep(self):

        if(self.meshing_control_is_time):
            #print( str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME])+">"+ str(self.next_meshing) )
            return ( self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] >= self.next_meshing )
        else:
            return ( self.step_count >= self.next_meshing )
