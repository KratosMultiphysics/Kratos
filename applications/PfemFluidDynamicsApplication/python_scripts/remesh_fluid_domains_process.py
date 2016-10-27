from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
KratosMultiphysics.CheckForPreviousImport()

from multiprocessing import Pool

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return RemeshFluidDomainsProcess(Model, settings["Parameters"])


class RemeshFluidDomainsProcess(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)
        
        self.main_model_part = Model[custom_settings["model_part_name"].GetString()]
    
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"       : "Fluid Domain",
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
        self.domain_size       = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
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

        self.neighbours_search_performed = False
        self.step_count   = 1
        self.counter      = 1
        self.next_meshing = 0.0
        self.meshing_before_output = self.settings["meshing_before_output"].GetBool()
     
    #
    def ExecuteInitialize(self):

        print("::[Remesh_Fluid_Domains]:: Execute Initialize ")

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
            self.meshing_output = self.meshing_frequency


        self.main_model_part.ProcessInfo.SetValue(KratosPfemBase.INITIALIZED_DOMAINS, False);

        # initialize modeler 
        if( self.remesh_domains_active ):    

            self.InitializeDomains()

            for domain in self.meshing_domains:
                domain.Initialize()
                print( "ComputeInitialAverageMeshParameters")
                domain.ComputeInitialAverageMeshParameters()      
                #domain.Check()
                

    #
    def InitializeDomains(self):

        # initialize the modeler 
        print("::[Remesh_Fluid_Domains]:: Initialize Domains ")
            
        import domain_utilities
        domain_utils = domain_utilities.DomainUtilities()
        
        # find node neighbours
        domain_utils.SearchNodeNeighbours(self.main_model_part, self.echo_level)
            
        # find element neighbours
        domain_utils.SearchElementNeighbours(self.main_model_part, self.echo_level)
            
        # set neighbour search performed
        self.neighbour_search_performed = True

        # set modeler utilities
        self.modeler_utils = KratosPfemBase.ModelerUtilities()
        
        # set the domain labels to conditions
        self.modeler_utils.SetModelPartNameToConditions(self.main_model_part)
        
        # find skin and boundary normals
        if(self.restart == False):
            self.BuildMeshBoundaryForFluids()
            #domain_utils.ConstructModelPartBoundary(self.main_model_part, self.echo_level)

            # search nodal h
            if(self.neighbour_search_performed):
                domain_utils.SearchNodalH(self.main_model_part, self.echo_level)
                                       
        # set the domain labels to nodes
        self.modeler_utils.SetModelPartNameToNodes(self.main_model_part)

        self.main_model_part.ProcessInfo.SetValue(KratosPfemBase.INITIALIZED_DOMAINS, True)

        print(self.main_model_part)
            
    def BuildMeshBoundaryForFluids(self):

        print("::[Remesh_Fluid_Domains_Process]:: Build Mesh Boundary for fluids ")
        # set building options:
        mesh_id = 0

        # define building utility
        #skin_build = KratosPfemFluid.BuildMeshBoundaryForFluids(self.main_model_part, self.echo_level, mesh_id)
        skin_build = KratosPfemBase.BuildMeshBoundary(self.main_model_part, mesh_id, self.echo_level)
 
        # execute building:
        skin_build.Execute()

        print("::[Remesh_Fluid_Domains_Process]:: Mesh Boundary Build executed ")

    ###

    #
    def ExecuteInitializeSolutionStep(self):

        self.step_count += 1
        for domain in self.meshing_domains:
            domain.ComputeAverageMeshParameters()  

        if(self.remesh_domains_active):
            if( self.meshing_before_output ):
                if(self.IsMeshingStep()):
                    print("::[Remesh_Fluid_Domains_Process]:: RemeshFluidDomains ")
                    self.RemeshFluidDomains()
   

    #
    def ExecuteBeforeOutputStep(self):
        
        pass
        #if(self.remesh_domains_active):
             #if( self.meshing_before_output ):
                # if(self.IsMeshingStep()):
                   #  print("::[Remesh_Fluid_Domains_Process]:: RemeshFluidDomains ")
                    # self.RemeshFluidDomains()
        
    #
    def ExecuteAfterOutputStep(self):
        
        if(self.remesh_domains_active):
            if( not self.meshing_before_output ):
                if(self.IsMeshingStep()):
                    self.RemeshFluidDomains()

    ###

    #
    def ExecuteMeshing(domain):
        domain.ExecuteMeshing()

    #
    def RemeshDomains(self):

        if( self.echo_level > 0 ):
            print("::[Meshing_Process]:: MESHING DOMAIN...( call:", self.counter,")")
            
        meshing_options = KratosMultiphysics.Flags()
        #self.model_meshing = KratosPfemBase.ModelMeshing(self.main_model_part, meshing_options, self.echo_level)
        self.model_meshing = KratosPfemFluid.ModelMeshingForFluids(self.main_model_part, meshing_options, self.echo_level)
        
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


        # schedule next meshing
        if(self.meshing_frequency > 0.0): # note: if == 0 always active
            if(self.meshing_control_is_time):
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

    #

    def RemeshFluidDomains(self):

        if(self.remesh_domains_active):
           # if(self.contact_search):
            #    self.ContactTransfer()

            if( self.echo_level > 0 ):
                print("::[Remesh_fluid_domains_process]:: MESH DOMAIN...", self.counter)

            meshing_options = KratosMultiphysics.Flags()
            self.modeler_utils = KratosPfemBase.ModelerUtilities()


            meshing_options.Set(self.modeler_utils.KEEP_ISOLATED_NODES, True)

            #self.model_meshing =  KratosPfemBase.ModelMeshing(self.main_model_part, meshing_options, self.echo_level)
            self.model_meshing =  KratosPfemFluid.ModelMeshingForFluids(self.main_model_part, meshing_options, self.echo_level)

            self.model_meshing.ExecuteInitialize()
         
            print("::[Remesh_fluid_domains_process]:: BEFORE LOOP", self.counter)

            id = 0
            for domain in self.meshing_domains:

                print("::[Remesh_fluid_domains_process]:: IN THE LOOP")

                domain.ExecuteMeshing();

                self.remesh_executed = True

                id+=1

            print("::[Remesh_fluid_domains_process]:: ExecuteFinalize")

            self.model_meshing.ExecuteFinalize()

            self.counter += 1 
