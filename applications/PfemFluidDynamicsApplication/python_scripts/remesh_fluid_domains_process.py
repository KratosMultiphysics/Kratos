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
                if(domain.Active()):
                    domain.ComputeInitialAverageMeshParameters()      
                #domain.Check()
                

    #
    def InitializeDomains(self):

        # initialize the modeler 
        if(self.echo_level>1):
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

        if(self.echo_level>1):
            print(self.main_model_part)
            
    def BuildMeshBoundaryForFluids(self):

              # set building options:
        mesh_id = 0

        # define building utility
        #skin_build = KratosPfemFluid.BuildMeshBoundaryForFluids(self.main_model_part, self.echo_level, mesh_id)
        model_part_name = self.settings["model_part_name"].GetString()
        skin_build = KratosPfemBase.BuildModelPartBoundary(self.main_model_part, model_part_name, self.echo_level)
 
        # execute building:
        skin_build.Execute()

    ###

    #
    def ExecuteInitializeSolutionStep(self):

        self.step_count += 1
        currentTime=self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        currentStep=self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]                

        if(currentStep > 1):
            for domain in self.meshing_domains:
                if(domain.Active()):
                    domain.ComputeAverageMeshParameters()  
                    meanVolumeBeforeMeshing=domain.GetMeanVolume()
                    totalVolumeBeforeMeshing=domain.GetTotalVolume()

                    fileTotalVolume = open("totalVolumeBeforeMeshing.ods", 'a')
                    if(currentStep==2):
                        fileTotalVolume.seek(0)
                        fileTotalVolume.truncate()

                    outstring = str(currentTime) + " " +  str(totalVolumeBeforeMeshing) + " "
                    fileTotalVolume.write(outstring)    
                    fileTotalVolume.close

        volume_acceleration=self.main_model_part.ProcessInfo[KratosMultiphysics.GRAVITY]
        if(currentStep == 1):
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION,volume_acceleration)
  
        if(self.remesh_domains_active):
            if( self.meshing_before_output ):
                if(self.IsMeshingStep()):
                    if(self.echo_level>1):
                        print("::[Remesh_Fluid_Domains_Process]:: RemeshFluidDomains ")
                    self.RemeshFluidDomains()

        if(currentStep > 1):
            for domain in self.meshing_domains:
                if(domain.Active()):
                    domain.ComputeAverageMeshParameters()  
                    meanVolumeAfterMeshing=domain.GetMeanVolume()
                    totalVolumeAfterMeshing=domain.GetTotalVolume()
                    diffMeanVolume=meanVolumeAfterMeshing-meanVolumeBeforeMeshing
                    diffTotalVolume=totalVolumeAfterMeshing-totalVolumeBeforeMeshing
                    fileTotalVolume = open("totalVolumeBeforeMeshing.ods", 'a')
                    
                    outstring =  str(totalVolumeAfterMeshing) + " " +  str(diffTotalVolume) + "\n"
                    fileTotalVolume.write(outstring)    
                    fileTotalVolume.close

      #if(self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] == 1):
          #  for node in self.main_model_part.Nodes:
           #     if (node.Is(KratosMultiphysics.FLUID)):
            #        if(node.IsNot(KratosMultiphysics.RIGID)):
             #           volume_acceleration=node.GetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION)
             #           break
          #  for node in self.main_model_part.Nodes:
            #    if (node.Is(KratosMultiphysics.RIGID)):
             #       node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION,volume_acceleration)
   

    #
    def ExecuteBeforeOutputStep(self):
        
        pass
        #if(self.remesh_domains_active):
             #if( self.meshing_before_output ):
                # if(self.IsMeshingStep()):
                   #  print("::[Remesh_Fluid_Domains_Process]:: RemeshFluidDomains ")
                    # self.RemeshFluidDomains()
        

    def NodalChecksAndAssignations(self):

        numFluid=0
        numRigid=0
        numRigidFluid=0
        numRigidNotFluid=0
        numBoundary=0
        numIsolated=0
        numFreeSurface=0
        numBlocked=0
        mean_nodal_h=0
        for node in self.main_model_part.Nodes:
            # adding dofs
            if (node.Is(KratosMultiphysics.FLUID)):
                numFluid+=1

                nodal_h=node.GetSolutionStepValue(KratosMultiphysics.NODAL_H)
                mean_nodal_h+=nodal_h

                #density=node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
                #print("density ",density)
                
                #viscosity=node.GetSolutionStepValue(KratosMultiphysics.VISCOSITY)
                #print("VISCOSITY ",viscosity)
                
                #bulk_modulus=node.GetSolutionStepValue(KratosMultiphysics.BULK_MODULUS)
                #print("bulk_modulus ",bulk_modulus)

            if (node.Is(KratosMultiphysics.RIGID)):
                numRigid+=1
                if (node.Is(KratosMultiphysics.FLUID)):
                    numRigidFluid+=1
                else:
                    numRigidNotFluid+=1   
                    node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0.0)
                if (node.Is(KratosMultiphysics.BOUNDARY)):
                    numBoundary+=1
                if (node.Is(KratosMultiphysics.ISOLATED)):
                    numIsolated+=1
                    node.SetSolutionStepValue(KratosMultiphysics.PRESSURE,0.0)

                if (node.Is(KratosMultiphysics.FREE_SURFACE)):
                    numFreeSurface+=1
                if (node.Is(KratosMultiphysics.BLOCKED)):
                    numBlocked+=1
 
        mean_nodal_h*=1.0/numFluid;
        if(self.echo_level>1):
            print("nodal_h is  ",nodal_h)
            print("numFluid ",numFluid)
            print("numRigid ",numRigid)
            print("numRigidFluid ",numRigidFluid)
            print("numRigidNotFluid ",numRigidNotFluid)
            print("numBoundary ",numBoundary)
            print("numIsolated ",numIsolated)
            print("numFreeSurface ",numFreeSurface)
            print("numBlocked ",numBlocked)

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

        if( self.echo_level > 1 ):
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

            if( self.echo_level > 1 ):
                print("::[Remesh_fluid_domains_process]:: MESH DOMAIN...", self.counter)

            meshing_options = KratosMultiphysics.Flags()
            self.modeler_utils = KratosPfemBase.ModelerUtilities()


            meshing_options.Set(self.modeler_utils.KEEP_ISOLATED_NODES, True)

            #self.model_meshing =  KratosPfemBase.ModelMeshing(self.main_model_part, meshing_options, self.echo_level)
            self.model_meshing =  KratosPfemFluid.ModelMeshingForFluids(self.main_model_part, meshing_options, self.echo_level)

            self.model_meshing.ExecuteInitialize()
         
            id = 0

            for domain in self.meshing_domains:

                domain.ExecuteMeshing();

                self.remesh_executed = True

                id+=1


            self.model_meshing.ExecuteFinalize()

            self.counter += 1 
