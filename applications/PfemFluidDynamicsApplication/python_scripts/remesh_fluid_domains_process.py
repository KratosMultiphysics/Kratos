from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.PfemApplication as KratosPfem
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
KratosMultiphysics.CheckForPreviousImport()

import remesh_domains_process

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return RemeshFluidDomainsProcess(Model, settings["Parameters"])


#class RemeshFluidDomainsProcess(KratosMultiphysics.Process):
class RemeshFluidDomainsProcess(remesh_domains_process.RemeshDomainsProcess):
    #

     
    #
    def ExecuteInitialize(self):

        self.fileTotalVolume = None
        self.probe1isolated = None
        self.probe1 = None
        self.probe2 = None
        self.probe3 = None
        self.probe4 = None
        self.probe5 = None
        self.probe6 = None
        self.probe7 = None
        self.probe8 = None
        self.probe9 = None

        print("::[FLUID Meshing_Process]:: meshing frequency", self.meshing_frequency)

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
            self.step_count = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
            if self.meshing_control_is_time:
                self.next_meshing  = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] + self.meshing_frequency
            else:
                self.next_meshing = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + self.meshing_frequency
            #self.meshing_output = self.meshing_frequency


        self.main_model_part.ProcessInfo.SetValue(KratosPfem.INITIALIZED_DOMAINS, False);

        # initialize modeler 
        if( self.remesh_domains_active ):    

            self.InitializeDomains()

            for domain in self.meshing_domains:
                domain.SetEchoLevel(self.echo_level)
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
        self.modeler_utils = KratosPfem.ModelerUtilities()
        
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

        self.main_model_part.ProcessInfo.SetValue(KratosPfem.INITIALIZED_DOMAINS, True)

        if(self.echo_level>1):
            print(self.main_model_part)
            
    def BuildMeshBoundaryForFluids(self):

        # set building options:
        
        # define building utility
        model_part_name = self.settings["model_part_name"].GetString()
        skin_build = KratosPfem.BuildModelPartBoundary(self.main_model_part, model_part_name, self.echo_level)
 
        # execute building:
        skin_build.Execute()

    ###

    #
    def ExecuteInitializeSolutionStep(self):

        self.step_count += 1
        currentTime=self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        currentStep=self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]                

        if currentStep >= 2 and self.fileTotalVolume is None:
            self.fileTotalVolume = open("totalVolumeBeforeMeshing.txt",'w')
            self.probe1isolated = open("probe1isolated.txt",'w')
            self.probe1 = open("probe1.txt",'w')
            self.probe2 = open("probe2.txt",'w')
            self.probe3 = open("probe3.txt",'w')
            self.probe4 = open("probe4.txt",'w')
            self.probe5 = open("probe5.txt",'w')
            self.probe6 = open("probe6.txt",'w')
            self.probe7 = open("probe7.txt",'w')
            self.probe8 = open("probe8.txt",'w')
            self.probe9 = open("probe9.txt",'w')

        if(currentStep > 1 and self.fileTotalVolume is not None):
            maxYprobe1isolated=0.1
            maxYprobe1=0.1
            maxYprobe2=0.1
            maxYprobe3=0.1
            maxYprobe4=0.1
            maxYprobe5=0.1
            maxYprobe6=0.1
            maxYprobe7=0.1
            maxYprobe8=0.1
            maxYprobe9=0.1
            for node in self.main_model_part.Nodes:
                if(node.X>1.87 and node.X<1.93):
                    if(node.Y>maxYprobe1isolated):
                        maxYprobe1isolated=node.Y
                if(node.IsNot(KratosMultiphysics.ISOLATED)):
                    if(node.X>1.87 and node.X<1.93):
                        if(node.Y>maxYprobe1):
                            maxYprobe1=node.Y
                    if(node.X>3.07 and node.X<3.13):
                        if(node.Y>maxYprobe2):
                            maxYprobe2=node.Y
                    if(node.X>5.97 and node.X<6.03):
                        if(node.Y>maxYprobe3):
                            maxYprobe3=node.Y
                    if(node.X>8.77 and node.X<8.83):
                        if(node.Y>maxYprobe4):
                            maxYprobe4=node.Y
                    if(node.X>11.27 and node.X<11.33):
                        if(node.Y>maxYprobe5):
                            maxYprobe5=node.Y
                    if(node.X>16.77 and node.X<16.83):
                        if(node.Y>maxYprobe6):
                            maxYprobe6=node.Y
                    if(node.X>21.97 and node.X<22.03):
                        if(node.Y>maxYprobe7):
                            maxYprobe7=node.Y
                    if(node.X>26.57 and node.X<26.63):
                        if(node.Y>maxYprobe8):
                            maxYprobe8=node.Y
                    if(node.X>32.47 and node.X<32.53):
                        if(node.Y>maxYprobe9):
                            maxYprobe9=node.Y

            outstring = str(currentTime) + " " +  str(maxYprobe1isolated) + "\n"
            self.probe1isolated.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe1) + "\n"
            self.probe1.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe2) + "\n"
            self.probe2.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe3) + "\n"
            self.probe3.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe4) + "\n"
            self.probe4.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe5) + "\n"
            self.probe5.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe6) + "\n"
            self.probe6.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe7) + "\n"
            self.probe7.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe8) + "\n"
            self.probe8.write(outstring)
            outstring = str(currentTime) + " " +  str(maxYprobe9) + "\n"
            self.probe9.write(outstring)
            
            for domain in self.meshing_domains:
                if(domain.Active()):
                    domain.ComputeAverageMeshParameters()  
                    meanVolumeBeforeMeshing=domain.GetMeanVolume()
                    totalVolumeBeforeMeshing=domain.GetTotalVolume()
                    outstring = str(currentTime) + " " +  str(totalVolumeBeforeMeshing) + " "
                    self.fileTotalVolume.write(outstring)
                    #fileTotalVolume = open("totalVolumeBeforeMeshing.txt", 'a')
                    #if(currentStep==2):
                        #fileTotalVolume.seek(0)
                        #fileTotalVolume.truncate()

                    #fileTotalVolume.write(outstring)    
                    #fileTotalVolume.close

        volume_acceleration=self.main_model_part.ProcessInfo[KratosMultiphysics.GRAVITY]
        if(currentStep == 1):
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION,volume_acceleration)
  
        if(self.remesh_domains_active or currentStep == 1):
            if( self.meshing_before_output or currentStep == 1):
                if(self.IsMeshingStep() or currentStep == 1):
                    if(self.echo_level>1 or self.meshing_frequency>2 or self.meshing_control_is_time):
                        print("--> Remesh Fluid Domain at", currentTime, 's')
                    self.RemeshFluidDomains()

        if(currentStep > 1 and self.fileTotalVolume is not None):
            for domain in self.meshing_domains:
                if(domain.Active()):
                    domain.ComputeAverageMeshParameters()  
                    meanVolumeAfterMeshing=domain.GetMeanVolume()
                    totalVolumeAfterMeshing=domain.GetTotalVolume()
                    diffMeanVolume=meanVolumeAfterMeshing-meanVolumeBeforeMeshing
                    diffTotalVolume=totalVolumeAfterMeshing-totalVolumeBeforeMeshing
                    #fileTotalVolume = open("totalVolumeBeforeMeshing.txt", 'a')
                    
                    outstring =  str(totalVolumeAfterMeshing) + " " +  str(diffTotalVolume) + "\n"
                    #fileTotalVolume.write(outstring)    
                    #fileTotalVolume.close
                    self.fileTotalVolume.write(outstring)
        if self.fileTotalVolume is not None:
            self.fileTotalVolume.flush()
            self.probe1isolated.flush()
            self.probe1.flush()
            self.probe2.flush()
            self.probe3.flush()
            self.probe4.flush()
            self.probe5.flush()
            self.probe6.flush()
            self.probe7.flush()
            self.probe8.flush()
            self.probe9.flush()


    def ExecuteFinalize(self):
        if self.fileTotalVolume is not None:
            self.fileTotalVolume.close()
            self.probe1isolated.close()
            self.probe1.close()
            self.probe2.close()
            self.probe3.close()
            self.probe4.close()
            self.probe5.close()
            self.probe6.close()
            self.probe7.close()
            self.probe8.close()
            self.probe9.close()
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
                self.main_model_part.ProcessInfo[KratosPfem.MESHING_STEP_PERFORMED] = False
                if(self.IsMeshingStep()):
                    self.RemeshFluidDomains()

   # 
    def RemeshFluidDomains(self):

        if(self.remesh_domains_active or self.counter==1):
           # if(self.contact_search):
            #    self.ContactTransfer()

            if( self.echo_level > 1 ):
                print("::[Remesh_fluid_domains_process]:: MESH DOMAIN...", self.counter)

            meshing_options = KratosMultiphysics.Flags()
            self.modeler_utils = KratosPfem.ModelerUtilities()
            meshing_options.Set(self.modeler_utils.KEEP_ISOLATED_NODES, True)
            #meshing_options.Set(self.modeler_utils.KEEP_ISOLATED_NODES, False)

            #self.model_meshing =  KratosPfem.ModelMeshing(self.main_model_part, meshing_options, self.echo_level)
            self.model_meshing =  KratosPfemFluid.ModelMeshingForFluids(self.main_model_part, meshing_options, self.echo_level)

            self.model_meshing.ExecuteInitialize()
         

            for domain in self.meshing_domains:
                domain.ExecuteMeshing();
                self.remesh_executed = True

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
