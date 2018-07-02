from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
KratosMultiphysics.CheckForPreviousImport()

import remesh_domains_process

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return RemeshFluidDomainsProcess(Model, settings["Parameters"])


class RemeshFluidDomainsProcess(remesh_domains_process.RemeshDomainsProcess):
    #
    def __init__(self, Model, custom_settings ):

        super(RemeshFluidDomainsProcess,self).__init__(Model,custom_settings)

    #
    def ExecuteInitialize(self):

        self.fileTotalVolume = None

        remesh_domains_process.RemeshDomainsProcess.ExecuteInitialize(self)

        # compute initial average parameters
        if( self.remesh_domains_active ):

            for domain in self.meshing_domains:
                if(domain.Active()):
                    domain.ComputeInitialAverageMeshParameters()


    #
    def InitializeDomains(self):

        # initialize the mesher
        import domain_utilities
        domain_utils = domain_utilities.DomainUtilities()

        # find node neighbours
        domain_utils.SearchNodeNeighbours(self.main_model_part, self.echo_level)

        # find element neighbours
        domain_utils.SearchElementNeighbours(self.main_model_part, self.echo_level)

        # set neighbour search performed
        self.neighbour_search_performed = True

        # set mesher utilities
        self.mesher_utils = KratosDelaunay.MesherUtilities()

        # set the domain labels to conditions
        self.mesher_utils.SetModelPartNameToConditions(self.main_model_part)

        # find skin and boundary normals
        if(self.restart == False):
            self.BuildMeshBoundaryForFluids()
            #domain_utils.ConstructModelPartBoundary(self.main_model_part, self.echo_level)

            # search nodal h
            if(self.neighbour_search_performed):
                domain_utils.SearchNodalH(self.main_model_part, self.echo_level)

        # set the domain labels to nodes
        self.mesher_utils.SetModelPartNameToNodes(self.main_model_part)

        self.main_model_part.ProcessInfo.SetValue(KratosDelaunay.INITIALIZED_DOMAINS, True)

        if(self.echo_level>1):
            print(self.main_model_part)


    def BuildMeshBoundaryForFluids(self):

        # set building options:

        # define building utility
        model_part_name = self.settings["model_part_name"].GetString()
        skin_build = KratosDelaunay.BuildModelPartBoundary(self.main_model_part, model_part_name, self.echo_level)

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
            #self.probe1 = open("probe1.txt",'w')
            #self.probe2 = open("probe2.txt",'w')
            #self.probe3 = open("probe3.txt",'w')

        if(currentStep > 1 and self.fileTotalVolume is not None):
            #maxYprobe1=0.1
            #maxYprobe2=0.1
            #maxYprobe3=0.1
            #for node in self.main_model_part.Nodes:
                #if(node.IsNot(KratosMultiphysics.ISOLATED)):
                    #if(node.X>5.9 and node.X<6.1):
                        #if(node.Y>maxYprobe1):
                            #maxYprobe1=node.Y
                    #if(node.X>8.9 and node.X<9.1):
                        #if(node.Y>maxYprobe2):
                            #maxYprobe2=node.Y
                    #if(node.X>11.9 and node.X<12.1):
                        #if(node.Y>maxYprobe3):
                            #maxYprobe3=node.Y

            #outstring = str(currentTime) + " " +  str(maxYprobe1) + "\n"
            #self.probe1.write(outstring)
            #outstring = str(currentTime) + " " +  str(maxYprobe2) + "\n"
            #self.probe2.write(outstring)
            #outstring = str(currentTime) + " " +  str(maxYprobe3) + "\n"
            #self.probe3.write(outstring)

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

        if(self.remesh_domains_active):
            if( self.meshing_before_output ):
                if(self.IsMeshingStep()):
                    if(self.echo_level>1):
                        print(self._class_prefix()+" RemeshFluidDomains ")
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
            #self.probe1.flush()
            #self.probe2.flush()
            #self.probe3.flush()


    def ExecuteFinalize(self):
        if self.fileTotalVolume is not None:
            self.fileTotalVolume.close()
            #self.probe1.close()
            #self.probe2.close()
            #self.probe3.close()


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
                   #  print(self._class_prefix()+" RemeshFluidDomains ")
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
                if(self.IsMeshingStep()):
                    self.RemeshFluidDomains()

    #
    def ExecuteMeshing(domain):
        domain.ExecuteMeshing()

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
                print(self._class_prefix()+" MESH DOMAIN...", self.counter)

            meshing_options = KratosMultiphysics.Flags()
            self.mesher_utils = KratosDelaunay.MesherUtilities()


            meshing_options.Set(self.mesher_utils.KEEP_ISOLATED_NODES, True)

            #self.model_structure =  KratosDelaunay.ModelStructure(self.main_model_part, meshing_options, self.echo_level)
            self.model_structure =  KratosPfem.FluidModelStructure(self.main_model_part, meshing_options, self.echo_level)

            self.model_structure.ExecuteInitialize()

            id = 0

            for domain in self.meshing_domains:

                domain.ExecuteMeshing();

                self.remesh_executed = True

                id+=1


            self.model_structure.ExecuteFinalize()

            self.counter += 1

    #
    def _class_prefix(self):
        header = "::[---Meshing_Fluid---]::"
        return header
