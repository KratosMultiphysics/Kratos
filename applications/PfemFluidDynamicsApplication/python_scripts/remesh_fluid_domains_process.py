
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
from importlib import import_module

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return RemeshFluidDomainsProcess(Model, settings["Parameters"])


class RemeshFluidDomainsProcess(KratosMultiphysics.Process):
    """The base class for the RemeshFluidDomainsProcess
    """
    def __init__(self, Model, custom_settings ):
        """The constructor of the RemeshFluidDomainsProcess-Object.
        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        parameters -- The ProjectParameters used
        """
        KratosMultiphysics.Process.__init__(self)
        self.main_model_part = Model[custom_settings["model_part_name"].GetString()]

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""{
            "echo_level"                        : 0,
            "model_part_name"                   : "Fluid Domain",
            "meshing_control_type"              : "step",
            "meshing_frequency"                 : 1.0,
            "meshing_before_output"             : true,
            "meshing_domains"                   : [],
            "update_conditions_on_free_surface" : {},
            "write_totalVolumeBeforeMeshing"    : true
        }""")

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.echo_level        = self.settings["echo_level"].GetInt()
        self.dimension         = self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]
        self.meshing_frequency = self.settings["meshing_frequency"].GetDouble()
        self.write_total_volume = self.settings["write_totalVolumeBeforeMeshing"].GetBool()

        self.meshing_control_is_time = False
        meshing_control_type   = self.settings["meshing_control_type"].GetString()
        if meshing_control_type == "time":
            self.meshing_control_is_time = True
        elif meshing_control_type == "step":
            self.meshing_control_is_time = False

        #construct meshing domains
        self.meshing_domains = []
        domains_list = self.settings["meshing_domains"]
        self.number_of_domains = domains_list.size()
        for i in range(0, self.number_of_domains):
            item = domains_list[i]
            python_module_name = "KratosMultiphysics.PfemFluidDynamicsApplication"
            full_module_name = python_module_name + "." + item["python_module"].GetString()
            domain_module = import_module(full_module_name)
            domain = domain_module.CreateMeshingDomain(self.main_model_part,item)
            self.meshing_domains.append(domain)

        # mesh mesher initial values
        self.remesh_domains_active = False
        for domain in self.meshing_domains:
            if domain.Active():
                self.remesh_domains_active = True

        self.neighbours_search_performed = False
        self.step_count   = 1
        self.counter      = 1
        self.next_meshing = 0.0
        self.meshing_before_output = self.settings["meshing_before_output"].GetBool()
        self.update_conditions_on_free_surface = self.settings["update_conditions_on_free_surface"]["update_conditions"].GetBool()

    def ExecuteInitialize(self):
        """This function performs the initialize of the process
        """
        self.fileTotalVolume = None

        if self.meshing_control_is_time:
            self.next_meshing  = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] + self.meshing_frequency
        else:
            self.next_meshing = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + self.meshing_frequency

        self.main_model_part.ProcessInfo.SetValue(KratosDelaunay.INITIALIZED_DOMAINS, False);

        # initialize mesher
        if self.remesh_domains_active:
            self.InitializeDomains()

            for domain in self.meshing_domains:
                domain.SetEchoLevel(self.echo_level)
                domain.Initialize()
                if(domain.Active()):
                    domain.ComputeInitialAverageMeshParameters()
                    domain.SetTimeDataOnProcessInfo()


    def InitializeDomains(self):
        """This function Initializes the Domains
        """
        # initialize the mesher
        if self.echo_level>1:
            print("::[Remesh_Fluid_Domains]:: Initialize Domains ")

        from KratosMultiphysics.DelaunayMeshingApplication import domain_utilities
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
        self.BuildMeshBoundaryForFluids()

            # search nodal h
        if self.neighbour_search_performed:
            self.SearchNodalH(self.main_model_part, self.echo_level)

        # set the domain labels to nodes
        self.mesher_utils.SetModelPartNameToNodes(self.main_model_part)
        self.main_model_part.ProcessInfo.SetValue(KratosDelaunay.INITIALIZED_DOMAINS, True)
        if self.echo_level>1:
            print(self.main_model_part)

    def BuildMeshBoundaryForFluids(self):
        """This function Builds the mesh boundaries for the fluids
        """
        # define building utility
        model_part_name = self.settings["model_part_name"].GetString()

        """
        choose just one of the following two options:
        use this if you want conditions
        ATTENTION: this is slow, and must be used together with ModelMeshingWithConditionsForFluids and GenerateNewConditionsForFluids
        skin_build = KratosDelaunay.BuildModelPartBoundary(self.main_model_part, model_part_name, self.echo_level)

        if you use the following, you will not use/build/compute conditions
        ATTENTION: it must be used together with ModelMeshingForFluids and BuildMeshBoundaryForFluids
        """
        skin_build = KratosPfemFluid.BuildModelPartBoundaryForFluids(self.main_model_part, model_part_name, self.echo_level)
        skin_build.Execute()


    def ExecuteInitializeSolutionStep(self):
        """This function executes the Initialize solution step of the process
        """
        self.step_count += 1
        currentTime=self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        currentStep=self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

        if currentStep >= 2 and self.fileTotalVolume is None and self.write_total_volume:
            self.fileTotalVolume = open("totalVolumeBeforeMeshing.txt",'w')

        if(currentStep > 1 and self.fileTotalVolume is not None):
            for domain in self.meshing_domains:
                if(domain.Active()):
                    domain.ComputeAverageMeshParameters()
                    meanVolumeBeforeMeshing=domain.GetMeanVolume()
                    totalVolumeBeforeMeshing=domain.GetTotalVolume()
                    outstring = str(currentTime) + " " +  str(totalVolumeBeforeMeshing) + " "
                    self.fileTotalVolume.write(outstring)

        volume_acceleration=self.main_model_part.ProcessInfo[KratosMultiphysics.GRAVITY]
        variable_utils = KratosMultiphysics.VariableUtils()
        if(currentStep == 1):
            variable_utils.SetVectorVar(KratosMultiphysics.VOLUME_ACCELERATION, volume_acceleration, self.main_model_part.Nodes)

        if self.remesh_domains_active:
            if self.meshing_before_output:
                if self.IsMeshingStep():
                    if self.echo_level>1:
                        print("::[Remesh_Fluid_Domains_Process]:: RemeshFluidDomains ")
                    self.RemeshFluidDomains()

                    # Updating conditions on the free surface
                    if self.update_conditions_on_free_surface:
                        if self.echo_level > 1:
                            print("::[Remesh_Fluid_Domains_Process]:: UpdateConditionsOnFreeSurface ")
                        KratosPfemFluid.UpdateConditionsOnFreeSurfaceProcess(self.main_model_part, \
                            self.settings["update_conditions_on_free_surface"]).Execute()


        if currentStep > 1 and self.fileTotalVolume is not None:
            for domain in self.meshing_domains:
                if domain.Active():
                    domain.ComputeAverageMeshParameters()
                    meanVolumeAfterMeshing=domain.GetMeanVolume()
                    totalVolumeAfterMeshing=domain.GetTotalVolume()
                    diffMeanVolume=meanVolumeAfterMeshing-meanVolumeBeforeMeshing
                    diffTotalVolume=totalVolumeAfterMeshing-totalVolumeBeforeMeshing
                    outstring =  str(totalVolumeAfterMeshing) + " " +  str(diffTotalVolume) + "\n"
                    self.fileTotalVolume.write(outstring)
        if self.fileTotalVolume is not None:
            self.fileTotalVolume.flush()

    def ExecuteFinalize(self):
        """This function executes the Finalize of the process
        """
        if self.fileTotalVolume is not None:
            self.fileTotalVolume.close()

    #
    def SearchNodalH(self, model_part, echo_level):
        # define search utility
        nodal_h_search = KratosMultiphysics.FindNodalHProcess(model_part)
        # execute search:
        nodal_h_search.Execute()

        # to fix the nodal_h computation of the walls and to avoid setting a small and not representative nodal.H
        # the previous function takes the minimum length of the neighbor fluid elements, this considers only the rigid ones. This avoids fluid leakage through the walls.
        nodal_h_search_for_rigid_walls = KratosPfemFluid.FindNodalHForRigidWallsProcess(model_part)
        # execute search:
        nodal_h_search_for_rigid_walls.Execute()

        if( echo_level > 0 ):
            print("::[Remesh_Fluid_Domains_Process]:: Nodal H Search executed ")

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
                    if self.echo_level>1:
                        print("::[Remesh_Fluid_Domains_Process]:: RemeshFluidDomains ")
                    self.RemeshFluidDomains()

                    # Updating conditions on the free surface
                    if self.update_conditions_on_free_surface:
                        if self.echo_level > 1:
                            print("::[Remesh_Fluid_Domains_Process]:: UpdateConditionsOnFreeSurface ")
                        KratosPfemFluid.FreeSurfaceTrackingProcess(self.main_model_part, self.settings["update_conditions_on_free_surface"]).Execute()


    #
    def ExecuteMeshing(domain):
        domain.ExecuteMeshing()

    #
    def GetMeshingStep(self):
        return self.counter

    #
    def IsMeshingStep(self):

        currentTime=self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        deltaTime=self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        currentStep=self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

        if (self.meshing_control_is_time):
            tolerance=deltaTime*0.000001
            timeDifference=abs(currentTime-self.next_meshing)
            if (tolerance >= timeDifference):
                self.next_meshing=currentTime + self.meshing_frequency
                return True
            else:
                return False
        else:
            tolerance=0.1
            stepDifference=abs(currentStep-self.next_meshing)
            if (tolerance >= stepDifference):
                self.next_meshing= currentStep + self.meshing_frequency
                return True
            else:
                return False

    #

    def RemeshFluidDomains(self):

        if(self.remesh_domains_active):
           # if(self.contact_search):
            #    self.ContactTransfer()

            if( self.echo_level > 1 ):
                print("::[Remesh_fluid_domains_process]:: MESH DOMAIN...", self.counter)

            meshing_options = KratosMultiphysics.Flags()
            self.mesher_utils = KratosDelaunay.MesherUtilities()


            meshing_options.Set(self.mesher_utils.KEEP_ISOLATED_NODES, True)

            ############ choose just one of the following two options: ############
            ## use this if you want conditions
            ## ATTENTION: this is slow, and must be used together with GenerateNewConditionsForFluids and BuildModelPartBoundary
            #self.model_meshing =  KratosPfemFluid.ModelMeshingWithConditionsForFluids(self.main_model_part, meshing_options, self.echo_level)

            ## if you use the following, you will not use/build/compute conditions
            ## ATTENTION: it must be used together with BuildMeshBoundaryForFluids and BuildModelPartBoundaryForFluids
            self.model_meshing =  KratosPfemFluid.ModelMeshingForFluids(self.main_model_part, meshing_options, self.echo_level)
            #######################################################################

            self.model_meshing.ExecuteInitialize()

            id = 0

            for domain in self.meshing_domains:

                domain.ExecuteMeshing();

                self.remesh_executed = True

                id+=1


            self.model_meshing.ExecuteFinalize()

            self.counter += 1
