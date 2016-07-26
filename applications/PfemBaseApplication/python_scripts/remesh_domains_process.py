from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return RemeshDomainsProcess(Model, settings["Parameters"])


class RemeshDomainsProcess(KratosMultiphysics.Process):
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
        self.domain_size       = self.model_part.ProcessInfo[DOMAIN_SIZE]
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
        if( self.number_of_domains ):
            self.remesh_domains_active = True

        self.neighbours_search_performed = False
        self.step_count   = 1
        self.counter      = 1
        self.next_meshing = 0.0
        self.meshing_before_output = self.settings["meshing_before_output"].GetBool()
                       
    #
    def ExecuteInitialize(self):

        self.restart = False
        if( self.model_part.ProcessInfo[IS_RESTARTED] ):
            self.restart = True
        
        # initialize the modeler 
        if( self.remesh_domains_active ):        
            print("::[Modeler_Utility]:: Initialize Domains ")
            
            # find node neighbours
            self.SearchNodeNeighbours()
            
            # find element neighbours
            self.SearchElementNeighbours()
            
            # set neighbour search performed
            self.neighbour_search_performed = True

            # find skin and boundary normals
            if(self.restart == False):
                self.BuildMeshBoundary()

                # search nodal h
                # self.SearchNodalH() #now done from main script
            
                
            # set modeler utilities
            self.modeler_utils = KratosPfemBase.ModelerUtilities()

            # set the domain labels to mesh modeler
            self.modeler_utils.SetDomainLabels(self.model_part)

            for domain in self.meshing_domains:
                domain.Check();

    #
    def SearchNodeNeighbours(self):

        mesh_id = 0

        # set search options:
        number_of_avg_elems = 10
        number_of_avg_nodes = 10

        # define search utility
        nodal_neighbour_search = KratosPfemBase.NodalNeighboursSearch(self.model_part, self.echo_level, number_of_avg_elems, number_of_avg_nodes, mesh_id)

        # execute search:
        nodal_neighbour_search.Execute()

        print("::[Modeler_Utility]:: Nodal Search executed ")

    #
    def SearchElementNeighbours(self):

        mesh_id = 0

        # set search options:
        number_of_avg_elems = 10
         
        # define search utility
        elemental_neighbour_search = KratosPfemBase.ElementalNeighboursSearch(self.model_part, self.domain_size, self.echo_level, number_of_avg_elems, mesh_id)

        # execute search:
        elemental_neighbour_search.Execute()

        if( self.echo_level > 0 ):
            print("::[Modeler_Utility]:: Elemental Search executed ")


    #
    def BuildMeshBoundary(self):

        mesh_id = 0

        print("::[Modeler_Utility]:: Build Mesh Boundary ")
        # set building options:
        

        # define building utility
        skin_build = KratosPfemBase.BuildMeshBoundary(self.model_part, mesh_id, self.echo_level)

        # execute building:
        skin_build.Execute()

        if( self.echo_level > 0 ):
            print("::[Modeler_Utility]:: Mesh Boundary Build executed ")


    ###

    #
    def SearchNodalH(self):

        if(self.neighbour_search_performed):
            # define search utility
            nodal_h_search = KratosMultiphysics.FindNodalHProcess(self.model_part)
            # execute search:
            nodal_h_search.Execute()

            # for node in self.model_part.Nodes:
                # nodal_h  = node.GetSolutionStepValue(NODAL_H);
                # print "nodal_h:",nodal_h

            if( self.echo_level > 0 ):
                print("::[Modeler_Utility]:: Nodal H Search executed ")

    #
    def ComputeBoundaryNormals(self):

        # define calculation utility
        normals_calculation = KratosPfemBase.BoundaryNormalsCalculation()

        # execute calculation:
        #(scaled normals)
        normals_calculation.CalculateBoundaryNormals(self.model_part, self.echo_level)
        #(unit normals)
        # normals_calculation.CalculateBoundaryUnitNormals(model_part, self.echo_level)

        if( self.echo_level > 0 ):
            print("::[Modeler_Utility]:: Boundary Normals computed ")


    ###

    #
    def ExecuteInitializeSolutionStep(self):

        self.step_count += 1

    #
    def ExecuteBeforeOutputStep(self):
        
        if(self.remesh_domains_active):
            if( self.meshing_before_output ):
                if(self.IsMeshingStep):
                    self.RemeshDomains()
        
    #
    def ExecuteAfterOutputStep(self):
        
        if(self.remesh_domains_active):
            if( !self.meshing_before_output ):
                if(self.IsMeshingStep):
                    self.RemeshDomains()

    ###

    #
    def RemeshDomains(self):

        if( self.echo_level > 0 ):
            print("::[Meshing_Process]:: MESH DOMAIN...", self.counter)
            
        meshing_options = KratosMultiphysics.Flags()
        self.model_meshing = KratosPfemBase.ModelMeshing(self.model_part, meshing_options, self.echo_level)
        
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
                        

    #
    def GetMeshingStep(self):
        return self.counter

    #
    def IsMeshingStep(self):

        if(self.meshing_control_is_time):
            #print( str(self.model_part.ProcessInfo[TIME])+">"+ str(self.next_meshing) )
            return ( self.model_part.ProcessInfo[TIME] > self.next_meshing )
        else:
            return ( self.step_count >= self.next_meshing )
