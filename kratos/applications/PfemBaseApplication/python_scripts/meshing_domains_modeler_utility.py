from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
KratosMultiphysics.CheckForPreviousImport()

#TODO: It must be build as a process
class ModelerUtility:
    #

    def __init__(self, main_model_part, meshing_domains, domain_size):


        self.echo_level = 1        
        self.model_part  = main_model_part
        self.domain_size = domain_size

        # mesh modeler vector
        self.meshing_domains = meshing_domains        
        self.remesh_domains_active = False
        self.neighbours_search_performed = False
        self.counter = 1

        # check if some of them is active:
        if( self.meshing_domains.size() ):
            self.remesh_domains_active = True

                       
    #
    def Initialize(self, restart=False ):
        
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
            if(restart == False):
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
    def RemeshDomains(self):

        if(self.remesh_domains_active):

            if( self.echo_level > 0 ):
                print("::[Modeler_Utility]:: MESH DOMAIN...", self.counter)

            meshing_options = KratosMultiphysics.Flags()
            self.model_meshing =  KratosPfemBase.ModelMeshing(self.model_part, meshing_options, self.echo_level)

            self.model_meshing.ExecuteInitialize()

            for domain in self.meshing_domains:

                domain.ExecuteMeshing();

                self.remesh_executed = True

 
            self.model_meshing.ExecuteFinalize()

            self.counter += 1 


    #
    def GetMeshingStep(self):
        return self.counter

 
