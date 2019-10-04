from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay


class ModelerUtility:
    #

    def __init__(self, model_part, dimension, remesh_domains):

        self.echo_level = 1        
        self.model_part = model_part
        self.dimension = dimension

        # set remesh flags
        self.modeler_active = False
        self.remesh_domains = remesh_domains
        self.neighbours_set = False

        # mesh modeler vector
        self.counter = 1
        self.mesh_modelers = []

        # mesh modeler parameters
        self.alpha_shape        = 1.4
        self.h_factor           = 0.5
        self.avoid_tip_elements = False
        self.offset_factor      = 0
        self.remesh_frequencies = []

        # time step meshing control parameters
        self.remesh_executed = False

    #
    def Initialize(self):

        self.remesh_executed = False
        
    #
    def InitializeDomains(self, ReloadFile = False):

        if( self.modeler_active ):        
            print("::[Modeler_Utility]:: Initialize Domains ")
            
            # set active search
            self.search_active = False

            if(self.remesh_domains):
                self.search_active = True

            self.neighbours_set = False
            if(self.search_active):
                # find neighbours
                self.SearchNeighbours()
                # find skin and boundary normals
                if(ReloadFile == False):
                    self.BuildMeshBoundaryForFluids()
                    #self.BuildMeshBoundary()
                self.neighbours_set = True

    #
    def SearchNeighbours(self):

        print("::[Modeler_Utility]:: Search Neighbours ")

        self.SearchNodeNeighbours()
        self.SearchElementNeighbours()

    #
    def SearchNodeNeighbours(self):

        # set search options:
        number_of_avg_elems = 10
        number_of_avg_nodes = 10

        # define search utility
        nodal_neighbour_search = KratosDelaunay.NodalNeighboursSearch(self.model_part, self.echo_level, number_of_avg_elems, number_of_avg_nodes)

        # execute search:
        nodal_neighbour_search.Execute()

        print("::[Modeler_Utility]:: Nodal Search executed ")

    #
    def SearchElementNeighbours(self):

        # set search options:
        number_of_avg_elems = 10
         
        # define search utility
        elemental_neighbour_search = KratosDelaunay.ElementalNeighboursSearch(self.model_part, self.dimension, self.echo_level, number_of_avg_elems)

        # execute search:
        elemental_neighbour_search.Execute()

        print("::[Modeler_Utility]:: Elemental Search executed ")

    #
    def ComputeBoundaryNormals(self):

        # define calculation utility
        # normals_calculation = BoundaryNormalsCalculation()
        normals_calculation = KratosDelaunay.BoundaryNormalsCalculation()

        # execute calculation:
        #(scaled normals)
        normals_calculation.CalculateWeightedBoundaryNormals(self.model_part, self.echo_level)
        #(unit normals)
        # normals_calculation.CalculateUnitBoundaryNormals(model_part, self.echo_level)

        print("::[Modeler_Utility]:: Boundary Normals computed ")

    #
    def BuildMeshBoundary(self):

        print("::[Modeler_Utility]:: Build Mesh Boundary ")
        # set building options:

        # define building utility
        # skin_build = BuildMeshBoundary(self.model_part, self.dimension, self.echo_level)
        skin_build = KratosDelaunay.BuildMeshBoundary(self.model_part, self.echo_level)

        # execute building:
        skin_build.Execute()

        print("::[Modeler_Utility]:: Mesh Boundary Build executed ")

    #
    def SearchNodalH(self):

        if(self.neighbours_set):
            # define search utility
            nodal_h_search = KratosMultiphysics.FindNodalHProcess(self.model_part)
            # execute search:
            nodal_h_search.Execute()

            # for node in self.model_part.Nodes:
                # nodal_h  = node.GetSolutionStepValue(NODAL_H);
                # print "nodal_h:",nodal_h

            print("::[Modeler_Utility]:: Nodal H Search executed ")

        #
    def ComputeAverageMeshParameters(self):
     
        for domain in self.meshing_domains:
            if(domain.Active()):
                domain.ComputeAverageMeshParameters()       
#
    def ComputeInitialAverageMeshParameters(self):     

        for domain in self.meshing_domains:
            if(domain.Active()):
                domain.ComputeInitialAverageMeshParameters()  
                domain.SetTimeDataOnProcessInfo()     
#

    def GetRemeshFrequency(self):
        
        remesh_frequency = 0
        for size in range(0,len(self.remesh_frequencies)):
            if((remesh_frequency > self.remesh_frequencies[size]) or remesh_frequency == 0):
                remesh_frequency = self.remesh_frequencies[size]
        
        return remesh_frequency
            

    #
    def InitializeStep(self):

        self.remesh_executed = False

    #
