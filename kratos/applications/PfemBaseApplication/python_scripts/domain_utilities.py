from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics 
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
KratosMultiphysics.CheckForPreviousImport()

class DomainUtilities(object):

    #
    def __init__(self):
        pass

    #
    def InitializeDomains(self, model_part, echo_level):

        if( model_part.ProcessInfo[KratosPfemBase.INITIALIZED_DOMAINS] == False ):
            
            # initialize the modeler 
            print("::[Domain_Utilities]:: Initialize", model_part.Name)
                   
            # find node neighbours
            self.SearchNodeNeighbours(model_part, echo_level)
            
            # find element neighbours
            self.SearchElementNeighbours(model_part, echo_level)
            

            # set modeler utilities
            modeler_utils = KratosPfemBase.ModelerUtilities()
        
            # set the domain labels to conditions
            modeler_utils.SetModelPartNameToConditions(model_part)
        
            # find skin and boundary normals
            if( model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False ):
                # construct boundary of a volumetric body domain
                self.ConstructModelPartBoundary(model_part, echo_level)

                # search nodal h
                self.SearchNodalH(model_part, echo_level)
                                       
            # set the domain labels to nodes
            modeler_utils.SetModelPartNameToNodes(model_part)

            model_part.ProcessInfo.SetValue(KratosPfemBase.INITIALIZED_DOMAINS, True)

            print("::[Domain_Utilities]:: Resultant ModelPart")
            print(model_part)
            
        
    #
    def SearchNodeNeighbours(self, model_part, echo_level):

        mesh_id = 0

        # set search options:
        number_of_avg_elems = 10
        number_of_avg_nodes = 10

        # define search utility
        nodal_neighbour_search = KratosPfemBase.NodalNeighboursSearch(model_part, echo_level, number_of_avg_elems, number_of_avg_nodes, mesh_id)

        # execute search:
        nodal_neighbour_search.Execute()

        print("::[Domain_Utilities]:: Nodal Search executed ")

    #
    def SearchElementNeighbours(self, model_part, echo_level):

        mesh_id = 0
        
        domain_size = model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # set search options:
        number_of_avg_elems = 10
         
        # define search utility
        elemental_neighbour_search = KratosPfemBase.ElementalNeighboursSearch(model_part, domain_size, echo_level, number_of_avg_elems, mesh_id)

        # execute search:
        elemental_neighbour_search.Execute()

        if( echo_level > 0 ):
            print("::[Domain_Utilities]:: Elemental Search executed ")


    #
    def ConstructModelPartBoundary(self, model_part, echo_level):

        mesh_id = 0

        print("::[Domain_Utilities]:: Build Mesh Boundary ")
        # set building options:
        

        # define building utility
        skin_build = KratosPfemBase.ConstructModelPartBoundary(model_part, model_part.Name, echo_level)

        # execute building:
        skin_build.Execute()

        if( echo_level > 0 ):
            print("::[Domain_Utilities]:: Mesh Boundary Build executed ")


    ###

    #
    def SearchNodalH(self, model_part, echo_level):

        # define search utility
        nodal_h_search = KratosMultiphysics.FindNodalHProcess(model_part)
        # execute search:
        nodal_h_search.Execute()
        
        # for node in self.main_model_part.Nodes:
        # nodal_h  = node.GetSolutionStepValue(NODAL_H);
        # print "nodal_h:",nodal_h
        
        if( echo_level > 0 ):
            print("::[Domain_Utilities]:: Nodal H Search executed ")

    #
    def ComputeBoundaryNormals(self, model_part, echo_level):

        # define calculation utility
        normals_calculation = KratosPfemBase.BoundaryNormalsCalculation()

        # execute calculation:
        #(scaled normals)
        normals_calculation.CalculateWeightedBoundaryNormals(model_part, echo_level)
        #(unit normals)
        # normals_calculation.CalculateUnitBoundaryNormals(model_part, self.echo_level)

        if( echo_level > 0 ):
            print("::[Domain_Utilities]:: Boundary Normals computed ")


