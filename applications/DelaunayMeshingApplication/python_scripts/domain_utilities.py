from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
KratosMultiphysics.CheckForPreviousImport()

class DomainUtilities(object):

    #
    def __init__(self):
        pass

    #
    def InitializeDomains(self, model_part, echo_level):

        if( model_part.ProcessInfo[KratosDelaunay.INITIALIZED_DOMAINS] == False ):

            # initialize the mesher
            print("::[--Domain Utilities-]:: Initialize", model_part.Name)

            # find node neighbours
            self.SearchNodeNeighbours(model_part, echo_level)

            # find element neighbours
            self.SearchElementNeighbours(model_part, echo_level)

            # set mesher utilities
            mesher_utils = KratosDelaunay.MesherUtilities()

            # set the domain labels to conditions
            mesher_utils.SetModelPartNameToConditions(model_part)

            # find skin and boundary normals
            if( model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == False ):
                # build boundary of a volumetric body domain
                self.BuildModelPartBoundary(model_part, echo_level)

                # search nodal h
                self.SearchNodalH(model_part, echo_level)

                # add rigid and solid boundary nodes to fluid domains:
                self.AddBoundaryNodesToFluidDomains(model_part)

                # set the domain labels to nodes
                mesher_utils.SetModelPartNameToNodes(model_part)


            model_part.ProcessInfo.SetValue(KratosDelaunay.INITIALIZED_DOMAINS, True)

            if( echo_level > 0 ):
                print("::[--Domain Utilities-]:: Resultant ModelPart")
                print(model_part)


    #
    @classmethod
    def SearchNodeNeighbours(self, model_part, echo_level):


        # set search options:
        number_of_avg_elems = 10
        number_of_avg_nodes = 10

        # define search utility
        nodal_neighbour_search = KratosDelaunay.NodalNeighboursSearch(model_part, echo_level, number_of_avg_elems, number_of_avg_nodes)

        # execute search:
        nodal_neighbour_search.Execute()

        print("::[--Domain Utilities-]:: Nodal Search executed ")

    #
    @classmethod
    def SearchElementNeighbours(self, model_part, echo_level):

        dimension = model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]
        # set search options:
        number_of_avg_elems = 10

        # define search utility
        elemental_neighbour_search = KratosDelaunay.ElementalNeighboursSearch(model_part, dimension, echo_level, number_of_avg_elems)

        # execute search:
        elemental_neighbour_search.Execute()

        if( echo_level > 0 ):
            print("::[--Domain Utilities-]:: Elemental Search executed ")


    #
    @classmethod
    def BuildModelPartBoundary(self, model_part, echo_level):


        print("::[--Domain Utilities-]:: Build Mesh Boundary ")
        # set building options:


        # define building utility
        skin_build = KratosDelaunay.BuildModelPartBoundary(model_part, model_part.Name, echo_level)

        # execute building:
        skin_build.Execute()

        # search condition masters: (check)
        # skin_build.SearchConditionMasters()

        if( echo_level > 0 ):
            print("::[--Domain Utilities-]:: Mesh Boundary Build executed ")


    ###

    #
    @classmethod
    def SearchNodalH(self, model_part, echo_level):

        # define search utility
        nodal_h_search = KratosMultiphysics.FindNodalHProcess(model_part)
        # execute search:
        nodal_h_search.Execute()

        # for node in self.main_model_part.Nodes:
        # nodal_h  = node.GetSolutionStepValue(NODAL_H);
        # print "nodal_h:",nodal_h

        if( echo_level > 0 ):
            print("::[--Domain Utilities-]:: Nodal H Search executed ")

    #
    @classmethod
    def ComputeBoundaryNormals(self, model_part, echo_level):

        # define calculation utility
        normals_calculation = KratosDelaunay.BoundaryNormalsCalculation()

        # execute calculation:
        #(scaled normals)
        normals_calculation.CalculateWeightedBoundaryNormals(model_part, echo_level)
        #(unit normals)
        # normals_calculation.CalculateUnitBoundaryNormals(model_part, self.echo_level)

        if( echo_level > 0 ):
            print("::[--Domain Utilities-]:: Boundary Normals computed ")


    #
    @classmethod
    def AddBoundaryNodesToFluidDomains(self,model_part):

        exist_fluid_domain = False
        for part in model_part.SubModelParts:
            if part.Is(KratosMultiphysics.FLUID):
               exist_fluid_domain = True
               break

        if( exist_fluid_domain ):

            print("::[--Domain Utilities-]:: Add boundary nodes to fluid domains ")

            try:
                import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
            except:
                raise Exception("SolidMechanicsApplication not imported and needed in this operation")

            transfer_flags = [KratosMultiphysics.BOUNDARY,KratosMultiphysics.NOT_FLUID]
            entity_type = "Nodes"
            for fluid_part in model_part.SubModelParts:
                if (fluid_part.IsNot(KratosMultiphysics.ACTIVE) and fluid_part.Is(KratosMultiphysics.FLUID)):
                    for part in model_part.SubModelParts:
                        if part.IsNot(KratosMultiphysics.ACTIVE):
                            if( part.Is(KratosMultiphysics.SOLID) or part.Is(KratosMultiphysics.RIGID) ):
                                transfer_process = KratosSolid.TransferEntitiesProcess(fluid_part,part,entity_type,transfer_flags)
                                transfer_process.Execute()
