from Kratos import *
from KratosFluidDynamicsApplication import *
from KratosMeshingApplication import *

class DynamicSmagorinsky:

    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.domain_size = domain_size

        # Dynamic Smagorinsky helper class
        self.smagorinsky_util = DynamicSmagorinskyUtils(model_part,domain_size)

        # Mesh refinement tool instance
        if domain_size == 2:
            self.refinement_tool = LocalRefineTriangleMesh(model_part)
        else:
            self.refinement_tool = LocalRefineTetrahedraMesh(model_part)

        # Neigbour search tool instance
        AvgElemNum = 10
        AvgNodeNum = 10
        self.nodal_neighbour_search = FindNodalNeighboursProcess(model_part,\
                                                        AvgElemNum,AvgNodeNum)
        # Find neighbours
        self.nodal_neighbour_search.Execute()

    def Refine(self, num_refinements = 1):

        refine_on_reference = True;
        interpolate_internal_variables = False;

        refinement_step = 0

        while refinement_step < num_refinements:

            #update refinement counter
            refinement_step += 1

            # Mark all elements for refinement
            for element in self.model_part.Elements:
                element.SetValue(SPLIT_ELEMENT,True)

            # Store the mesh before the last refinement step
            # Will be used as coarse mesh in the variational Germano identity
            if refinement_step == num_refinements:
                self.smagorinsky_util.StoreCoarseMesh()

            # Refine
            self.refinement_tool.LocalRefineMesh(refine_on_reference, \
                                                 interpolate_internal_variables)

            # Update neighbours
            self.nodal_neighbour_search.Execute()

    def CorrectBoundary(self,boundary_var = FLAG_VARIABLE):

        # Ensure that the refined mesh has its boundary properly identified
        # (Use on bridge analysis problems only)
        self.smagorinsky_util.CorrectFlagValues(boundary_var)

    def UpdateC(self):
        # Update the Smagorinsky parameter
        self.smagorinsky_util.CalculateC()
