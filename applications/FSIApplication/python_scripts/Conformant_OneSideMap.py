from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
from KratosMultiphysics import *

# Import applications
from KratosMultiphysics.FSIApplication import *

def AddVariables(fluid_model_part, structure_model_part):
    fluid_model_part.AddNodalSolutionStepVariable(PRESSURE)
    fluid_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    fluid_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)

    structure_model_part.AddNodalSolutionStepVariable(PRESSURE)
    structure_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    structure_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)

    print("variables for the mesh solver added correctly")


class Conformant_OneSideMap:

    def __init__(self, fluid_model_part_nodes, structure_model_part_nodes):
        self.FluidToStructureMapper = SharedPointsMapper(fluid_model_part_nodes, structure_model_part_nodes, 1e-9)
        self.StructureToFluidMapper = SharedPointsMapper(structure_model_part_nodes, fluid_model_part_nodes, 1e-9)

    def StructureToFluid_VectorMap(self, VectorVar_Origin, VectorVar_Destination):
        (self.StructureToFluidMapper).VectorMap(VectorVar_Origin, VectorVar_Destination)

    def StructureToFluid_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination):
        (self.StructureToFluidMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination)

    def FluidToStructure_VectorMap(self, VectorVar_Origin, VectorVar_Destination):
        (self.FluidToStructureMapper).VectorMap(VectorVar_Origin, VectorVar_Destination)

    def FluidToStructure_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination):
        (self.FluidToStructureMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination)
