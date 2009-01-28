#importing the Kratos Library
from Kratos import *
from KratosFSIApplication import *

def AddVariables(fluid_model_part,structure_model_part):
    fluid_model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    fluid_model_part.AddNodalSolutionStepVariable(PRESSURE)
    fluid_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    fluid_model_part.AddNodalSolutionStepVariable(AUX)
    fluid_model_part.AddNodalSolutionStepVariable(VAUX)
    fluid_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    fluid_model_part.AddNodalSolutionStepVariable(NORMAL)

    structure_model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    structure_model_part.AddNodalSolutionStepVariable(PRESSURE)
    structure_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    structure_model_part.AddNodalSolutionStepVariable(AUX)
    structure_model_part.AddNodalSolutionStepVariable(VAUX)
    structure_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    structure_model_part.AddNodalSolutionStepVariable(NORMAL)
    
   

    print "variables for the mesh solver added correctly"


class NonConformant_OneSideMap:
    
    def __init__(self,fluid_model_part,structure_model_part):
        self.FluidToStructureMapper = AdvancedNMPointsMapper(fluid_model_part,structure_model_part);
        print "1111111111"
        self.StructureToFluidMapper = AdvancedNMPointsMapper(structure_model_part,fluid_model_part);

        search_radius_factor = 2.0
        (self.FluidToStructureMapper).FindNeighbours(search_radius_factor)
        (self.StructureToFluidMapper).FindNeighbours(search_radius_factor)
        self.it_max = 3
        self.tol = 1e-3

        print (self.FluidToStructureMapper)
        print (self.StructureToFluidMapper)

    def RecomputeTransferPairs(self,search_radius_factor):
        (self.FluidToStructureMapper).FindNeighbours(search_radius_factor)
        (self.StructureToFluidMapper).FindNeighbours(search_radius_factor)
        self.it_max = 3
        self.tol = 1e-3

    def StructureToFluid_VectorMap(self, VectorVar_Origin, VectorVar_Destination ):
        (self.StructureToFluidMapper).VectorMap(VectorVar_Origin, VectorVar_Destination,self.it_max,self.tol)

    def StructureToFluid_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination ):
        (self.StructureToFluidMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination,self.it_max,self.tol)
    
    def FluidToStructure_VectorMap(self, VectorVar_Origin, VectorVar_Destination ):
        (self.FluidToStructureMapper).VectorMap(VectorVar_Origin, VectorVar_Destination,self.it_max,self.tol)

    def FluidToStructure_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination ):
        (self.FluidToStructureMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination,self.it_max,self.tol)
