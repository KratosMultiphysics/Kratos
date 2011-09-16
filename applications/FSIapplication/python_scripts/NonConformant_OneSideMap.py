#importing the Kratos Library
from Kratos import *
from KratosFSIApplication import *

def AddVariables(fluid_model_part,structure_model_part):
    fluid_model_part.AddNodalSolutionStepVariable(NODAL_MAUX)# Stores Nodal Area
    fluid_model_part.AddNodalSolutionStepVariable(PRESSURE)
    fluid_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    fluid_model_part.AddNodalSolutionStepVariable(AUX)
    fluid_model_part.AddNodalSolutionStepVariable(VAUX)
    fluid_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    fluid_model_part.AddNodalSolutionStepVariable(NORMAL)

    structure_model_part.AddNodalSolutionStepVariable(NODAL_MAUX)
    structure_model_part.AddNodalSolutionStepVariable(PRESSURE)
    structure_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    structure_model_part.AddNodalSolutionStepVariable(AUX)
    structure_model_part.AddNodalSolutionStepVariable(VAUX)
    structure_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    structure_model_part.AddNodalSolutionStepVariable(NORMAL)
    
   

    print "variables for the mesh solver added correctly"


class NonConformant_OneSideMap:
    
    def __init__(self,fluid_model_part,structure_model_part,\
                 search_radius_factor = 2.0, it_max = 3,tol = 1e-3):

        # Error check
	if fluid_model_part.NumberOfConditions(0) < 1:
            raise ValueError("No conditions found in the fluid model part, please check that the interface is meshed using Condition3D")
        if structure_model_part.NumberOfConditions(0) < 1:
            raise ValueError("No conditions found in the structure model part, please check that the interface is meshed using Face3D3N")

        search_radius_factor = search_radius_factor
        self.it_max = it_max
        self.tol = tol

        self.Preprocess = InterfacePreprocess()
        self.fl_interface = ModelPart("fluid_interface")
        self.str_interface = ModelPart("structure_interface")

        print "Identifying fluid interface"
        self.Preprocess.GenerateInterfacePart(fluid_model_part,self.fl_interface)
        print "Identifying structure interface"
        self.Preprocess.GenerateInterfacePart(structure_model_part,self.str_interface)
        print "Interface identified"
        
        self.FluidToStructureMapper = AdvancedNMPointsMapper\
                                      (self.fl_interface,self.str_interface);
        self.StructureToFluidMapper = AdvancedNMPointsMapper\
                                      (self.str_interface,self.fl_interface);
        print "Interface Mappers created"

        (self.FluidToStructureMapper).FindNeighbours(search_radius_factor)
        (self.StructureToFluidMapper).FindNeighbours(search_radius_factor)

        print (self.FluidToStructureMapper)
        print (self.StructureToFluidMapper)

    def RecomputeTransferPairs(self,search_radius_factor):
        (self.FluidToStructureMapper).FindNeighbours(search_radius_factor)
        (self.StructureToFluidMapper).FindNeighbours(search_radius_factor)

    def StructureToFluid_VectorMap(self, VectorVar_Origin, VectorVar_Destination ):
        (self.StructureToFluidMapper).VectorMap(VectorVar_Origin, VectorVar_Destination,self.it_max,self.tol)

    def StructureToFluid_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination ):
        (self.StructureToFluidMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination,self.it_max,self.tol)
    
    def FluidToStructure_VectorMap(self, VectorVar_Origin, VectorVar_Destination ):
        (self.FluidToStructureMapper).VectorMap(VectorVar_Origin, VectorVar_Destination,self.it_max,self.tol)

    def FluidToStructure_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination ):
        (self.FluidToStructureMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination,self.it_max,self.tol)
