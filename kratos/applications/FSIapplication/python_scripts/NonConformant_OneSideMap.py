from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FSIApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(fluid_model_part, structure_model_part):
    fluid_model_part.AddNodalSolutionStepVariable(NODAL_MAUX)  # Stores Nodal Area
    fluid_model_part.AddNodalSolutionStepVariable(PRESSURE)
    fluid_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    fluid_model_part.AddNodalSolutionStepVariable(MAPPER_SCALAR_PROJECTION_RHS)
    fluid_model_part.AddNodalSolutionStepVariable(MAPPER_VECTOR_PROJECTION_RHS)
    fluid_model_part.AddNodalSolutionStepVariable(VAUX_EQ_TRACTION)
    fluid_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    fluid_model_part.AddNodalSolutionStepVariable(NORMAL)
    fluid_model_part.AddNodalSolutionStepVariable(SCALAR_PROJECTED)
    fluid_model_part.AddNodalSolutionStepVariable(VECTOR_PROJECTED)

    structure_model_part.AddNodalSolutionStepVariable(NODAL_MAUX)  # Stores Nodal Area
    structure_model_part.AddNodalSolutionStepVariable(PRESSURE)
    structure_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    structure_model_part.AddNodalSolutionStepVariable(MAPPER_SCALAR_PROJECTION_RHS)
    structure_model_part.AddNodalSolutionStepVariable(MAPPER_VECTOR_PROJECTION_RHS)
    structure_model_part.AddNodalSolutionStepVariable(VAUX_EQ_TRACTION)
    structure_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    structure_model_part.AddNodalSolutionStepVariable(NORMAL)
    structure_model_part.AddNodalSolutionStepVariable(SCALAR_PROJECTED)
    structure_model_part.AddNodalSolutionStepVariable(VECTOR_PROJECTED)

    print("Mapper variables added correctly.")


class NonConformant_OneSideMap:

    def __init__(self, fluid_model_part, structure_model_part,
                 search_radius_factor=2.0, it_max=25, tol=1e-5):

        # Error check
        if fluid_model_part.NumberOfConditions(0) < 1:
            raise ValueError("No conditions found in the fluid model part, please check that the interface is meshed using SurfaceCondition2D2N/SurfaceCondition3D3N")
        if structure_model_part.NumberOfConditions(0) < 1:
            raise ValueError("No conditions found in the structure model part, please check that the interface is meshed using SurfaceCondition2D2N/SurfaceCondition3D3N")

        search_radius_factor = search_radius_factor
        self.it_max = it_max
        self.tol = tol

        self.Preprocess = InterfacePreprocess()
        self.fl_interface = ModelPart("fluid_interface")
        self.str_interface = ModelPart("structure_interface")
        
        domain_size_fl = fluid_model_part.ProcessInfo[DOMAIN_SIZE]
        domain_size_str = structure_model_part.ProcessInfo[DOMAIN_SIZE]
        
        if (domain_size_fl != domain_size_fl):
          raise ValueError("Domain sizes from two model parts are not compatible")

        print("Identifying fluid interface...")
        if (domain_size_fl == 3):
          self.Preprocess.GenerateTriangleInterfacePart(fluid_model_part, self.fl_interface)
        else:
          self.Preprocess.GenerateLineInterfacePart(fluid_model_part, self.fl_interface)

        print("Identifying structure interface...")
        if (domain_size_fl == 3):
          self.Preprocess.GenerateTriangleInterfacePart(structure_model_part, self.str_interface)
        else:
          self.Preprocess.GenerateLineInterfacePart(structure_model_part, self.str_interface)
        print("Fluid and structure interfaces identified.")
        
        # Interface mappers construction
        self.FluidToStructureMapper = AdvancedNMPointsMapper\
            (self.fl_interface, self.str_interface)
        self.StructureToFluidMapper = AdvancedNMPointsMapper\
            (self.str_interface, self.fl_interface)
        print("Interface Mappers created.")

        # Neighbour search
        (self.FluidToStructureMapper).FindNeighbours(search_radius_factor)
        (self.StructureToFluidMapper).FindNeighbours(search_radius_factor)
        print("Neighbours search finished.")

    def RecomputeTransferPairs(self, search_radius_factor):
        (self.FluidToStructureMapper).FindNeighbours(search_radius_factor)
        (self.StructureToFluidMapper).FindNeighbours(search_radius_factor)

    # Standard mappers
    def StructureToFluid_VectorMap(self, VectorVar_Origin, VectorVar_Destination, sign_pos, distributed):
        (self.StructureToFluidMapper).VectorMap(VectorVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos, distributed)

    def StructureToFluid_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination, sign_pos):
        (self.StructureToFluidMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    def FluidToStructure_VectorMap(self, VectorVar_Origin, VectorVar_Destination, sign_pos, distributed):
        (self.FluidToStructureMapper).VectorMap(VectorVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos, distributed)

    def FluidToStructure_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination, sign_pos):
        (self.FluidToStructureMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)
        
    # Normal vectors
    def StructureToFluid_ScalarToNormalVectorMap(self, ScalarVar_Origin, VectorVar_Destination, sign_pos):
        (self.StructureToFluidMapper).ScalarToNormalVectorMap(ScalarVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def StructureToFluid_NormalVectorToScalarMap(self, VectorVar_Origin, ScalarVar_Destination, sign_pos):
        (self.StructureToFluidMapper).NormalVectorToScalarMap(VectorVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    def FluidToStructure_ScalarToNormalVectorMap(self, ScalarVar_Origin, VectorVar_Destination, sign_pos):
        (self.FluidToStructureMapper).ScalarToNormalVectorMap(ScalarVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def FluidToStructure_NormalVectorToScalarMap(self, VectorVar_Origin, ScalarVar_Destination, sign_pos):
        (self.FluidToStructureMapper).NormalVectorToScalarMap(VectorVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)
