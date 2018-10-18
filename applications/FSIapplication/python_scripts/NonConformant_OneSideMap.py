from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7


# Importing the Kratos Library
from KratosMultiphysics import *

# Check that applications were imported in the main scriptº
CheckRegisteredApplications("FSIApplication")

# Import applications
from KratosMultiphysics.FSIApplication import *

def AddVariables(fluid_model_part, structure_model_part):
    fluid_model_part.AddNodalSolutionStepVariable(PRESSURE)
    fluid_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    fluid_model_part.AddNodalSolutionStepVariable(MAPPER_SCALAR_PROJECTION_RHS)
    fluid_model_part.AddNodalSolutionStepVariable(MAPPER_VECTOR_PROJECTION_RHS)
    fluid_model_part.AddNodalSolutionStepVariable(VAUX_EQ_TRACTION)
    fluid_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    fluid_model_part.AddNodalSolutionStepVariable(NORMAL)
    fluid_model_part.AddNodalSolutionStepVariable(SCALAR_PROJECTED)
    fluid_model_part.AddNodalSolutionStepVariable(VECTOR_PROJECTED)

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
        self.fl_interface = fluid_model_part.GetOwnerModel().CreateModelPart("fluid_interface")
        self.str_interface = structure_model_part.GetOwnerModel().CreateModelPart("structure_interface")

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
        if distributed:
            VariableRedistributionUtility.DistributePointValues(self.str_interface,VectorVar_Origin,VAUX_EQ_TRACTION,self.tol,self.it_max)
            self.StructureToFluidMapper.VectorMap(VAUX_EQ_TRACTION, VAUX_EQ_TRACTION, self.it_max, self.tol, sign_pos)
            VariableRedistributionUtility.ConvertDistributedValuesToPoint(self.fl_interface,VAUX_EQ_TRACTION,VectorVar_Destination)
        else:
            self.StructureToFluidMapper.VectorMap(VectorVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def StructureToFluid_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination, sign_pos):
        (self.StructureToFluidMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    def FluidToStructure_VectorMap(self, VectorVar_Origin, VectorVar_Destination, sign_pos, distributed):
        if distributed:
            VariableRedistributionUtility.DistributePointValues(self.fl_interface,VectorVar_Origin,VAUX_EQ_TRACTION,self.tol,self.it_max)
            self.FluidToStructureMapper.VectorMap(VAUX_EQ_TRACTION, VAUX_EQ_TRACTION, self.it_max, self.tol, sign_pos)
            VariableRedistributionUtility.ConvertDistributedValuesToPoint(self.str_interface,VAUX_EQ_TRACTION,VectorVar_Destination)
        else:
            self.FluidToStructureMapper.VectorMap(VectorVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

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


class NonConformantTwoFaces_OneSideMap:

    def __init__(self, fluid_model_part_positive, fluid_model_part_negative, structure_model_part,
                 search_radius_factor=2.0, it_max=25, tol=1e-5):

        # Error check
        if fluid_model_part_positive.NumberOfConditions(0) < 1:
            raise ValueError("No conditions found in the positive interface fluid model part, please check that the interface is meshed.")
        if fluid_model_part_negative.NumberOfConditions(0) < 1:
            raise ValueError("No conditions found in the negative interface fluid model part, please check that the interface is meshed.")
        if structure_model_part.NumberOfConditions(0) < 1:
            raise ValueError("No conditions found in the structure interface model part, please check that the interface is meshed.")

        self.search_radius_factor = search_radius_factor
        self.it_max = it_max
        self.tol = tol

        self.fluid_model_part_positive = fluid_model_part_positive
        self.fluid_model_part_negative = fluid_model_part_negative
        self.structure_model_part = structure_model_part

        self.PositiveFluidToStructureMapper = AdvancedNMPointsMapper(fluid_model_part_positive, structure_model_part)
        self.NegativeFluidToStructureMapper = AdvancedNMPointsMapper(fluid_model_part_negative, structure_model_part)
        self.StructureToPositiveFluidMapper = AdvancedNMPointsMapper(structure_model_part, fluid_model_part_positive)
        self.StructureToNegativeFluidMapper = AdvancedNMPointsMapper(structure_model_part, fluid_model_part_negative)
        print("Interface Mappers created.")

        # Neighbour search
        (self.PositiveFluidToStructureMapper).FindNeighbours(self.search_radius_factor)
        (self.NegativeFluidToStructureMapper).FindNeighbours(self.search_radius_factor)
        (self.StructureToPositiveFluidMapper).FindNeighbours(self.search_radius_factor)
        (self.StructureToNegativeFluidMapper).FindNeighbours(self.search_radius_factor)
        print("Neighbours search finished.")

    def RecomputeTransferPairs(self, search_radius_factor=2.0):
        (self.PositiveFluidToStructureMapper).FindNeighbours(search_radius_factor)
        (self.NegativeFluidToStructureMapper).FindNeighbours(search_radius_factor)
        (self.StructureToPositiveFluidMapper).FindNeighbours(search_radius_factor)
        (self.StructureToNegativeFluidMapper).FindNeighbours(search_radius_factor)

    # Standard mappers for positive fluid face
    def StructureToPositiveFluid_VectorMap(self, VectorVar_Origin, VectorVar_Destination, sign_pos, distributed):
        if distributed:
            VariableRedistributionUtility.DistributePointValues(self.structure_model_part,VectorVar_Origin,VAUX_EQ_TRACTION,self.tol,self.it_max)
            self.StructureToPositiveFluidMapper.VectorMap(VAUX_EQ_TRACTION, VAUX_EQ_TRACTION, self.it_max, self.tol, sign_pos)
            VariableRedistributionUtility.ConvertDistributedValuesToPoint(self.positive_fluid,VAUX_EQ_TRACTION,VectorVar_Destination)
        else:
            self.StructureToPositiveFluidMapper.VectorMap(VectorVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def StructureToPositiveFluid_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination, sign_pos):
        (self.StructureToPositiveFluidMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    def PositiveFluidToStructure_VectorMap(self, VectorVar_Origin, VectorVar_Destination, sign_pos, distributed):
        if distributed:
            VariableRedistributionUtility.DistributePointValues(self.fluid_model_part_positive,VectorVar_Origin,VAUX_EQ_TRACTION,self.tol,self.it_max)
            self.PositiveFluidToStructureMapper.VectorMap(VAUX_EQ_TRACTION, VAUX_EQ_TRACTION, self.it_max, self.tol, sign_pos)
            VariableRedistributionUtility.ConvertDistributedValuesToPoint(self.structure_model_part,VAUX_EQ_TRACTION,VectorVar_Destination)
        else:
            self.PositiveFluidToStructureMapper.VectorMap(VectorVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def PositiveFluidToStructure_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination, sign_pos):
        (self.PositiveFluidToStructureMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    # Standard mappers for negative fluid face
    def StructureToNegativeFluid_VectorMap(self, VectorVar_Origin, VectorVar_Destination, sign_pos, distributed):
        if distributed:
            VariableRedistributionUtility.DistributePointValues(self.structure_model_part,VectorVar_Origin,VAUX_EQ_TRACTION,self.tol,self.it_max)
            self.StructureToNegativeFluidMapper.VectorMap(VAUX_EQ_TRACTION, VAUX_EQ_TRACTION, self.it_max, self.tol, sign_pos)
            VariableRedistributionUtility.ConvertDistributedValuesToPoint(self.negative_fluid,VAUX_EQ_TRACTION,VectorVar_Destination)
        else:
            self.StructureToNegativeFluidMapper.VectorMap(VectorVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def StructureToNegativeFluid_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination, sign_pos):
        (self.StructureToNegativeFluidMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    def NegativeFluidToStructure_VectorMap(self, VectorVar_Origin, VectorVar_Destination, sign_pos, distributed):
        if distributed:
            VariableRedistributionUtility.DistributePointValues(self.fluid_model_part_negative,VectorVar_Origin,VAUX_EQ_TRACTION,self.tol,self.it_max)
            self.NegativeFluidToStructureMapper.VectorMap(VAUX_EQ_TRACTION, VAUX_EQ_TRACTION, self.it_max, self.tol, sign_pos)
            VariableRedistributionUtility.ConvertDistributedValuesToPoint(self.structure_model_part,VAUX_EQ_TRACTION,VectorVar_Destination)
        else:
            self.NegativeFluidToStructureMapper.VectorMap(VectorVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def NegativeFluidToStructure_ScalarMap(self, ScalarVar_Origin, ScalarVar_Destination, sign_pos):
        (self.NegativeFluidToStructureMapper).ScalarMap(ScalarVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    # Normal vectors mappers for positive fluid face
    def StructureToPositiveFluid_ScalarToNormalVectorMap(self, ScalarVar_Origin, VectorVar_Destination, sign_pos):
        (self.StructureToPositiveFluidMapper).ScalarToNormalVectorMap(ScalarVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def StructureToPositiveFluid_NormalVectorToScalarMap(self, VectorVar_Origin, ScalarVar_Destination, sign_pos):
        (self.StructureToPositiveFluidMapper).NormalVectorToScalarMap(VectorVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    def PositiveFluidToStructure_ScalarToNormalVectorMap(self, ScalarVar_Origin, VectorVar_Destination, sign_pos):
        (self.PositiveFluidToStructureMapper).ScalarToNormalVectorMap(ScalarVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def PositiveFluidToStructure_NormalVectorToScalarMap(self, VectorVar_Origin, ScalarVar_Destination, sign_pos):
        (self.PositiveFluidToStructureMapper).NormalVectorToScalarMap(VectorVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    # Normal vectors mappers for negative fluid face
    def StructureToNegativeFluid_ScalarToNormalVectorMap(self, ScalarVar_Origin, VectorVar_Destination, sign_pos):
        (self.StructureToNegativeFluidMapper).ScalarToNormalVectorMap(ScalarVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def StructureToNegativeFluid_NormalVectorToScalarMap(self, VectorVar_Origin, ScalarVar_Destination, sign_pos):
        (self.StructureToNegativeFluidMapper).NormalVectorToScalarMap(VectorVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)

    def NegativeFluidToStructure_ScalarToNormalVectorMap(self, ScalarVar_Origin, VectorVar_Destination, sign_pos):
        (self.NegativeFluidToStructureMapper).ScalarToNormalVectorMap(ScalarVar_Origin, VectorVar_Destination, self.it_max, self.tol, sign_pos)

    def NegativeFluidToStructure_NormalVectorToScalarMap(self, VectorVar_Origin, ScalarVar_Destination, sign_pos):
        (self.NegativeFluidToStructureMapper).NormalVectorToScalarMap(VectorVar_Origin, ScalarVar_Destination, self.it_max, self.tol, sign_pos)
