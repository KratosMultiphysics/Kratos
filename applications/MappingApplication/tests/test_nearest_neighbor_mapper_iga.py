import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication # registering the mappers
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import ParallelEnvironment, IsDistributedRun
import basic_mapper_tests
data_comm = KM.Testing.GetDefaultDataCommunicator()
if data_comm.IsDistributed():
    from KratosMultiphysics.MappingApplication import MPIExtension as MappingMPIExtension
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Additional imports
import mapper_test_case
import os
from pathlib import Path


def GetFilePath(file_name):
    return Path(__file__).resolve().parent / "mdpa_files" / file_name

class BasicTestsLineMappingIGAFEM(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        if not CheckIfApplicationsAvailable("IgaApplication"):
            raise KratosUnittest.SkipTest("The IgaApplication is not available!")
        import KratosMultiphysics.IgaApplication as Iga
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor_iga",
            "interface_submodel_part_origin": "neumann_boundary_iga",
            "interface_submodel_part_destination": "dirichlet_boundary_fem",
            "is_origin_iga": true,
            "echo_level" : 0
        }""")
        cls.setUpMapper(mapper_params)

    @classmethod
    def setUpMapper(cls, mapper_parameters, switch_sides=False):
        super().setUpModelParts("origin_line_interface_iga", "destination_line_interface_fem")

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()
        cls.mapper_parameters = mapper_parameters.Clone()

        if mapper_parameters.Has("interface_submodel_part_origin"):
            cls.interface_model_part_origin = cls.model_part_origin.GetSubModelPart(
                mapper_parameters["interface_submodel_part_origin"].GetString())
        else:
            cls.interface_model_part_origin = cls.model_part_origin

        if mapper_parameters.Has("interface_submodel_part_destination"):
            cls.interface_model_part_destination = cls.model_part_destination.GetSubModelPart(
                mapper_parameters["interface_submodel_part_destination"].GetString())
        else:
            cls.interface_model_part_destination = cls.model_part_destination

        cls.mapper = KM.MapperFactory.CreateMapper(
                cls.model_part_origin, cls.model_part_destination, mapper_parameters)
        
        
    @classmethod
    def setUpModelParts(cls, input_file_name_origin, input_file_name_destination):
        cls.input_file_origin      = GetFilePath(input_file_name_origin)
        cls.input_file_destination = GetFilePath(input_file_name_destination)

        cls.current_model = KM.Model()
        cls.model_part_origin = cls.current_model.CreateModelPart("origin")
        cls.model_part_destination = cls.current_model.CreateModelPart("destination")

        # list of variables involved in the Mapper-Tests
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.PRESSURE)
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.FORCE)
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

        cls.model_part_destination.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.VELOCITY)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.REACTION)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.MESH_DISPLACEMENT)

        cls.model_part_origin.ProcessInfo[KM.DOMAIN_SIZE] = 3 # needed for the partitioner!
        cls.model_part_destination.ProcessInfo[KM.DOMAIN_SIZE] = 3 # needed for the partitioner!
        cls.model_part_origin.ProcessInfo[KM.TIME] = 0.0 # needed for the check-processes
        cls.model_part_destination.ProcessInfo[KM.TIME] = 0.0 # needed for the check-processes
        cls.model_part_origin.ProcessInfo[KM.DELTA_TIME] = 1.0 # needed for the check-processes
        cls.model_part_destination.ProcessInfo[KM.DELTA_TIME] = 1.0 # needed for the check-processes

        cls.ReadModelParts()

    @classmethod
    def ReadModelParts(cls):
        # Read the origin input and create the elements and conditions
        KM.CadJsonInput(cls.input_file_origin).ReadModelPart(cls.model_part_origin)
        origin_interface_brep_curve = cls.model_part_origin.GetGeometry(4)
        origin_interface_quadrature_point_geometries = KM.GeometriesVector()
        origin_interface_brep_curve.CreateQuadraturePointGeometries(origin_interface_quadrature_point_geometries, 3)
        origin_interface_sub_model_part = cls.model_part_origin.CreateSubModelPart("neumann_boundary_iga")

        # Create properties for the elements
        if not cls.model_part_origin.HasProperties(1):
            cls.model_part_origin.CreateNewProperties(1)

        shell_properties = cls.model_part_origin.GetProperties()[1]
        condition_id = 1
        for i in range(0, len(origin_interface_quadrature_point_geometries)):
            origin_interface_sub_model_part.CreateNewCondition('LoadCondition', condition_id, origin_interface_quadrature_point_geometries[i], shell_properties)
            condition_id += 1

        for node in cls.model_part_origin.Nodes:
            origin_interface_sub_model_part.AddNode(node)

        # Read the destination model part 
        KM.ModelPartIO(cls.input_file_destination).ReadModelPart(cls.model_part_destination)

    # No inverse mapping for this mapper is defined as the origin must always be IGA
    def test_InverseMap_constant_scalar(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_InverseMap_constant_vector(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_Map_non_constant_scalar(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_InverseMap_non_constant_scalar(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_Map_non_constant_vector(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_InverseMap_non_constant_vector(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_SWAP_SIGN_InverseMap_scalar(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_SWAP_SIGN_InverseMap_vector(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_ADD_VALUES_InverseMap_scalar(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_ADD_VALUES_InverseMap_vector(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_Map_USE_TRANSPOSE_constant_scalar(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_Map_USE_TRANSPOSE_constant_vector(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_InverseMap_constant_scalar_TO_NON_HISTORICAL(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_InverseMap_constant_scalar_FROM_NON_HISTORICAL(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_InverseMap_constant_scalar_both_non_historical(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_Map_USE_TRANSPOSE_constant_scalar_TO_NON_HISTORICAL(self):
        self.skipTest("Not implemented for this mapper")

    def test_Map_USE_TRANSPOSE_constant_scalar_FROM_NON_HISTORICAL(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_Map_USE_TRANSPOSE_constant_scalar_both_non_historical(self):
        self.skipTest("Not implemented for this mapper")

class BasicTestsSurfaceMappingIGAFEM(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        if not CheckIfApplicationsAvailable("IgaApplication"):
            raise KratosUnittest.SkipTest("The IgaApplication is not available!")
        import KratosMultiphysics.IgaApplication as Iga
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor_iga",
            "interface_submodel_part_origin": "neumann_boundary_iga",
            "interface_submodel_part_destination": "dirichlet_boundary_fem",
            "is_origin_iga": true,
            "echo_level" : 0
        }""")
        cls.setUpMapper(mapper_params)

    @classmethod
    def setUpMapper(cls, mapper_parameters, switch_sides=False):
        super().setUpModelParts("origin_surface_iga", "destination_surface_fem")

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()
        cls.mapper_parameters = mapper_parameters.Clone()

        if mapper_parameters.Has("interface_submodel_part_origin"):
            cls.interface_model_part_origin = cls.model_part_origin.GetSubModelPart(
                mapper_parameters["interface_submodel_part_origin"].GetString())
        else:
            cls.interface_model_part_origin = cls.model_part_origin

        if mapper_parameters.Has("interface_submodel_part_destination"):
            cls.interface_model_part_destination = cls.model_part_destination.GetSubModelPart(
                mapper_parameters["interface_submodel_part_destination"].GetString())
        else:
            cls.interface_model_part_destination = cls.model_part_destination

        cls.mapper = KM.MapperFactory.CreateMapper(
                cls.model_part_origin, cls.model_part_destination, mapper_parameters)
        
        
    @classmethod
    def setUpModelParts(cls, input_file_name_origin, input_file_name_destination):
        cls.input_file_origin      = GetFilePath(input_file_name_origin)
        cls.input_file_destination = GetFilePath(input_file_name_destination)

        cls.current_model = KM.Model()
        cls.model_part_origin = cls.current_model.CreateModelPart("origin")
        cls.model_part_destination = cls.current_model.CreateModelPart("destination")

        # list of variables involved in the Mapper-Tests
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.PRESSURE)
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.FORCE)
        cls.model_part_origin.AddNodalSolutionStepVariable(KM.DISPLACEMENT)

        cls.model_part_destination.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.VELOCITY)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.REACTION)
        cls.model_part_destination.AddNodalSolutionStepVariable(KM.MESH_DISPLACEMENT)

        cls.model_part_origin.ProcessInfo[KM.DOMAIN_SIZE] = 3 # needed for the partitioner!
        cls.model_part_destination.ProcessInfo[KM.DOMAIN_SIZE] = 3 # needed for the partitioner!
        cls.model_part_origin.ProcessInfo[KM.TIME] = 0.0 # needed for the check-processes
        cls.model_part_destination.ProcessInfo[KM.TIME] = 0.0 # needed for the check-processes
        cls.model_part_origin.ProcessInfo[KM.DELTA_TIME] = 1.0 # needed for the check-processes
        cls.model_part_destination.ProcessInfo[KM.DELTA_TIME] = 1.0 # needed for the check-processes

        cls.ReadModelParts()

    @classmethod
    def ReadModelParts(cls):
        # Read the origin input and create the elements and conditions
        KM.CadJsonInput(cls.input_file_origin).ReadModelPart(cls.model_part_origin)
        origin_interface_brep_surface = cls.model_part_origin.GetGeometry(2)
        origin_interface_quadrature_point_geometries = KM.GeometriesVector()
        origin_interface_brep_surface.CreateQuadraturePointGeometries(origin_interface_quadrature_point_geometries, 3)
        origin_interface_sub_model_part = cls.model_part_origin.CreateSubModelPart("neumann_boundary_iga")

        # Create properties for the elements
        if not cls.model_part_origin.HasProperties(1):
            cls.model_part_origin.CreateNewProperties(1)

        shell_properties = cls.model_part_origin.GetProperties()[1]
        condition_id = 1
        for i in range(0, len(origin_interface_quadrature_point_geometries)):
            origin_interface_sub_model_part.CreateNewCondition('LoadCondition', condition_id, origin_interface_quadrature_point_geometries[i], shell_properties)
            condition_id += 1

        for node in cls.model_part_origin.Nodes:
            origin_interface_sub_model_part.AddNode(node)

        # Read the destination model part 
        KM.ModelPartIO(cls.input_file_destination).ReadModelPart(cls.model_part_destination)

    # No inverse mapping for this mapper is defined as the origin must always be IGA
    def test_InverseMap_constant_scalar(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_InverseMap_constant_vector(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_Map_non_constant_scalar(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_InverseMap_non_constant_scalar(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_Map_non_constant_vector(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_InverseMap_non_constant_vector(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_SWAP_SIGN_InverseMap_scalar(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_SWAP_SIGN_InverseMap_vector(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_ADD_VALUES_InverseMap_scalar(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_ADD_VALUES_InverseMap_vector(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_Map_USE_TRANSPOSE_constant_scalar(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_Map_USE_TRANSPOSE_constant_vector(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_InverseMap_constant_scalar_TO_NON_HISTORICAL(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_InverseMap_constant_scalar_FROM_NON_HISTORICAL(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_InverseMap_constant_scalar_both_non_historical(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_Map_USE_TRANSPOSE_constant_scalar_TO_NON_HISTORICAL(self):
        self.skipTest("Not implemented for this mapper")

    def test_Map_USE_TRANSPOSE_constant_scalar_FROM_NON_HISTORICAL(self):
        self.skipTest("Not implemented for this mapper")
    
    def test_Map_USE_TRANSPOSE_constant_scalar_both_non_historical(self):
        self.skipTest("Not implemented for this mapper")

def SetHistoricalNonUniformSolutionScalar(conditions, variable):
    for condition in conditions:
        geometry = condition.GetGeometry()
        val = 12*sin(geometry.Center().X) + geometry.Center().Y*15 + 22*geometry.Center().Z
        condition.SetSolutionStepValue(variable, val)

def SetHistoricalNonUniformSolutionVector(nodes, variable):
    for node in nodes:
        val_1 = 12*sin(node.X0) + node.Y0*15
        val_2 = 33*cos(node.X0) + node.Y0*5
        val_3 = 12*sin(node.Y0) + 22*node.X0
        node.SetSolutionStepValue(variable, KM.Vector([val_1, val_2, val_3]))

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()
