import KratosMultiphysics
import KratosMultiphysics.MappingApplication # registers the "dual_mortar" mapper
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics import from_json_check_result_process
from KratosMultiphysics.MappingApplication import python_mapper_factory

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestDualMortarMapper(KratosUnittest.TestCase):
    """Tests for the mortar mapping.

    The mapping is exercised both through the legacy ``KratosMultiphysics.SimpleMortarMapperProcess``
    (the self-contained in-core operator) and, for the explicit cases, directly through the
    standalone ``dual_mortar`` Mapper (``DualMortarMapper``) created via the MapperFactory. Both
    share the same operator and must produce identical results.
    """

    def setUp(self):
        pass

    def __base_test_mapping(self, input_filename, num_nodes, master_num_nodes, pure_implicit, inverted, discontinuous, origin_are_conditions, destination_are_conditions, consider_tessellation):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        self.model = KratosMultiphysics.Model()

        self.main_model_part = self.model.CreateModelPart("Main",2)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        self.main_model_part.CloneTimeStep(1.01)

        KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)

        if inverted:
            self.model_part_slave = self.main_model_part.GetSubModelPart("Parts_Parts_Auto2")
            self.model_part_master = self.main_model_part.GetSubModelPart("Parts_Parts_Auto1")
        else:
            self.model_part_slave = self.main_model_part.GetSubModelPart("Parts_Parts_Auto1")
            self.model_part_master = self.main_model_part.GetSubModelPart("Parts_Parts_Auto2")

        for node in self.model_part_master.Nodes:
            x = node.X
            y = node.Y
            z = node.Z
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, z)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, x)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, y)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, z)
        del(node)

        self.map_parameters = KratosMultiphysics.Parameters("""
        {
            "echo_level"                       : 0,
            "consider_tessellation"            : false,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2,
            "origin_variable"                  : "TEMPERATURE",
            "discontinuous_interface"          : false,
            "origin_are_conditions"            : true,
            "destination_are_conditions"       : true
        }
        """)
        self.map_parameters["discontinuous_interface"].SetBool(discontinuous)
        self.map_parameters["origin_are_conditions"].SetBool(origin_are_conditions)
        self.map_parameters["destination_are_conditions"].SetBool(destination_are_conditions)
        self.map_parameters["consider_tessellation"].SetBool(consider_tessellation)
        self.pure_implicit = pure_implicit

    def __execute_with_simple_mortar_mapper_process(self):
        # The legacy entry point: the self-contained in-core mortar operator (independent of the
        # MappingApplication). Supports both the explicit and the pure implicit (linear solver) modes.
        if self.pure_implicit:
            linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        else:
            linear_solver = None

        self.map_parameters["origin_variable"].SetString("TEMPERATURE")
        KratosMultiphysics.SimpleMortarMapperProcess(self.model_part_master, self.model_part_slave, self.map_parameters, linear_solver).Execute()
        self.map_parameters["origin_variable"].SetString("DISPLACEMENT")
        KratosMultiphysics.SimpleMortarMapperProcess(self.model_part_master, self.model_part_slave, self.map_parameters, linear_solver).Execute()

    def __execute_with_dual_mortar_mapper_factory(self):
        # The MappingApplication-native entry point: the "dual_mortar" Mapper through the MapperFactory.
        mapper_settings = self.map_parameters.Clone()
        mapper_settings.AddEmptyValue("mapper_type").SetString("dual_mortar")
        mapper = python_mapper_factory.CreateMapper(self.model_part_master, self.model_part_slave, mapper_settings)
        mapper.Map(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.TEMPERATURE)
        mapper.Map(KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.DISPLACEMENT)

    def _mapper_tests(self, input_filename, num_nodes, master_num_nodes, pure_implicit = False, inverted = False, discontinuous = False, origin_are_conditions = True, destination_are_conditions = True, consider_tessellation = False, tolerance_factor = 1.0, use_factory = False):

        self.__base_test_mapping(input_filename, num_nodes, master_num_nodes, pure_implicit, inverted, discontinuous, origin_are_conditions, destination_are_conditions, consider_tessellation)

        if use_factory:
            self.__execute_with_dual_mortar_mapper_factory()
        else:
            self.__execute_with_simple_mortar_mapper_process()

        check_parameters = KratosMultiphysics.Parameters("""
        {
            "check_variables"      : ["TEMPERATURE","DISPLACEMENT"],
            "input_file_name"      : "",
            "model_part_name"      : "Main",
            "sub_model_part_name"  : "Parts_Parts_Auto1",
            "tolerance"            : 1e-3,
            "relative_tolerance"   : 1e-6
        }
        """)
        check_parameters["tolerance"].SetDouble(tolerance_factor * check_parameters["tolerance"].GetDouble())
        check_parameters["relative_tolerance"].SetDouble(tolerance_factor * check_parameters["relative_tolerance"].GetDouble())

        if inverted:
            check_parameters["input_file_name"].SetString(input_filename+"_inverted.json")
        else:
            check_parameters["input_file_name"].SetString(input_filename+".json")

        check = from_json_check_result_process.FromJsonCheckResultProcess(self.model, check_parameters)
        check.ExecuteInitialize()
        check.ExecuteBeforeSolutionLoop()
        check.ExecuteFinalizeSolutionStep()

    def test_less_basic_mortar_mapping_triangle_pure_implicit(self):
        input_filename = GetFilePath("mortar_mapper_python_tests/test_integration_several_triangles")
        self._mapper_tests(input_filename, 3, 3, True, False, False, True, True, False)
        self._mapper_tests(input_filename, 3, 3, True, False, False, True, True, True)

    def test_less_basic_mortar_mapping_triangle(self):
        input_filename = GetFilePath("mortar_mapper_python_tests/test_integration_several_triangles")
        self._mapper_tests(input_filename, 3, 3, False, False, False, True, True, False)
        self._mapper_tests(input_filename, 3, 3, False, False, False, True, True, True)
        # Also exercise the dual_mortar Mapper directly through the factory
        self._mapper_tests(input_filename, 3, 3, False, False, False, True, True, False, use_factory=True)

    def test_simple_curvature_mortar_mapping_triangle(self):
        input_filename = GetFilePath("mortar_mapper_python_tests/test_simple_curvature")
        self._mapper_tests(input_filename, 3, 3, False, False, False, True, True, False)
        self._mapper_tests(input_filename, 3, 3, False, False, False, True, True, True)

    def test_mortar_mapping_triangle(self):
        input_filename = GetFilePath("mortar_mapper_python_tests/test_double_curvature_integration_triangle")
        self._mapper_tests(input_filename, 3, 3, False, False, False, True, True, False)
        self._mapper_tests(input_filename, 3, 3, False, False, False, True, True, True)

    def test_mortar_mapping_triangle_discontinous_interface(self):
        input_filename = GetFilePath("mortar_mapper_python_tests/test_double_curvature_integration_triangle_discontinous_interface")
        self._mapper_tests(input_filename, 3, 3, False, False, True, True, True, False)
        self._mapper_tests(input_filename, 3, 3, False, False, True, True, True, True)

    def test_mortar_mapping_quad(self):
        input_filename = GetFilePath("mortar_mapper_python_tests/test_double_curvature_integration_quadrilateral")
        self._mapper_tests(input_filename, 4, 4, False, False, False, True, True, False)
        self._mapper_tests(input_filename, 4, 4, False, False, False, True, True, True, 15.0)

    def test_mortar_mapping_quad_tri(self):
        input_filename = GetFilePath("mortar_mapper_python_tests/test_double_curvature_integration_triangle_quadrilateral")
        self._mapper_tests(input_filename, 4, 3, False, False, False, False, True, False)
        self._mapper_tests(input_filename, 4, 3, False, False, False, False, True, True, 10.0)

    def test_mortar_mapping_tri_quad(self):
        input_filename = GetFilePath("mortar_mapper_python_tests/test_double_curvature_integration_triangle_quadrilateral")
        self._mapper_tests(input_filename, 3, 4, False, True, False, True, False, False)
        self._mapper_tests(input_filename, 3, 4, False, True, False, True, False, True)

if __name__ == '__main__':
    KratosUnittest.main()
