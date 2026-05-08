import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.mpm_analysis import MpmAnalysis


class TestCauchyStressProjectionMPM(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KratosMultiphysics.Model()

    def test_cauchy_stress_projection_updated_lagrangian(self):
        self._RunCauchyStressProjectionTest("quadrilateral", False)

    def test_cauchy_stress_projection_updated_lagrangian_up(self):
        self._RunCauchyStressProjectionTest("triangle", True)

    def _RunCauchyStressProjectionTest(self, geometry_type, pressure_dofs):
        project_parameters = KratosMultiphysics.Parameters("""{
            "problem_data" : {
                "problem_name"  : "test_nodal_cauchy_stress",
                "parallel_type" : "OpenMP",
                "echo_level"    : 0,
                "start_time"    : 0.0,
                "end_time"      : 0.1
            },
            "solver_settings" : {
                "time_stepping"                   : { "time_step" : 0.1 },
                "solver_type"                     : "Dynamic",
                "model_part_name"                 : "MPM_Material",
                "domain_size"                     : 2,
                "echo_level"                      : 0,
                "analysis_type"                   : "non_linear",
                "time_integration_method"         : "implicit",
                "scheme_type"                     : "newmark",
                "model_import_settings"           : { "input_type" : "use_input_model_part" },
                "grid_model_import_settings"      : { "input_type" : "use_input_model_part" },
                "compute_reactions"               : false,
                "compute_nodal_cauchy_stress"     : true,
                "convergence_criterion"           : "residual_criterion",
                "displacement_relative_tolerance" : 0.0001,
                "displacement_absolute_tolerance" : 1e-9,
                "residual_relative_tolerance"     : 0.0001,
                "residual_absolute_tolerance"     : 1e-9,
                "max_iteration"                   : 10,
                "pressure_dofs"                   : false,
                "linear_solver_settings"          : { "solver_type" : "LinearSolversApplication.sparse_lu" },
                "auxiliary_variables_list"        : ["NORMAL","IS_STRUCTURE"]
            },
            "processes" : {
                "constraints_process_list"     : [],
                "loads_process_list"           : [],
                "list_other_processes"         : [],
                "initial_conditions_processes" : [],
                "gravity"                      : []
            },
            "output_processes" : {}
        }""")

        project_parameters["solver_settings"]["pressure_dofs"].SetBool(pressure_dofs)

        simulation = CauchyStressTestSimulation(self.model, project_parameters, geometry_type, pressure_dofs)
        simulation.Run()

        grid_model_part = self.model.GetModelPart("Background_Grid")
        expected_nodal_stress = KratosMultiphysics.Vector([10.0, 5.0, 2.0])

        for node in grid_model_part.Nodes:
            nodal_stress = node.GetSolutionStepValue(KratosMPM.NODAL_CAUCHY_STRESS_VECTOR)
            self.assertVectorAlmostEqual(nodal_stress, expected_nodal_stress)


class CauchyStressTestSimulation(MpmAnalysis):

    def __init__(self, model, parameters, geometry_type, pressure_dofs):
        self.model = model
        self.project_parameters = parameters
        self.geometry_type = geometry_type
        self.pressure_dofs = pressure_dofs

        super().__init__(self.model, self.project_parameters)

    def _PrepareModelPart(self):
        self._CreateBackgroundGridModelPart()
        self._CreateInitialMaterialModelPart()

        super()._PrepareModelPart()

    def _CreateBackgroundGridModelPart(self):
        grid_model_part = self.model.GetModelPart("Background_Grid")
        grid_sub_model_part = grid_model_part.CreateSubModelPart("Grid_Auto1")

        self._CreateNodes(grid_sub_model_part)

        grid_properties = grid_sub_model_part.GetProperties()[1]
        if self.geometry_type == "triangle":
            grid_sub_model_part.CreateNewElement("Element2D3N", 1, [3, 4, 1], grid_properties)
            grid_sub_model_part.CreateNewElement("Element2D3N", 2, [4, 2, 1], grid_properties)
            grid_sub_model_part.CreateNewElement("Element2D3N", 3, [5, 6, 3], grid_properties)
            grid_sub_model_part.CreateNewElement("Element2D3N", 4, [6, 4, 3], grid_properties)
        else:
            grid_sub_model_part.CreateNewElement("Element2D4N", 1, [3, 4, 2, 1], grid_properties)
            grid_sub_model_part.CreateNewElement("Element2D4N", 2, [5, 6, 4, 3], grid_properties)

    def _CreateInitialMaterialModelPart(self):
        initial_mesh_model_part = self.model.GetModelPart("Initial_MPM_Material")
        material_sub_model_part = initial_mesh_model_part.CreateSubModelPart("Material_domain_Auto1")

        properties = material_sub_model_part.GetProperties()[1]
        if self.pressure_dofs:
            properties.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, KratosMPM.HyperElasticPlaneStrainUP2DLaw())
        else:
            properties.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, KratosMPM.LinearElasticIsotropicPlaneStrain2DLaw())
        properties.SetValue(KratosMultiphysics.THICKNESS, 1.0)
        properties.SetValue(KratosMPM.MATERIAL_POINTS_PER_ELEMENT, 1)
        properties.SetValue(KratosMultiphysics.DENSITY, 7850.0)
        properties.SetValue(KratosMultiphysics.YOUNG_MODULUS, 206900000000.0)
        properties.SetValue(KratosMultiphysics.POISSON_RATIO, 0.29)

        self._CreateNodes(material_sub_model_part)

        if self.geometry_type == "triangle":
            material_sub_model_part.CreateNewElement("MPMUpdatedLagrangian2D3N", 1, [3, 4, 1], properties)
            material_sub_model_part.CreateNewElement("MPMUpdatedLagrangian2D3N", 2, [4, 2, 1], properties)
            material_sub_model_part.CreateNewElement("MPMUpdatedLagrangian2D3N", 3, [5, 6, 3], properties)
            material_sub_model_part.CreateNewElement("MPMUpdatedLagrangian2D3N", 4, [6, 4, 3], properties)
        else:
            material_sub_model_part.CreateNewElement("MPMUpdatedLagrangian2D4N", 1, [3, 4, 2, 1], properties)
            material_sub_model_part.CreateNewElement("MPMUpdatedLagrangian2D4N", 2, [5, 6, 4, 3], properties)

    @staticmethod
    def _CreateNodes(model_part):
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        model_part.CreateNewNode(2, 0.0, 1.0, 0.0)
        model_part.CreateNewNode(3, 1.0, 0.0, 0.0)
        model_part.CreateNewNode(4, 1.0, 1.0, 0.0)
        model_part.CreateNewNode(5, 2.0, 0.0, 0.0)
        model_part.CreateNewNode(6, 2.0, 1.0, 0.0)

    def ModifyBeforeSolutionLoop(self):
        super().ModifyBeforeSolutionLoop()

        mp_model_part = self.model.GetModelPart("MPM_Material")
        cauchy_stress_vector = [KratosMultiphysics.Vector([10.0, 5.0, 2.0])]

        for mp in mp_model_part.Elements:
            mp.SetValuesOnIntegrationPoints(
                KratosMPM.MP_CAUCHY_STRESS_VECTOR,
                cauchy_stress_vector,
                0,
                mp_model_part.ProcessInfo)


if __name__ == '__main__':
    KratosUnittest.main()
