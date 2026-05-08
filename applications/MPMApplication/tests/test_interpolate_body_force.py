import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.mpm_analysis import MpmAnalysis


class TestBodyForceInterpolationMPM(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KratosMultiphysics.Model()

    def test_body_force_interpolation(self):
        project_parameters = KratosMultiphysics.Parameters("""{
            "problem_data"     : {
                "problem_name"  : "test_body_force",
                "parallel_type" : "OpenMP",
                "echo_level"    : 0,
                "start_time"    : 0.0,
                "end_time"      : 0.1
            },
            "solver_settings"  : {
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
                "convergence_criterion"           : "residual_criterion",
                "displacement_relative_tolerance" : 0.0001,
                "displacement_absolute_tolerance" : 1e-9,
                "residual_relative_tolerance"     : 0.0001,
                "residual_absolute_tolerance"     : 1e-9,
                "max_iteration"                   : 10,
                "pressure_dofs"                   : false,
                "linear_solver_settings"          : { "solver_type" : "LinearSolversApplication.sparse_lu" },
                "auxiliary_variables_list"        : ["NORMAL","IS_STRUCTURE","BODY_FORCE"]
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

        simulation = BodyForceTestSimulation(self.model, project_parameters)
        simulation.Run()

        # --- Get MP model part ---
        mp_model_part = self.model.GetModelPart("MPM_Material")

        self.assertEqual(len(mp_model_part.Elements), 1)

        # Take first MP
        mp = next(iter(mp_model_part.Elements))
        variable_disp = KratosMultiphysics.KratosGlobals.GetVariable( "MP_DISPLACEMENT" )
        variable_body_force = KratosMultiphysics.KratosGlobals.GetVariable( "MP_BODY_FORCE" )
        displacement = mp.CalculateOnIntegrationPoints(variable_disp,mp_model_part.ProcessInfo)[0]
        body_force = mp.CalculateOnIntegrationPoints(variable_body_force,mp_model_part.ProcessInfo)[0]

        #print(body_force)
        #print(displacement)


        # Expected (constant field)
        self.assertAlmostEqual(displacement[0], 0.005, places=12)
        self.assertAlmostEqual(displacement[1], -0.0075, places=12)
        self.assertAlmostEqual(displacement[2], 0.0, places=12)

        self.assertAlmostEqual(body_force[0], 2, places=12)
        self.assertAlmostEqual(body_force[1], -3, places=12)
        self.assertAlmostEqual(body_force[2], 0, places=12)


class BodyForceTestSimulation(MpmAnalysis):

    def __init__(self, model, parameters):
        self.model = model
        self.project_parameters = parameters

        super().__init__(self.model, self.project_parameters)

    def _PrepareModelPart(self):
        self._CreateBackgroundGridModelPart()
        self._CreateInitialMaterialModelPart()

        super()._PrepareModelPart()

    def _CreateBackgroundGridModelPart(self):
        grid_model_part = self.model.GetModelPart("Background_Grid")
        grid_sub_model_part = grid_model_part.CreateSubModelPart("Grid_Auto1")

        grid_sub_model_part.CreateNewNode(1, 3.0, 0.0, 0.0)
        grid_sub_model_part.CreateNewNode(2, 3.0, 1.0, 0.0)
        grid_sub_model_part.CreateNewNode(3, 2.0, 0.0, 0.0)
        grid_sub_model_part.CreateNewNode(4, 2.0, 1.0, 0.0)
        grid_sub_model_part.CreateNewNode(6, 1.0, 0.0, 0.0)
        grid_sub_model_part.CreateNewNode(8, 1.0, 1.0, 0.0)
        grid_sub_model_part.CreateNewNode(9, 0.0, 0.0, 0.0)
        grid_sub_model_part.CreateNewNode(11, 0.0, 1.0, 0.0)

        grid_properties = grid_sub_model_part.GetProperties()[1]
        grid_sub_model_part.CreateNewElement("Element2D4N", 1, [6, 8, 11, 9], grid_properties)
        grid_sub_model_part.CreateNewElement("Element2D4N", 2, [3, 4, 8, 6], grid_properties)
        grid_sub_model_part.CreateNewElement("Element2D4N", 3, [1, 2, 4, 3], grid_properties)

    def _CreateInitialMaterialModelPart(self):
        initial_mesh_model_part = self.model.GetModelPart("Initial_MPM_Material")
        material_sub_model_part = initial_mesh_model_part.CreateSubModelPart("Material_domain_Auto1")

        properties = material_sub_model_part.GetProperties()[1]
        properties.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, KratosMPM.LinearElasticIsotropicPlaneStrain2DLaw())
        properties.SetValue(KratosMultiphysics.THICKNESS, 1.0)
        properties.SetValue(KratosMPM.MATERIAL_POINTS_PER_ELEMENT, 1)
        properties.SetValue(KratosMultiphysics.DENSITY, 7850.0)
        properties.SetValue(KratosMultiphysics.YOUNG_MODULUS, 206900000000.0)
        properties.SetValue(KratosMultiphysics.POISSON_RATIO, 0.29)

        material_sub_model_part.CreateNewNode(5, 1.0, 0.0, 0.0)
        material_sub_model_part.CreateNewNode(7, 1.0, 1.0, 0.0)
        material_sub_model_part.CreateNewNode(10, 0.0, 0.0, 0.0)
        material_sub_model_part.CreateNewNode(12, 0.0, 1.0, 0.0)

        material_sub_model_part.CreateNewElement("MPMUpdatedLagrangian2D4N", 4, [5, 7, 12, 10], properties)

    def ModifyBeforeSolutionLoop(self):
        super().ModifyBeforeSolutionLoop()

        grid_model_part = self.model.GetModelPart("Background_Grid")

        # Apply constant BODY_FORCE
        for node in grid_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_X, 2.0)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y, -3.0)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Z, 0.0)


if __name__ == '__main__':
    KratosUnittest.main()
