import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.mpm_analysis import MpmAnalysis


class TestBodyForceInterpolationMPM(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.work_folder = "test_interpolate_body_force"
        self.parameters_file = "ProjectParameters.json"

    def test_body_force_interpolation(self):

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            with open(self.parameters_file, 'r') as parameter_file:
                project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

            simulation = BodyForceTestSimulation(self.model, project_parameters)
            simulation.Run()

        # --- Get MP model part ---
        mp_model_part = self.model.GetModelPart("MPM_Material")

        self.assertEqual(len(mp_model_part.Elements), 1)

        # Take first MP
        mp = next(iter(mp_model_part.Elements))
        variable_disp = KratosMultiphysics.KratosGlobals.GetVariable( "MP_DISPLACEMENT" )
        variable_acce = KratosMultiphysics.KratosGlobals.GetVariable( "MP_VOLUME_ACCELERATION" )
        displacement = mp.CalculateOnIntegrationPoints(variable_disp,mp_model_part.ProcessInfo)[0]
        acceleration = mp.CalculateOnIntegrationPoints(variable_acce,mp_model_part.ProcessInfo)[0]

        print(acceleration)
        print(displacement)


        # Expected (constant field)
        self.assertAlmostEqual(displacement[0], 0.005, places=12)
        self.assertAlmostEqual(displacement[1], -0.0075, places=12)
        self.assertAlmostEqual(displacement[2], 0.0, places=12)


class BodyForceTestSimulation(MpmAnalysis):

    def __init__(self, model, parameters):
        self.model = model
        self.project_parameters = parameters

   
        #self.project_parameters["problem_data"]["end_time"].SetDouble(0.1)
        #self.project_parameters["solver_settings"]["time_stepping"]["time_step"].SetDouble(0.1)

        super().__init__(self.model, self.project_parameters)

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
