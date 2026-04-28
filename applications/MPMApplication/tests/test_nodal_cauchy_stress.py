import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.MPMApplication as KratosMPM
from KratosMultiphysics.MPMApplication.particle_mechanics_analysis import ParticleMechanicsAnalysis


class TestCauchyStressProjectionMPM(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KratosMultiphysics.Model()
        self.work_folder = "test_nodal_cauchy_stress"
        self.parameters_file = "ProjectParameters.json"

    def test_cauchy_stress_projection(self):
    
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
           with open(self.parameters_file, 'r') as parameter_file:
              project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

           simulation = CauchyStressTestSimulation(self.model, project_parameters)
           simulation.Run()

        # --- Get grid model part ---
        grid_model_part = self.model.GetModelPart("Background_Grid")

        self.assertGreater(len(grid_model_part.Nodes), 0)

        #print("---- NODE RESULTS ----")
        for node in grid_model_part.Nodes:
            nodal_stress = node.GetSolutionStepValue(KratosMPM.NODAL_CAUCHY_STRESS_VECTOR)

            print(f"Node {node.Id}: {nodal_stress}")

            self.assertAlmostEqual(nodal_stress[0], 10.0, places=6)
            self.assertAlmostEqual(nodal_stress[1], 5.0, places=6)
            self.assertAlmostEqual(nodal_stress[2], 2.0, places=6)


# ---------------------------------------------------------------------

class CauchyStressTestSimulation(ParticleMechanicsAnalysis):

    def __init__(self, model, parameters):
        self.model = model
        self.project_parameters = parameters

        # Ajuste temporal mínimo para que haya un step real
        self.project_parameters["problem_data"]["end_time"].SetDouble(0.1)
        self.project_parameters["solver_settings"]["time_stepping"]["time_step"].SetDouble(0.1)

        super().__init__(self.model, self.project_parameters)

    # -----------------------------------------------------------------

    def Initialize(self):
        super().Initialize()

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

    def Predict(self):
        super().Predict()

    def SolveSolutionStep(self):
        super().SolveSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

    # -----------------------------------------------------------------

    def ModifyAfterSolverInitialize(self):
        super().ModifyAfterSolverInitialize()

        #print(">>> MODIFY AFTER SOLVER INITIALIZE (IF CALLED)")

        mp_model_part = self.model.GetModelPart("MPM_Material")

        for mp in mp_model_part.Elements:
        
            cauchy_stress_vector = [KratosMultiphysics.Vector([10,5,2])]	

            mp.SetValuesOnIntegrationPoints(KratosMPM.MP_CAUCHY_STRESS_VECTOR, cauchy_stress_vector,0, mp_model_part.ProcessInfo)
            

if __name__ == '__main__':
    KratosUnittest.main()
