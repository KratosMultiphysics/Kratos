import os
import KratosMultiphysics as Kratos
from Kratos import Logger
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.DEMApplication.DEM_analysis_stage as dem_analysis

# This test consists in a system with a single already existing particle and an inlet that injects a few
# particles during the simulation, which consists in letting the particle fall under gravity.
# The bounding box, which has its bottom placed at z=0 is set to mark the particles to be erased when they
# cross this limit. Depending on the delay imposed on the destruction of the particles after they are marked,
# a different number of particles is recovered at the end of the simulation (more delay should lead
# to equal or greater number of particles at the end).

debug_mode = False


class TestDEMEraseParticles(dem_analysis.DEMAnalysisStage):

    @staticmethod
    def StaticGetMainPath():
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "erase_particles_test_files")

    def GetMainPath(self):
        return self.StaticGetMainPath()

    def Finalize(self):
        self.number_of_particles_by_the_end = len(self.spheres_model_part.GetElements())
        parent_return = super().Finalize()

class TestDEMEraseParticlesWithDelay(KratosUnittest.TestCase):

    def setUp(self):
        self.parameters_file_name = 'ProjectParametersDEMWithNoDelay.json'
        self.path = TestDEMEraseParticles.StaticGetMainPath()
        self.analysis = TestDEMEraseParticles

    def test_erase_particles_no_delay(self):
        project_parameters_file_name = 'ProjectParametersDEMWithNoDelay.json'
        expected_number_of_particles = 27
        self.RunTest(project_parameters_file_name, expected_number_of_particles)

    def test_erase_particles_little_delay(self):
        project_parameters_file_name = 'ProjectParametersDEMWithLittleDelay.json'
        expected_number_of_particles = 27
        self.RunTest(project_parameters_file_name, expected_number_of_particles)

    def test_erase_particles_with_delay(self):
        project_parameters_file_name = 'ProjectParametersDEMWithDelay.json'
        expected_number_of_particles = 29
        self.RunTest(project_parameters_file_name, expected_number_of_particles)

    def RunTest(self, project_parameters_file_name, expected_number_of_particles):
        parameters_file_path = os.path.join(self.path, project_parameters_file_name)
        model = Kratos.Model()

        with open(parameters_file_path, 'r') as parameter_file:
            project_parameters = Kratos.Parameters(parameter_file.read())

        analysis = self.analysis(model, project_parameters)
        analysis.Run()
        self.assertEqual(expected_number_of_particles, analysis.number_of_particles_by_the_end)

if __name__ == "__main__":
    if debug_mode:
        Logger.GetDefaultOutput().SetSeverity(Logger.Severity.INFO)
    else:
        Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)

    KratosUnittest.main()
