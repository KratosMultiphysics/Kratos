import os
import KratosMultiphysics as Kratos
from KratosMultiphysics import Logger
Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.DEMApplication.DEM_analysis_stage

import auxiliary_functions_for_tests

this_working_dir_backup = os.getcwd()

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class PropertiesMeasureUtilityTestSolution(KratosMultiphysics.DEMApplication.DEM_analysis_stage.DEMAnalysisStage, KratosUnittest.TestCase):

    @classmethod
    def GetMainPath(self):
        return os.path.join(os.path.dirname(os.path.realpath(__file__)), "properties_measure_utility_tests_files")

    def GetProblemNameWithPath(self):
        return os.path.join(self.main_path, self.DEM_parameters["problem_name"].GetString())

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        tolerance = 1e-8

        if self.time >= 100e-8 and self.time < 102e-8:

            self._GetSolver().PrepareContactElementsForPrinting()
            
            #averaged stress tensor for the whole packing
            stress_tensor = self.MeasureSphereForGettingGlobalStressTensor()
            mean_stress = (stress_tensor[0][0]+stress_tensor[1][1]+stress_tensor[2][2])/3
            expected_value_mean_stress = 25.108942763374966
            self.assertAlmostEqual(mean_stress, expected_value_mean_stress, delta=tolerance)

            #define the measured region
            center_x, center_y, center_z = 0.0, 0.0, 0.0
            side_length = 0.00134
            
            #porosity
            measured_porosity = self.MeasureSphereForGettingPackingProperties((side_length/2), center_x, center_y, center_z, 'porosity')
            expected_value_porosity = 0.40930151899675904
            self.assertAlmostEqual(measured_porosity, expected_value_porosity, delta=tolerance)

            #averaged coordination number
            measured_averaged_coordination_number = self.MeasureSphereForGettingPackingProperties((side_length/2), center_x, center_y, center_z, 'averaged_coordination_number')
            expected_value_averaged_coordination_number = 4.2594936708860756
            self.assertAlmostEqual(measured_averaged_coordination_number, expected_value_averaged_coordination_number, delta=tolerance)

            #fabric tensor
            eigenvalues, second_invariant_of_deviatoric_tensor, measured_fabric_tensor = self.MeasureSphereForGettingPackingProperties((side_length/2), center_x, center_y, center_z, 'fabric_tensor')
            expected_value_Eigenvalues = [0.34737273267669533, 0.3348419604048057, 0.31778530691849943]
            expected_value_second_invariant_of_deviatoric_tensor = 0.05940514637933302
            expected_value_fabric_tensor_0 = [ 0.33633556, -0.00391174,  0.00121594]
            expected_value_fabric_tensor_1 = [-0.00391174,  0.34145668, -0.01021716]
            expected_value_fabric_tensor_2 = [ 0.00121594, -0.01021716,  0.32220776]
            self.assertAlmostEqual(eigenvalues[0], expected_value_Eigenvalues[0], delta=tolerance)
            self.assertAlmostEqual(eigenvalues[1], expected_value_Eigenvalues[1], delta=tolerance)
            self.assertAlmostEqual(eigenvalues[2], expected_value_Eigenvalues[2], delta=tolerance)
            self.assertAlmostEqual(second_invariant_of_deviatoric_tensor, expected_value_second_invariant_of_deviatoric_tensor, delta=tolerance)
            self.assertAlmostEqual(measured_fabric_tensor[0][0], expected_value_fabric_tensor_0[0], delta=tolerance)
            self.assertAlmostEqual(measured_fabric_tensor[0][1], expected_value_fabric_tensor_0[1], delta=tolerance)
            self.assertAlmostEqual(measured_fabric_tensor[0][2], expected_value_fabric_tensor_0[2], delta=tolerance)
            self.assertAlmostEqual(measured_fabric_tensor[1][0], expected_value_fabric_tensor_1[0], delta=tolerance)
            self.assertAlmostEqual(measured_fabric_tensor[1][1], expected_value_fabric_tensor_1[1], delta=tolerance)
            self.assertAlmostEqual(measured_fabric_tensor[1][2], expected_value_fabric_tensor_1[2], delta=tolerance)
            self.assertAlmostEqual(measured_fabric_tensor[2][0], expected_value_fabric_tensor_2[0], delta=tolerance)
            self.assertAlmostEqual(measured_fabric_tensor[2][1], expected_value_fabric_tensor_2[1], delta=tolerance)
            self.assertAlmostEqual(measured_fabric_tensor[2][2], expected_value_fabric_tensor_2[2], delta=tolerance)

            #conductivity tensor
            particle_number_inside, measured_non_homogenized_conductivity_tensor_diag, conductivity_tensor_trace, angles_xy, angles_xz, angles_yz = self.MeasureSphereForGettingPackingProperties((side_length/2), center_x, center_y, center_z, 'conductivity_tensor')
            expected_value_particle_number_inside = 158
            expected_value_non_homogenized_conductivity_tensor_diag = [0.0001868674957271866, 0.00019180225820702713, 0.00016420820441842714]
            expected_value_conductivity_tensor_trace = 0.0001809593194508803
            self.assertAlmostEqual(particle_number_inside, expected_value_particle_number_inside, delta=tolerance)
            self.assertAlmostEqual(measured_non_homogenized_conductivity_tensor_diag[0], expected_value_non_homogenized_conductivity_tensor_diag[0], delta=tolerance)
            self.assertAlmostEqual(measured_non_homogenized_conductivity_tensor_diag[1], expected_value_non_homogenized_conductivity_tensor_diag[1], delta=tolerance)
            self.assertAlmostEqual(measured_non_homogenized_conductivity_tensor_diag[2], expected_value_non_homogenized_conductivity_tensor_diag[2], delta=tolerance)
            self.assertAlmostEqual(conductivity_tensor_trace, expected_value_conductivity_tensor_trace, delta=tolerance)

            #stress tensor
            stress_tensor = self.MeasureSphereForGettingPackingProperties((side_length/2), center_x, center_y, center_z, 'stress_tensor')
            expected_value_stress_tensor_0 = [27.41480411, -1.8162967,   0.33081428]
            expected_value_stress_tensor_1 = [-1.74243884, 30.82999517, -3.33601006]
            expected_value_stress_tensor_2 = [ 0.41156874, -3.38081926, 25.73509936]
            self.assertAlmostEqual(stress_tensor[0][0], expected_value_stress_tensor_0[0], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[0][1], expected_value_stress_tensor_0[1], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[0][2], expected_value_stress_tensor_0[2], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[1][0], expected_value_stress_tensor_1[0], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[1][1], expected_value_stress_tensor_1[1], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[1][2], expected_value_stress_tensor_1[2], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[2][0], expected_value_stress_tensor_2[0], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[2][1], expected_value_stress_tensor_2[1], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[2][2], expected_value_stress_tensor_2[2], delta=tolerance)

            #unbalanced force
            measured_unbalanced_force = self.MeasureSphereForGettingPackingProperties((side_length/2), center_x, center_y, center_z, 'unbalanced_force')
            expected_value_unbalanced_force = 0.10717813239546194
            self.assertAlmostEqual(measured_unbalanced_force, expected_value_unbalanced_force, delta=tolerance)

            #stress tensor
            stress_tensor = self.MeasureCubicForGettingPackingProperties((side_length/2), center_x, center_y, center_z, 'stress_tensor')
            expected_value_stress_tensor_0 = [35.79565787, -5.48469616,  4.9171505 ]
            expected_value_stress_tensor_1 = [-5.41156554, 34.72392181, -3.67371185]
            expected_value_stress_tensor_2 = [ 4.86777429, -3.59185829, 27.24868094]
            self.assertAlmostEqual(stress_tensor[0][0], expected_value_stress_tensor_0[0], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[0][1], expected_value_stress_tensor_0[1], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[0][2], expected_value_stress_tensor_0[2], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[1][0], expected_value_stress_tensor_1[0], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[1][1], expected_value_stress_tensor_1[1], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[1][2], expected_value_stress_tensor_1[2], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[2][0], expected_value_stress_tensor_2[0], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[2][1], expected_value_stress_tensor_2[1], delta=tolerance)
            self.assertAlmostEqual(stress_tensor[2][2], expected_value_stress_tensor_2[2], delta=tolerance)

    def Finalize(self):
        self.procedures.RemoveFoldersWithResults(str(self.main_path), str(self.problem_name), '')
        super().Finalize()

class TestPropertiesMeasureUtility(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def test_PropertiesMeasureUtility(self):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "properties_measure_utility_tests_files")
        parameters_file_name = os.path.join(path, "ProjectParametersDEM.json")
        model = Kratos.Model()
        auxiliary_functions_for_tests.CreateAndRunStageInSelectedNumberOfOpenMPThreads(PropertiesMeasureUtilityTestSolution, model, parameters_file_name, 1)

if __name__ == "__main__":
    Kratos.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()
