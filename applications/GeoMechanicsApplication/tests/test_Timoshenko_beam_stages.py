import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import test_helper

class KratosGeoMechanicsTimoshenkoBeamStagedTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """
    def test_cantilever_Timoshenko3D3N(self):
        """
        Tests deflection of cantilever beam alogn y axis in 2 stages
        stage 1: z direction tip load is applied
        stage 2: z direction load erased, x direction tip load applied
        """

        # calculate deflection for pure bending beam ( the shear areas have been 
        F = 0.1
        EI_33 = 0.083333333 * 1000.0
        EI_22 = 0.008333333 * 1000.0
        L = 10.0

        uz_tip = (F*L*L*L)/(3*EI_33)
        ux_tip = (F*L*L*L)/(3*EI_22)

        # get stages
        test_name = 'timoshenko3D3N_staged'
        project_path = test_helper.get_file_path(os.path.join('.', test_name))
        n_stages = 2

        original_directory = os.getcwd()
        parameter_file_paths = [os.path.join(project_path, f'ProjectParameters_Stage{i+1}.json') for i in range(n_stages)]

        # set stage parameters
        parameters_stages = []
        os.chdir(project_path)
        for parameter_file_path in parameter_file_paths:
            with open(parameter_file_path, 'r') as parameter_file:
                parameters_stages.append(KratosMultiphysics.Parameters(parameter_file.read()))

        model = KratosMultiphysics.Model()
        displacement_stages = []

        # run stages and get results
        for stage_parameters in parameters_stages:
            stage = analysis.GeoMechanicsAnalysis(model, stage_parameters)
            stage.Run()
            displacement_stages.append(test_helper.get_displacement(stage))

        os.chdir(original_directory)

        # Assert
        stage_nr = 1
        self.assertAlmostEqual(displacement_stages[stage_nr-1][10][2], uz_tip, msg = f"u_z at node 11 in stage {stage_nr}")

        stage_nr = 2
        self.assertAlmostEqual(displacement_stages[stage_nr-1][10][2], -uz_tip, msg = f"u_z at node 11 in stage {stage_nr}")
        self.assertAlmostEqual(displacement_stages[stage_nr-1][10][0], ux_tip, msg = f"u_x at node 11 in stage {stage_nr}")

if __name__ == '__main__':
    KratosUnittest.main()
