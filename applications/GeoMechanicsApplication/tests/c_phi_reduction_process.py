# from KratosMultiphysics import * as Kratos

import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis

import test_helper

class KratosGeoMechanicsCPhiReductionProcess(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the original solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_c_phi_reduction_process(self):

        # get the parameter file names for all stages
        test_name = 'C-Phi_reduction_process'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        n_stages = 2
        parameter_file_names = [os.path.join(file_path, 'ProjectParameters_stage' + str(i + 1) + '.json') for i in
                            range(n_stages)]
                            
        # change to project directory
        os.chdir(file_path)
        
        #setup stages from parameterfiles
        parameters_stages = [None] * n_stages
        
        for idx, parameter_file_name in enumerate(parameter_file_names):
            with open(parameter_file_name, 'r') as parameter_file:
                parameters_stages[idx] = Kratos.Parameters(parameter_file.read())

        model = Kratos.Model()
        stages = [GeoMechanicsAnalysis(model, stage_parameters) for stage_parameters in parameters_stages]

        # execute the stages
        [stage.Run() for stage in stages]
    
        # read results
        reader = test_helper.GiDOutputFileReader()
        reader.read_output_from(os.path.join(file_path, "stage2.post.res"))
        displacement = test_helper.get_displacement(stages[1])
        self.assertAlmostEqual(-0.001288561508339424, displacement[12][0])
        self.assertAlmostEqual(-0.001688453602754448, displacement[27][0])
        self.assertAlmostEqual(-0.01123562533250766, displacement[150][0])


if __name__ == '__main__':
    KratosUnittest.main()
