import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication
import KratosMultiphysics.IgaApplication
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

import math


class ConvectionDiffusionAnalysisWithFlush(ConvectionDiffusionAnalysis):
    def __init__(self, model, parameters):
        super().__init__(model, parameters)


@KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication")
@KratosUnittest.skipIfApplicationsNotAvailable("IgaApplication")
class testIGASBMLaplacian(KratosUnittest.TestCase):

    def setUp(self):
        self.work_folder = "test_iga_sbm_laplacian_2D"
        self.print_output = False
        self.print_reference_values = False

    def _readAndCustomizeTestSettings(self):
        # Read the simulation settings
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(self.file_name,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
    def _CheckError(self, DiamondSimulation, toll, expected_results, id_nodes):

        main_model_part = DiamondSimulation._GetSolver().main_model_part

        # need to change this (?)
        i = 0
        for id_node in id_nodes:
            solution_at_control_point = main_model_part.GetNode(id_node).GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            if (abs(solution_at_control_point - expected_results[i]) > toll):
                print('Error in test .....')
                exit()
            i=i+1


    def tearDown(self):
        if not self.print_output:
            with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
                KratosUtilities.DeleteFileIfExisting("square.post.bin")

    # 2D example with dimond hole
    def testIGASBMLaplacianDiamond(self):
        self.file_name = "ProjectParameters.json"

        self._readAndCustomizeTestSettings()

        # test solver without distance modification and without MLS constraints
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            diamond_simulation = ConvectionDiffusionAnalysisWithFlush(self.model, self.parameters)
            diamond_simulation.Run()

            # Check three results in 3 control points
            toll = 1e-10
            expected_results = [2.8860512930390803, 0.6318146904973432, 2.7991964140515844]
            id_nodes = [4216, 4185, 4036]

            # Check error simulation
            self._CheckError(diamond_simulation, toll, expected_results, id_nodes)
    
    # 2D example with dimond hole and external boundary
    def testIGASBMLaplacianDiamondPlusExternal(self):
        self.file_name = "ProjectParameters_plus_external.json"

        self._readAndCustomizeTestSettings()

        # test solver without distance modification and without MLS constraints
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            diamondPlusExternal_simulation = ConvectionDiffusionAnalysisWithFlush(self.model, self.parameters)
            diamondPlusExternal_simulation.Run()

            # Check three results in 3 control points
            toll = 1e-10
            expected_results = [0.14098261189579964, 0.1170606198299867, 0.23531882700536202]
            id_nodes = [600, 700, 800]

            # Check error simulation
            self._CheckError(diamondPlusExternal_simulation, toll, expected_results, id_nodes)



if __name__ == "__main__":
    KratosUnittest.main()
