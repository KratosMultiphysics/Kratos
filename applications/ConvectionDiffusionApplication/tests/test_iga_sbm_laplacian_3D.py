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
        self.work_folder = "test_iga_sbm_laplacian_3D"
        self.print_output = False
        self.print_reference_values = False

    def _readAndCustomizeTestSettings(self):
        # Read the simulation settings
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(self.file_name,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())
    
    def _CheckError(self, SphereSimulation):
        # Check L2 error via midpoint rule
        toll = 1e-10
        main_model_part = SphereSimulation._GetSolver().main_model_part

        # Check three results in 3 control points
        expected_results = [0.9438305876574609, 0.9437424537274997, 1.0028950383638127]
        id_nodes = [500, 600, 700]

        # need to change this
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

    def testIGASBMLaplacianSphere(self):
        self.file_name = "ProjectParameters.json"

        self._readAndCustomizeTestSettings()

        # test solver without distance modification and without MLS constraints
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            sphere_simulation = ConvectionDiffusionAnalysisWithFlush(self.model, self.parameters)
            sphere_simulation.Run()

            # Check error simulation
            self._CheckError(sphere_simulation)



if __name__ == "__main__":
    KratosUnittest.main()
