import KratosMultiphysics as kratos
import KratosMultiphysics.ConvectionDiffusionApplication as convdiff

import KratosMultiphysics.KratosUnittest as unittest

from convection_diffusion_analysis import ConvectionDiffusionAnalysis

class AdjointHeatDiffusionTest(unittest.TestCase):

    def setUp(self):
        self.input_file = "adjoint_test"
        #self.reference_file = "reference10_qasgs"
        self.work_folder = "adjoint_diffusion_test"
        self.primal_parameter_file_name = "ProjectParametersPrimal.json"
        self.adjoint_parameter_file_name = "ProjectParametersAdjoint.json"

        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False

    def testAdjointHeatDiffusion(self):
        model = kratos.Model()
        settings = kratos.Parameters(r'''{}''')

        with unittest.WorkFolderScope(self.work_folder, __file__):
            with open(self.primal_parameter_file_name,'r') as primal_parameter_file:
                settings.AddValue("primal_settings", kratos.Parameters(primal_parameter_file.read()))

            primal_analysis = ConvectionDiffusionAnalysis(model,settings["primal_settings"])
            primal_analysis.Run()

            with open(self.adjoint_parameter_file_name,'r') as adjoint_parameter_file:
                settings.AddValue("adjoint_settings", kratos.Parameters(adjoint_parameter_file.read()))

            adjoint_analysis = ConvectionDiffusionAnalysis(model,settings["adjoint_settings"])
            adjoint_analysis.Run()

if __name__ == '__main__':
    unittest.main()