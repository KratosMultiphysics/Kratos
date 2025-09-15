import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utilities

from KratosMultiphysics.StructuralMechanicsApplication.rve_analysis import RVEAnalysis

class TestRVESimplestTest(KratosUnittest.TestCase):

    def test_rve_computation_block_version_2d(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("rve_test/smallest_rve_test_parameters.json", 'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            parameters["solver_settings"]["domain_size"].SetInt(2)
            parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("rve_test/smallest_rve_test_2D")
            parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString("rve_test/smallest_rve_test_materials_2D.json")
            parameters["solver_settings"]["builder_and_solver_settings"]["use_block_builder"].SetBool(True)
            parameters["solver_settings"]["multi_point_constraints_used"].SetBool(True)

            self._aux_rve_computation(parameters)

    def test_rve_computation_block_version_3d(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("rve_test/smallest_rve_test_parameters.json", 'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            parameters["solver_settings"]["domain_size"].SetInt(3)
            parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("rve_test/smallest_rve_test_3D")
            parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString("rve_test/smallest_rve_test_materials_3D.json")
            parameters["solver_settings"]["builder_and_solver_settings"]["use_block_builder"].SetBool(True)
            parameters["solver_settings"]["multi_point_constraints_used"].SetBool(True)

            self._aux_rve_computation(parameters)

    # FIXME: Commenting until fixed random fail on elimination B&S
    #def test_rve_computation_elimination_version(self):
        ## Within this location context:
        #with KratosUnittest.WorkFolderScope(".", __file__):
            #with open("rve_test/smallest_rve_test_parameters.json", 'r') as parameter_file:
                #parameters = KratosMultiphysics.Parameters(parameter_file.read())

            #parameters["solver_settings"]["block_builder"].SetBool(False)
            #parameters["solver_settings"]["multi_point_constraints_used"].SetBool(True)

            #self._aux_rve_computation(parameters)

    def _aux_rve_computation(self, parameters):

        domain_size = parameters["solver_settings"]["domain_size"].GetInt()
        if parameters["rve_settings"]["print_rve_post"].GetBool():
            output_settings = KratosMultiphysics.Parameters(R'''[{
                "python_module" : "gid_output_process",
                "kratos_module" : "KratosMultiphysics",
                "process_name"  : "GiDOutputProcess",
                "help"          : "This process writes postprocessing files for GiD",
                "Parameters"    : {
                    "model_part_name"        : "Structure",
                    "output_name"            : "TO_BE_SET",
                    "postprocess_parameters" : {
                        "result_file_configuration" : {
                            "gidpost_flags"       : {
                                "GiDPostMode"           : "GiD_PostBinary",
                                "WriteDeformedMeshFlag" : "WriteDeformed",
                                "WriteConditionsFlag"   : "WriteConditions",
                                "MultiFileFlag"         : "SingleFile"
                            },
                            "file_label"          : "step",
                            "output_control_type" : "step",
                            "output_interval"     : 1,
                            "body_output"         : true,
                            "node_output"         : false,
                            "skin_output"         : false,
                            "plane_output"        : [],
                            "nodal_results"       : ["DISPLACEMENT","REACTION"],
                            "gauss_point_results" : ["PK2_STRESS_TENSOR"]
                        },
                        "point_data_configuration"  : []
                    }
                }
            }]''')
            output_settings[0]["Parameters"]["output_name"].SetString("rve_test/smallest_rve_test_{}D".format(domain_size))
            parameters["output_processes"].AddValue("gid_output", output_settings)

        model = KratosMultiphysics.Model()
        simulation = RVEAnalysis(model, parameters)
        simulation.Run()

        model_part = model["Structure"]

        # Compare C
        Canalytic, l, G = self._GetAnalyticElasticityTensor(domain_size)
        Cestimated = model_part.GetValue(StructuralMechanicsApplication.ELASTICITY_TENSOR)

        for i in range(0, Cestimated.Size1()):
            for j in range(0, Cestimated.Size2()):
                self.assertAlmostEqual(
                    abs(Cestimated[i, j] - Canalytic[i, j])/(l+2*G), 0.0, 5)

        if not parameters["rve_settings"]["print_rve_post"].GetBool():
            kratos_utilities.DeleteFileIfExisting("rve_elasticity_tensor_{}D.txt".format(domain_size))

    def _GetAnalyticElasticityTensor(self, domain_size):
        E = 1e6
        nu = 0.3
        l = E*nu/((1+nu)*(1-2*nu))
        G = E/(2.0*(1.0+nu))

        if domain_size == 3:
            Canalytic = KratosMultiphysics.Matrix(6, 6)
            Canalytic.fill(0.0)

            Canalytic[0, 0] = l+2*G
            Canalytic[0, 1] = l
            Canalytic[0, 2] = l

            Canalytic[1, 0] = l
            Canalytic[1, 1] = l+2*G
            Canalytic[1, 2] = l

            Canalytic[2, 0] = l
            Canalytic[2, 1] = l
            Canalytic[2, 2] = l+2*G

            Canalytic[3, 3] = G
            Canalytic[4, 4] = G
            Canalytic[5, 5] = G
        else:
            Canalytic = KratosMultiphysics.Matrix(3, 3)
            Canalytic.fill(0.0)

            Canalytic[0, 0] = l+2*G
            Canalytic[0, 1] = l

            Canalytic[1, 0] = l
            Canalytic[1, 1] = l+2*G

            Canalytic[2, 2] = G

        return Canalytic, l, G


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
