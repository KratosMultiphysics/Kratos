from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utilities

from KratosMultiphysics.StructuralMechanicsApplication.rve_analysis import RVEAnalysis

class TestRVESimplestTest(KratosUnittest.TestCase):

    def test_rve_computation_block_version(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("rve_test/smallest_rve_test_parameters.json", 'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            parameters["solver_settings"]["block_builder"].SetBool(True)
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

        if parameters["rve_settings"]["print_rve_post"].GetBool():
            parameters["output_processes"].AddValue("gid_output", KratosMultiphysics.Parameters(R'''[{
                "python_module" : "gid_output_process",
                "kratos_module" : "KratosMultiphysics",
                "process_name"  : "GiDOutputProcess",
                "help"          : "This process writes postprocessing files for GiD",
                "Parameters"    : {
                    "model_part_name"        : "Structure",
                    "output_name"            : "rve_test/smallest_rve_test",
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
                            "output_frequency"    : 1,
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
            }]'''))

        model = KratosMultiphysics.Model()
        simulation = RVEAnalysis(model, parameters)
        simulation.Run()

        model_part = model["Structure"]

        # Compare C
        Cestimated = model_part.GetValue(StructuralMechanicsApplication.ELASTICITY_TENSOR)

        Canalytic = KratosMultiphysics.Matrix(6, 6)
        Canalytic.fill(0.0)
        E = 1e6
        nu = 0.3
        l = E*nu/((1+nu)*(1-2*nu))
        G = E/(2.0*(1.0+nu))
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

        for i in range(0, Cestimated.Size1()):
            for j in range(0, Cestimated.Size2()):
                self.assertAlmostEqual(
                    abs(Cestimated[i, j] - Canalytic[i, j])/(l+2*G), 0.0, 5)

        if not parameters["rve_settings"]["print_rve_post"].GetBool():
            kratos_utilities.DeleteFileIfExisting("rve_elasticity_tensor.txt")


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
