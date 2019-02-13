from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utilities

import RVEAnalysis

import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestRVESimplestTest(KratosUnittest.TestCase):

    def test_rve_computation_block_version(self):
        with open(GetFilePath("rve_test/smallest_rve_test_parameters.json"), 'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        parameters["solver_settings"]["block_builder"].SetBool(True)
        parameters["solver_settings"]["multi_point_constraints_used"].SetBool(True)

        self._aux_rve_computation(parameters)

    #@KratosUnittest.skip("Reactions and displacements not 100% accurate")
    def test_rve_computation_elimination_version(self):
        with open(GetFilePath("rve_test/smallest_rve_test_parameters.json"), 'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        parameters["solver_settings"]["block_builder"].SetBool(False)
        parameters["solver_settings"]["multi_point_constraints_used"].SetBool(True)

        self._aux_rve_computation(parameters)

    def _aux_rve_computation(self, parameters):

        model = KratosMultiphysics.Model()
        simulation = RVEAnalysis.RVEAnalysis(model, parameters)
        simulation.Run()

        model_part = model["Structure.computing_domain"]

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
            kratos_utilities.DeleteFileIfExisting("smallest_test.post.bin")
            kratos_utilities.DeleteFileIfExisting("rve_test.post.lst")
            kratos_utilities.DeleteFileIfExisting("rve_elasticity_tensor.txt")


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
