from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

import math

class BFECCConvectionTest(UnitTest.TestCase):
    def setUp(self):
        self.dt = 0.05
        self.end_time = 10.0
        self.bfecc_substeps = 10.0

    def runTest(self):
        # Create the test model part
        self.model = KratosMultiphysics.Model()
        main_model_part = self.model.CreateModelPart("MainModelPart")
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

        # Generate the problem domain
        problem_domain = KratosMultiphysics.Quadrilateral2D4(
            KratosMultiphysics.Node(1, 0.0, 0.0, 0.0),
            KratosMultiphysics.Node(2, 0.0, 1.0, 0.0),
            KratosMultiphysics.Node(3, 1.0, 1.0, 0.0),
            KratosMultiphysics.Node(4, 1.0, 0.0, 0.0))

        mesh_parameters = KratosMultiphysics.Parameters("{}")
        mesh_parameters.AddEmptyValue("element_name").SetString("Element2D3N")
        mesh_parameters.AddEmptyValue("condition_name").SetString("Condition2D2N")
        mesh_parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
        mesh_parameters.AddEmptyValue("number_of_divisions").SetInt(100)

        KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, main_model_part, mesh_parameters).Execute()

        # Set buffer size
        main_model_part.SetBufferSize(2)

        # Set the convection velocity and the initial temperature fields
        x_c = 0.2
        y_c = 0.2
        for node in main_model_part.Nodes:
            x_local = node.X - 0.5
            y_local = node.Y - 0.5
            r = math.sqrt(x_local**2 + y_local**2)
            if(r < 0.45):
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, [-y_local, x_local, 0.0])
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0, math.sqrt((x_local - x_c)**2 + (y_local - y_c)**2) - 0.1)

        # Create the search structure
        locator = KratosMultiphysics.BinBasedFastPointLocator2D(main_model_part)
        locator.UpdateSearchDatabase()

        # Construct the BFECC utility
        bfecc_utility = ConvectionDiffusionApplication.BFECCConvection2D(locator)
        main_model_part.CloneTimeStep(0.0)

        # Advance in time and convect
        t = 0.0
        while(t < self.end_time):
            t += self.dt
            main_model_part.CloneTimeStep(t)
            bfecc_utility.BFECCconvect(
                main_model_part,
                KratosMultiphysics.TEMPERATURE,
                KratosMultiphysics.VELOCITY,
                self.bfecc_substeps)

    def checkResults(self):
        pass

    def testBFECCConvection(self):
        self.runTest()
        self.checkResults()

    def tearDown(self):
        pass

if __name__ == '__main__':
    test = BFECCConvectionTest()
    test.setUp()
    test.testBFECCConvection()
    test.tearDown()

