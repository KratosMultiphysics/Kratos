from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

import math

class BFECCConvectionTest(UnitTest.TestCase):
    def setUp(self):
        self.dt = 0.05
        self.end_time = 10.0
        self.mesh_divisions = 50
        self.bfecc_substeps = 10.0
        self.print_output = False
        self.check_tolerance = 1.0e-8
        self.print_reference_values = False
        self.work_folder = "BFECCConvectionTest"
        self.reference_file = "bfecc_convection_test"

    def runTest(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            # Create the test model part
            self.model = KratosMultiphysics.Model()
            self.main_model_part = self.model.CreateModelPart("MainModelPart")
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

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
            mesh_parameters.AddEmptyValue("number_of_divisions").SetInt(self.mesh_divisions)

            KratosMultiphysics.StructuredMeshGeneratorProcess(problem_domain, self.main_model_part, mesh_parameters).Execute()

            # Set buffer size
            self.main_model_part.SetBufferSize(2)

            # Set the convection velocity and the initial temperature fields
            x_c = 0.2
            y_c = 0.2
            for node in self.main_model_part.Nodes:
                x_local = node.X - 0.5
                y_local = node.Y - 0.5
                r = math.sqrt(x_local**2 + y_local**2)
                if(r < 0.45):
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, [-y_local, x_local, 0.0])
                    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0, math.sqrt((x_local - x_c)**2 + (y_local - y_c)**2) - 0.1)

            # Create the search structure
            locator = KratosMultiphysics.BinBasedFastPointLocator2D(self.main_model_part)
            locator.UpdateSearchDatabase()

            # Construct the BFECC utility
            bfecc_utility = ConvectionDiffusionApplication.BFECCConvection2D(locator)
            self.main_model_part.CloneTimeStep(0.0)

            # Set output
            if (self.print_output):
                gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
                multifile = KratosMultiphysics.MultiFileFlag.SingleFile
                deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
                write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly
                gid_io = KratosMultiphysics.GidIO(self.reference_file, gid_mode, multifile, deformed_mesh_flag, write_conditions)

                mesh_name = 0.0
                gid_io.InitializeMesh(mesh_name)
                gid_io.WriteMesh(self.main_model_part.GetMesh())
                gid_io.FinalizeMesh()
                gid_io.InitializeResults(mesh_name,(self.main_model_part).GetMesh())

            # Advance in time and convect
            t = 0.0
            while(t < self.end_time):
                t += self.dt
                self.main_model_part.CloneTimeStep(t)

                bfecc_utility.BFECCconvect(
                    self.main_model_part,
                    KratosMultiphysics.TEMPERATURE,
                    KratosMultiphysics.VELOCITY,
                    self.bfecc_substeps)

                if (self.print_output):
                    gid_io.WriteNodalResults(KratosMultiphysics.VELOCITY, self.main_model_part.Nodes, t, 0)
                    gid_io.WriteNodalResults(KratosMultiphysics.TEMPERATURE, self.main_model_part.Nodes, t, 0)

            if (self.print_output):
                gid_io.FinalizeResults()

    def checkResults(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            if self.print_reference_values:
                with open(self.reference_file + '.csv','w') as reference_file:
                    reference_file.write("#ID, TEMPERATURE\n")
                    for node in self.main_model_part.Nodes:
                        reference_file.write("{0}, {1}\n".format(node.Id, node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0)))
                    reference_file.close()
            else:
                with open(self.reference_file + '_results' + '.csv','w') as reference_file:
                    reference_file.write("#ID, TEMPERATURE\n")
                    for node in self.main_model_part.Nodes:
                        reference_file.write("{0}, {1}\n".format(node.Id, node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0)))
                    reference_file.close()

                compare_files_settings = KratosMultiphysics.Parameters(r'''{
                    "reference_file_name"   : "bfecc_convection_test.csv",
                    "output_file_name"      : "bfecc_convection_test_results.csv",
                    "remove_output_file"    : true,
                    "comparison_type"       : "deterministic",
                    "tolerance"             : 1e-2,
                    "relative_tolerance"    : 1e-5,
                    "dimension"             : 2
                }''')
                CompareTwoFilesCheckProcess(compare_files_settings).Execute()

    def testBFECCConvection(self):
        self.setUp()
        self.runTest()
        self.checkResults()
        self.tearDown()

if __name__ == '__main__':
    test = BFECCConvectionTest()
    test.setUp()
    test.testBFECCConvection()
    test.tearDown()
