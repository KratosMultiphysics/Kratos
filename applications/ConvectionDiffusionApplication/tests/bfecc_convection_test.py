from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

from KratosMultiphysics.json_output_process import JsonOutputProcess
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

import os
import math

# from KratosMultiphysics.gid_output_process import GiDOutputProcess

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def BaseJumpedDistance(x, y, z):
    if (x >= 0.1 and x <= 0.3):
        return 1.0
    else:
        return 0.0

def ConvectionVelocity(x, y, z):
    vel = KratosMultiphysics.Vector(3, 0.0)
    vel[0] = 1.0
    return vel

class BFECCConvectionTest(UnitTest.TestCase):
    def setUp(self):
        self.dt = 0.05
        self.end_time = 10.0
        self.mesh_divisions = 50
        self.bfecc_substeps = 10.0
        self.print_output = False
        self.check_tolerance = 1.0e-8
        self.print_reference_values = False
        self.work_folder = "."#"BFECCConvectionTest"

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
            mesh_parameters.AddEmptyValue("condition_name").SetString("LineCondition2D2N")
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
                    aux_temp = math.sqrt((x_local - x_c)**2 + (y_local - y_c)**2) - 0.1
                    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0, aux_temp)
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, [-y_local, x_local, 0.0])
                    # If the limiter is used, set some initial overshoots
                    # if self.has_limiter:
                    #     if aux_temp <= 0.0:
                    #         node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0, 1.0)
                    #     else:
                    #         node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0, 0.0)

            # Create the search structure
            locator = KratosMultiphysics.BinBasedFastPointLocator2D(self.main_model_part)
            locator.UpdateSearchDatabase()

            # Construct the BFECC utility
            if self.has_limiter:
                bfecc_utility = ConvectionDiffusionApplication.BFECCConvection2D(locator)
            else:
                bfecc_utility = ConvectionDiffusionApplication.BFECCLimiterConvection2D(locator)

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
                json_output_settings = KratosMultiphysics.Parameters(r'''{
                        "output_variables": ["TEMPERATURE"],
                        "output_file_name": "",
                        "model_part_name": "MainModelPart",
                        "time_frequency": 0.0
                }''')
                json_output_settings["output_file_name"].SetString(GetFilePath("BFECCConvectionTest/" + self.reference_file + "_results.json"))
                json_output_process = JsonOutputProcess(self.model, json_output_settings)
                json_output_process.ExecuteInitialize()
                json_output_process.ExecuteBeforeSolutionLoop()
                json_output_process.ExecuteFinalizeSolutionStep()
            else:
                json_check_parameters = KratosMultiphysics.Parameters(r'''{
                    "check_variables"      : ["TEMPERATURE"],
                    "input_file_name"      : "",
                    "model_part_name"      : "MainModelPart",
                    "time_frequency"       : 0.0
                }''')
                json_check_parameters["input_file_name"].SetString(GetFilePath("BFECCConvectionTest/" + self.reference_file + "_results.json"))
                json_check_process = FromJsonCheckResultProcess(self.model, json_check_parameters)
                json_check_process.ExecuteInitialize()
                json_check_process.ExecuteBeforeSolutionLoop()
                json_check_process.ExecuteFinalizeSolutionStep()

    def testBFECCConvection(self):
        self.has_limiter = False
        self.reference_file = "bfecc_convection_test"
        self.setUp()
        self.runTest()
        self.checkResults()
        self.tearDown()

    def testBFECCElementalLimiterConvection(self):
        self.has_limiter = True
        self.reference_file = "bfecc_elemental_limiter_convection_test"
        self.setUp()
        self.runTest()
        self.checkResults()
        self.tearDown()

    def testBFECCNodalLimiterConvection(self):
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        KratosMultiphysics.ModelPartIO(GetFilePath("BFECCConvectionTest/Semi1D_1x02_50")).ReadModelPart(model_part)
        model_part.SetBufferSize(2)

        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 0.05)

        # Create the search structure
        locator = KratosMultiphysics.BinBasedFastPointLocator2D(model_part)
        locator.UpdateSearchDatabase()

        bfecc_utility = ConvectionDiffusionApplication.BFECCConvection2D(
            locator,
            False,  # No partial dt (splitting)
            True)  # Activating nodal limiter

        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, BaseJumpedDistance(node.X,node.Y,node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, ConvectionVelocity(node.X,node.Y,node.Z))
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 1, ConvectionVelocity(node.X,node.Y,node.Z))

        for node in model_part.Nodes:
            if node.X < 0.001:
                node.Fix(KratosMultiphysics.DISTANCE)

        kratos_comm  = KratosMultiphysics.DataCommunicator.GetDefault()
        KratosMultiphysics.FindGlobalNodalNeighboursProcess(
                kratos_comm, model_part).Execute()

        KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(
            model_part,
            KratosMultiphysics.DISTANCE,
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.NODAL_AREA).Execute()

        for i in range(12):
            bfecc_utility.CopyScalarVarToPreviousTimeStep(
               model_part,
               KratosMultiphysics.DISTANCE)
            bfecc_utility.BFECCconvect(
                model_part,
                KratosMultiphysics.DISTANCE,
                KratosMultiphysics.VELOCITY,
                5) # substeps
            dt = model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
            model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, float(i)*dt)

        max_distance = -1.0
        min_distance = +1.0
        for node in model_part.Nodes:
            d =  node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            max_distance = max(max_distance, d)
            min_distance = min(min_distance, d)

        # gid_output = GiDOutputProcess(model_part,
        #                            "BFECCConvection_test_2D",
        #                            KratosMultiphysics.Parameters("""
        #                                {
        #                                    "result_file_configuration" : {
        #                                        "gidpost_flags": {
        #                                            "GiDPostMode": "GiD_PostBinary",
        #                                            "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                            "WriteConditionsFlag": "WriteConditions",
        #                                            "MultiFileFlag": "SingleFile"
        #                                        },
        #                                        "nodal_results"       : ["DISTANCE","VELOCITY"]
        #                                    }
        #                                }
        #                                """)
        #                            )

        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()

        self.assertAlmostEqual(max_distance, 1.0000000000000013)
        self.assertAlmostEqual(min_distance, 0.0)

if __name__ == '__main__':
    test = BFECCConvectionTest()
    test.testBFECCConvection()
    # test.testBFECCElementalLimiterConvection()
    test.testBFECCNodalLimiterConvection()
