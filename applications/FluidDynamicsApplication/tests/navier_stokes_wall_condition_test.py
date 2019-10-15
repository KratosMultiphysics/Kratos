import KratosMultiphysics
import KratosMultiphysics.process_factory as process_factory
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid

import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

have_external_solvers = KratosUtilities.CheckIfApplicationsAvailable("ExternalSolversApplication")

@UnitTest.skipUnless(have_external_solvers,"Missing required application: ExternalSolversApplication")
class NavierStokesWallConditionTest(UnitTest.TestCase):
    def testNavierStokesWallCondition(self):
        self.setUp()
        self.setUpProblem()
        self.runTest()
        self.tearDown()
        self.checkResults()

    def setUp(self):
        self.ext_pres_res = 1000.0
        self.check_tolerance = 2e-2
        self.print_output = False
        self.work_folder = "NavierStokesWallConditionTest"
        self.settings = "NavierStokesWallConditionTestParameters.json"

    def tearDown(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting(
                self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()+'.time')

    def setUpProblem(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            with open(self.settings, 'r') as parameter_file:
                self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            self.model = KratosMultiphysics.Model()

            ## Solver construction
            self.solver = python_solvers_wrapper_fluid.CreateSolver(self.model, self.ProjectParameters)

            self.solver.AddVariables()

            ## Read the model - note that SetBufferSize is done here
            self.solver.ImportModelPart()
            self.solver.PrepareModelPart()

            ## Add AddDofs
            self.solver.AddDofs()

            ## Solver initialization
            self.solver.Initialize()

            ## Processes construction
            self.list_of_processes  = process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["processes"]["gravity"] )
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["processes"]["boundary_conditions_process_list"] )

            ## Processes initialization
            for process in self.list_of_processes:
                process.ExecuteInitialize()

            self.main_model_part = self.model.GetModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())

    def runTest(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            if (self.print_output):
                gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
                multifile = KratosMultiphysics.MultiFileFlag.SingleFile
                deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
                write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly
                gid_io = KratosMultiphysics.GidIO(self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString(),gid_mode,multifile,deformed_mesh_flag, write_conditions)

                mesh_name = 0.0
                gid_io.InitializeMesh( mesh_name)
                gid_io.WriteMesh( self.main_model_part.GetMesh() )
                gid_io.FinalizeMesh()
                gid_io.InitializeResults(mesh_name,(self.main_model_part).GetMesh())

            end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

            time = 0.0
            step = 0

            for process in self.list_of_processes:
                process.ExecuteBeforeSolutionLoop()

            while(time <= end_time):

                time = self.solver.AdvanceInTime(time)

                for process in self.list_of_processes:
                    process.ExecuteInitializeSolutionStep()

                self.solver.InitializeSolutionStep()
                self.solver.Predict()
                self.solver.SolveSolutionStep()
                self.solver.FinalizeSolutionStep()

                for process in self.list_of_processes:
                    process.ExecuteFinalizeSolutionStep()

                for process in self.list_of_processes:
                    process.ExecuteBeforeOutputStep()

                if (self.print_output):
                    gid_io.WriteNodalResults(KratosMultiphysics.VELOCITY,self.main_model_part.Nodes,time,0)
                    gid_io.WriteNodalResults(KratosMultiphysics.PRESSURE,self.main_model_part.Nodes,time,0)
                    gid_io.WriteNodalResults(KratosMultiphysics.DISTANCE,self.main_model_part.Nodes,time,0)

                for process in self.list_of_processes:
                    process.ExecuteAfterOutputStep()

            for process in self.list_of_processes:
                process.ExecuteFinalize()

            if (self.print_output):
                gid_io.FinalizeResults()

    def checkResults(self):
        react_x = 0.0

        for node in self.main_model_part.Nodes:
            react_x += node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)

        self.assertAlmostEqual(-react_x, self.ext_pres_res, delta=2e-2)

if __name__ == '__main__':
    test = NavierStokesWallConditionTest()
    test.setUp()
    test.setUpProblem()
    test.runTest()
    test.tearDown()
    test.checkResults()
