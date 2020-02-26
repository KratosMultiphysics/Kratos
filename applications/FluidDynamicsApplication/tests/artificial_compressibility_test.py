import KratosMultiphysics
import KratosMultiphysics.process_factory as process_factory
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics.KratosUnittest as UnitTest

from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid

class ArtificialCompressibilityTest(UnitTest.TestCase):
    def testArtificialCompressibility(self):
        self.setUp()
        self.setUpProblem()
        self.runTest()
        self.tearDown()
        self.checkResults()

    def setUp(self):
        self.check_tolerance = 1e-6
        self.print_output = True
        self.print_reference_values = False
        self.work_folder = "ArtificialCompressibilityTest"
        self.reference_file = "reference_cavity_compressibility"
        self.settings = "ArtificialCompressibilityTestParameters.json"

    def tearDown(self):
        with UnitTest.WorkFolderScope(self.work_folder,__file__):
            KratosUtilities.DeleteFileIfExisting(
                self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()+'.time')

    def setUpProblem(self):
        with UnitTest.WorkFolderScope(self.work_folder,__file__):
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

            self.main_model_part = self.model.GetModelPart(self.ProjectParameters["solver_settings"]["model_part_name"].GetString())

    def runTest(self):
        with UnitTest.WorkFolderScope(self.work_folder,__file__):
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
        with UnitTest.WorkFolderScope(self.work_folder,__file__):
            if self.print_reference_values:
                with open(self.reference_file+'.csv','w') as ref_file:
                    ref_file.write("#ID, VELOCITY_X, VELOCITY_Y\n")
                    for node in self.main_model_part.Nodes:
                        vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                        ref_file.write("{0}, {1}, {2}\n".format(node.Id, vel[0], vel[1]))
            else:
                with open(self.reference_file+'.csv','r') as reference_file:
                    reference_file.readline() # skip header
                    line = reference_file.readline()

                    for node in self.main_model_part.Nodes:
                        values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                        node_id = values[0]
                        reference_vel_x = values[1]
                        reference_vel_y = values[2]

                        velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                        self.assertAlmostEqual(reference_vel_x, velocity[0], delta = self.check_tolerance)
                        self.assertAlmostEqual(reference_vel_y, velocity[1], delta = self.check_tolerance)

                        line = reference_file.readline()
                    if line != '': # If we did not reach the end of the reference file
                        self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

if __name__ == '__main__':
    test = ArtificialCompressibilityTest()
    test.setUp()
    test.setUpProblem()
    test.runTest()
    test.tearDown()
    test.checkResults()
