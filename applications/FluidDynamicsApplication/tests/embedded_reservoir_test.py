import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True
except ImportError as e:
    have_external_solvers = False


import KratosMultiphysics.KratosUnittest as UnitTest

import os

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

@UnitTest.skipUnless(have_external_solvers,"Missing required application: ExternalSolversApplication")
class EmbeddedReservoirTest(UnitTest.TestCase):
    def testEmbeddedReservoir2D(self):
        self.distance = 0.5
        self.slip_level_set = False
        self.work_folder = "EmbeddedReservoirTest"
        self.reference_file = "reference_reservoir_2D"
        self.settings = "EmbeddedReservoir2DTest_parameters.json"
        self.ExecuteEmbeddedReservoirTest()

    def testEmbeddedReservoir3D(self):
        self.distance = 0.5
        self.slip_level_set = False
        self.work_folder = "EmbeddedReservoirTest"
        self.reference_file = "reference_reservoir_3D"
        self.settings = "EmbeddedReservoir3DTest_parameters.json"
        self.ExecuteEmbeddedReservoirTest()

    def testEmbeddedSlipReservoir2D(self):
        self.distance = 0.5
        self.slip_level_set = True
        self.work_folder = "EmbeddedReservoirTest"
        self.reference_file = "reference_slip_reservoir_2D"
        self.settings = "EmbeddedReservoir2DTest_parameters.json"
        self.ExecuteEmbeddedReservoirTest()

    def testEmbeddedSlipReservoir3D(self):
        self.distance = 0.5
        self.slip_level_set = True
        self.work_folder = "EmbeddedReservoirTest"
        self.reference_file = "reference_slip_reservoir_3D"
        self.settings = "EmbeddedReservoir3DTest_parameters.json"
        self.ExecuteEmbeddedReservoirTest()

    def ExecuteEmbeddedReservoirTest(self):
        with WorkFolderScope(self.work_folder):
            self.setUp()
            self.setUpProblem()
            self.setUpDistanceField()
            self.runTest()
            self.tearDown()
            self.checkResults()

    def setUp(self):
        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False

    def tearDown(self):
        with WorkFolderScope(self.work_folder):
            try:
                os.remove(self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()+'.time')
            except FileNotFoundError as e:
                pass

    def setUpProblem(self):
        with WorkFolderScope(self.work_folder):
            with open(self.settings, 'r') as parameter_file:
                self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            self.model = KratosMultiphysics.Model()

            ## Solver construction
            import python_solvers_wrapper_fluid
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
            import process_factory
            self.list_of_processes  = process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["gravity"] )
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["boundary_conditions_process_list"] )

            ## Processes initialization
            for process in self.list_of_processes:
                process.ExecuteInitialize()

            self.main_model_part = self.model.GetModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())

    def setUpDistanceField(self):
        # Set the distance function
        if (self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            for node in self.main_model_part.Nodes:
                distance = node.Y-self.distance
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, distance)
        elif (self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
            for node in self.main_model_part.Nodes:
                distance = node.Z-self.distance
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, distance)

        # Set the ELEMENTAL_DISTANCES value
        n_nodes = len(self.main_model_part.Elements[1].GetNodes())
        for element in self.main_model_part.Elements:
            elem_dist = KratosMultiphysics.Vector(n_nodes)
            elem_nodes = element.GetNodes()
            for i_node in range(0,n_nodes):
                elem_dist[i_node] = elem_nodes[i_node].GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            element.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, elem_dist)

        # If proceeds, set the SLIP flag
        if (self.slip_level_set):
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLIP, True, self.main_model_part.Elements)

    def runTest(self):
        with WorkFolderScope(self.work_folder):
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
        with WorkFolderScope(self.work_folder):
            if self.print_reference_values:
                with open(self.reference_file+'.csv','w') as ref_file:
                    ref_file.write("#ID, PRESSURE\n")
                    for node in self.main_model_part.Nodes:
                        pres = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
                        ref_file.write("{0}, {1}\n".format(node.Id, pres))
            else:
                with open(self.reference_file+'.csv','r') as reference_file:
                    reference_file.readline() # skip header
                    line = reference_file.readline()

                    for node in self.main_model_part.Nodes:
                        values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                        node_id = values[0]
                        reference_pres = values[1]

                        pres = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
                        self.assertAlmostEqual(reference_pres, pres, delta = self.check_tolerance)

                        line = reference_file.readline()
                    if line != '': # If we did not reach the end of the reference file
                        self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

if __name__ == '__main__':
    test = EmbeddedReservoirTest()
    test.setUp()
    test.distance = 0.5
    test.slip_level_set = False
    test.print_output = False
    test.print_reference_values = False
    test.work_folder = "EmbeddedReservoirTest"
    test.reference_file = "reference_slip_reservoir_2D"
    test.settings = "EmbeddedReservoir2DTest_parameters.json"
    test.setUpProblem()
    test.setUpDistanceField()
    test.runTest()
    test.tearDown()
    test.checkResults()

