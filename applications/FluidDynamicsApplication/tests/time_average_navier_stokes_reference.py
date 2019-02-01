import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True
except ImportError as e:
    have_external_solvers = False

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class TimeAveragedNavierStokesTest(UnitTest.TestCase):

    def testCylinderFlow2DReferenceWater(self):
        self.work_folder = "TimeAveragedNavierStokesTest/cylinder"
        self.settings = "cylinder_2d_reference_water.json"
        self.ExcecuteFlowTest()


    def testCylinderFlow2DReferenceAir(self):
        self.work_folder = "TimeAveragedNavierStokesTest/cylinder"
        self.settings = "cylinder_2d_reference_air.json"
        self.ExcecuteFlowTest()


    def testBackStepFlow2DReference(self):
        self.work_folder = "TimeAveragedNavierStokesTest/back_step"
        self.settings = "back_step_flow_reference_100.json"
        self.ExcecuteFlowTest()


    def testPipeFlow2DReference(self):
        self.work_folder = "TimeAveragedNavierStokesTest/pipe_flow"
        self.settings = "pipe_flow_reference.json"
        self.ExcecuteFlowTest()


    def ExcecuteFlowTest(self):
        with WorkFolderScope(self.work_folder):
            print(os.path.realpath(__file__))
            self.setUp()
            self.setUpProblem()
            self.runTest()
            self.tearDown()

            if self.print_reference_values == True:
                self.checkResults()


    def setUp(self):
        self.check_tolerance = 1e-6
        self.print_output = True
        self.print_reference_values = False


    def tearDown(self):
        try:
            os.remove(self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()+'.time')
        except FileNotFoundError as e:
            pass

    def setUpProblem(self):
        with open(self.settings, 'r') as parameter_file:
            self.ProjectParameters = Parameters(parameter_file.read())

        self.model = Model()
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

        ## Defining main model part
        self.main_model_part = self.model.GetModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        #self.InitializeModelPart()

        ## Processes construction
        import process_factory
        self.list_of_processes  = process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["gravity"] )
        self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["boundary_conditions_process_list"] )
        self.appy_time_averaging_process = process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["auxiliar_process_list"] )

        ## Processes initialization
        for process in self.list_of_processes:
            process.ExecuteInitialize()
        for process in self.appy_time_averaging_process:
            process.ExecuteInitialize()

        ## Set results file Configuration
        if (self.print_output):
            self.gid_output = self.setUpGiDOutput()
            self.gid_output.ExecuteInitialize()
            self.gid_output.ExecuteBeforeSolutionLoop()


    def runTest(self):

        end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

        time = 0.0
        step = 0

        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        while(time <= end_time):
            print ("TIME: ", str(time))
            print ("STEP: ", str(step))

            for process in self.list_of_processes:
                process.ExecuteInitializeSolutionStep()
            if (step >= 10):
                for process in self.appy_time_averaging_process:
                    process.ExecuteInitializeSolutionStep()

            self.gid_output.ExecuteInitializeSolutionStep()

            self.solver.InitializeSolutionStep()
            self.solver.Predict()
            self.solver.SolveSolutionStep()
            self.solver.FinalizeSolutionStep()

            for process in self.list_of_processes:
                process.ExecuteFinalizeSolutionStep()
            
            if (step >= 10):
                for process in self.appy_time_averaging_process:
                    process.ExecuteFinalizeSolutionStep()

            if (self.print_output):
                if self.gid_output.IsOutputStep():
                    self.gid_output.PrintOutput()
            
            time = self.solver.AdvanceInTime(time)
            step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

        for process in self.list_of_processes:
            process.ExecuteFinalize()

        if (self.print_output):
            self.gid_output.ExecuteFinalize()
        
        #self.solver.ExportModelPart()


    def setUpGiDOutput(self):
        self.parallel_type = self.ProjectParameters["problem_data"]["parallel_type"].GetString()
        if self.parallel_type == "OpenMP":
            from gid_output_process import GiDOutputProcess as OutputProcess
        elif self.parallel_type == "MPI":
            from gid_output_process_mpi import GiDOutputProcessMPI as OutputProcess

        gid_output = OutputProcess( self.main_model_part,
                                self.ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                self.ProjectParameters["output_configuration"])

        return gid_output

if __name__ == '__main__':
    TimeAveragedNavierStokesTest().testCylinderFlow2DReferenceWater()
    # TimeAveragedNavierStokesTest().testCylinderFlow2DReferenceAir()
    # TimeAveragedNavierStokesTest().testBackStepFlow2DReference()
    # TimeAveragedNavierStokesTest().testPipeFlow2DReference()
