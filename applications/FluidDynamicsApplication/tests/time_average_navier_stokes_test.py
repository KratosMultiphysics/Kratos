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

    def testCylinderFlow2DWater(self):
        self.work_folder = "TimeAveragedNavierStokesTest/cylinder"
        self.settings = "cylinder_2d_water.json"
        self.ExcecuteFlowTest()


    def testCylinderFlow2DAir(self):
        self.work_folder = "TimeAveragedNavierStokesTest/cylinder"
        self.settings = "cylinder_2d_air.json"
        self.ExcecuteFlowTest()


    def testBackStepFlow2D(self):
        self.work_folder = "TimeAveragedNavierStokesTest/back_step"
        self.settings = "back_step_flow_100.json"
        self.ExcecuteFlowTest()


    def testBackStepFlow2DReference(self):
        self.work_folder = "TimeAveragedNavierStokesTest/back_step"
        self.settings = "back_step_flow_reference_100.json"
        self.ExcecuteFlowTest()


    def testPipeFlow2D(self):
        self.work_folder = "TimeAveragedNavierStokesTest/pipe_flow"
        self.settings = "pipe_flow.json"
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


    def InitializeModelPart(self):
        current_model = Model()

        ##set the origin model part
        origin_model_part = current_model.CreateModelPart("fluid_computational_model_part")
        model_part_io = KratosMultiphysics.ModelPartIO("TimeAveragedNavierStokesTest/cylinder/cylinder_2d_water_initial_re1000")
        origin_model_part = model_part_io.ReadModelPart(origin_model_part)

        ##set the destination model part
        self.main_model_part = self.model.GetModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())

        #model_part_io = ModelPartIO(GetFilePath("TimeAveragedNavierStokesTest\cylinder\cylinder_2d_water_initial_re1000.mdpa"))

        ##copy the values to the destination model part
        VariableUtils().CopyModelPartNodalVar(VISCOSITY, origin_model_part, self.main_model_part, 0)
        VariableUtils().CopyModelPartNodalVar(DISPLACEMENT, origin_model_part, self.main_model_part, 0)

    
    def SetupMonolithicSolver(self):
        import json  
        with open(self.settings, 'r') as parameter_file:
            self.monolithic_json = json.load(parameter_file)
        self.monolithic_json["solver_settings"]["solver_type"] = "Monolithic"
        if "predictor_corrector" in self.monolithic_json["solver_settings"]:
            del self.monolithic_json["solver_settings"]["predictor_corrector"]
        if "time_averaging_acceleration" in self.monolithic_json["solver_settings"]:
            del self.monolithic_json["solver_settings"]["time_averaging_acceleration"]
        for process in self.monolithic_json["boundary_conditions_process_list"]:
            if process["Parameters"]["variable_name"] == "TIME_AVERAGED_VELOCITY":
                process["Parameters"]["variable_name"] = "VELOCITY"
            elif process["Parameters"]["variable_name"] == "TIME_AVERAGED_PRESSURE":
                process["Parameters"]["variable_name"] = "PRESSURE"
        self.monolithic_parameters = Parameters(json.dumps(self.monolithic_json))

        self.monolithic_model = Model()

        import python_solvers_wrapper_fluid
        self.monolithic_solver = python_solvers_wrapper_fluid.CreateSolver(self.monolithic_model, self.monolithic_parameters)
        self.monolithic_solver.AddVariables()
        self.monolithic_solver.ImportModelPart()
        self.monolithic_solver.PrepareModelPart()
        self.monolithic_solver.AddDofs()
        self.monolithic_solver.Initialize()

        self.monolithic_main_model_part = self.monolithic_model.GetModelPart(self.monolithic_parameters["problem_data"]["model_part_name"].GetString())
        
        import process_factory
        self.monolithic_list_of_processes  = process_factory.KratosProcessFactory(self.monolithic_model).ConstructListOfProcesses( self.monolithic_parameters["gravity"] )
        self.monolithic_list_of_processes += process_factory.KratosProcessFactory(self.monolithic_model).ConstructListOfProcesses( self.monolithic_parameters["boundary_conditions_process_list"] )
        for process in self.monolithic_list_of_processes:
            process.ExecuteInitialize()
            process.ExecuteBeforeSolutionLoop()

        monolithic_initialization_time = 10
        time = 0.0
        while(time <= monolithic_initialization_time):
            print ("Monolithic Initialization Time: ", str(time))
            for process in self.monolithic_list_of_processes:
                process.ExecuteInitializeSolutionStep()
            (self.monolithic_solver).InitializeSolutionStep()
            (self.monolithic_solver).Predict()
            (self.monolithic_solver).SolveSolutionStep()
            (self.monolithic_solver).FinalizeSolutionStep()
            for process in self.monolithic_list_of_processes:
                process.ExecuteFinalizeSolutionStep()
            time = self.monolithic_solver.AdvanceInTime(time)
        

    def InitializeMonolithicSolution(self):
        VariableUtils().CopyModelPartNodalVar(VELOCITY, TIME_AVERAGED_VELOCITY, self.monolithic_main_model_part, self.main_model_part, 0)
        VariableUtils().CopyModelPartNodalVar(PRESSURE, TIME_AVERAGED_PRESSURE, self.monolithic_main_model_part, self.main_model_part, 0)

        for i in range(0, self.solver.GetMinimumBufferSize()):
            self.main_model_part.CloneSolutionStep()
        

    def setUpProblem(self):
        with open(self.settings, 'r') as parameter_file:
            self.ProjectParameters = Parameters(parameter_file.read())
        
        self.SetupMonolithicSolver()

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
        if self.ProjectParameters.Has("auxiliar_process_list") == True:
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["auxiliar_process_list"] )

        ## Processes initialization
        for process in self.list_of_processes:
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

        self.InitializeMonolithicSolution()

        while(time <= end_time):
            print ("TIME: ", str(time))
            print ("STEP: ", str(step))

            for process in self.list_of_processes:
                process.ExecuteInitializeSolutionStep()

            self.gid_output.ExecuteInitializeSolutionStep()

            self.solver.InitializeSolutionStep()
            self.solver.Predict()
            self.solver.SolveSolutionStep()
            self.solver.FinalizeSolutionStep()

            for process in self.list_of_processes:
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
    # TimeAveragedNavierStokesTest().testCylinderFlow2DWater()
    # TimeAveragedNavierStokesTest().testCylinderFlow2DAir()
    # TimeAveragedNavierStokesTest().testBackStepFlow2D()
    # TimeAveragedNavierStokesTest().testPipeFlow2D()
