from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

class BaseKratosAnalysisStage(object):
    """The base class for the analysis classes in the applications
    """
    def __init__(self, project_parameters, external_model_part=None):
        """The constructor of the Analysis-Object.
        It obtains the project parameters used for the analysis
        This function is intended to be called from the constructor
        of deriving classes:
        super(DerivedAnalysis, self).__init__(project_parameters)

        It is intended that this constructor creates the solver and
        adds the Variables that are needed by the solver to the
        ModelPart

        An external ModelPart is used if the AnalysisStage is used
        in a larger context, e.g. in a MultiStage-Analysis or Optimization
        """
        if (type(project_parameters) == str): # a file name is provided
            with open(project_parameters,'r') as parameter_file:
                self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())
        elif (type(project_parameters) == KratosMultiphysics.Parameters): # a Parameters object is provided
            self.ProjectParameters = project_parameters
        else:
            raise Exception("Input is expected to be provided as a Kratos Parameters object or a file name")

        if external_model_part is not None:
            if (type(external_model_part) != KratosMultiphysics.ModelPart):
                raise Exception("Input is expected to be provided as a Kratos ModelPart object")
            self.using_external_model_part = True

        ## Get echo level and parallel type
        self.echo_level = self.ProjectParameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.ProjectParameters["problem_data"]["parallel_type"].GetString()

        ## Import parallel modules if needed
        if (self.parallel_type == "MPI"):
            import KratosMultiphysics.mpi as KratosMPI
            import KratosMultiphysics.MetisApplication as MetisApplication
            import KratosMultiphysics.TrilinosApplication as TrilinosApplication
            self.is_printing_rank = (KratosMPI.mpi.rank == 0)
        else:
            self.is_printing_rank = True

        ## model part definition
        if self.using_external_model_part:
            self.main_model_part = external_model_part
        else:
            main_model_part_name = self.ProjectParameters["problem_data"]["model_part_name"].GetString()
            self.main_model_part = KratosMultiphysics.ModelPart(main_model_part_name)
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE,
                                                      self.ProjectParameters["problem_data"]["domain_size"].GetInt())

        self._CreateSolver()

        self.solver.AddVariables()

        if not self.using_external_model_part:
            ## Read the model
            self._ReadModelPart()


    #### Public functions to run the Analysis ####
    def Run(self):
        """This function executes the entire analysis
        It is NOT intended to be overridden in deriving classes!
        """
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

    def RunMainTemporalLoop(self):
        """This function executes the temporal loop of the analysis
        It is NOT intended to be overridden in deriving classes!
        """
        while self.time < self.end_time:
            self.InitializeTimeStep()
            self.SolveTimeStep()
            self.FinalizeTimeStep()

    def Initialize(self):
        """This function initializes the analysis
        Usage: It is designed to be called ONCE, BEFORE the execution of the time-loop
        This function IS intended to be overridden in deriving classes!
        At the end of this function the StageAnalysis is ready for the time-loop
        It should be called AFTER the ModelPart used for this AnalysisStage is used
        """
        pass

    def Finalize(self):
        """This function finalizes the analysis
        Usage: It is designed to be called ONCE, AFTER the execution of the time-loop
        This function IS intended to be overridden in deriving classes!
        """
        pass

    def InitializeTimeStep(self):
        """This function initializes the time-step
        Usage: It is designed to be called once at the beginning of EACH time-step
        This function IS intended to be overridden in deriving classes!
        """
        pass

    def SolveStep(self):
        """This function solves one step
        It can be called several times during one time-step
        This is equivalent to calling "solving_strategy.Solve()" (without "Initialize")
        This function is NOT intended to be overridden in deriving classes!
        """
        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()

    def FinalizeTimeStep(self):
        """This function finalizes the time-step
        Usage: It is designed to be called once at the end of EACH time-step
        This function IS intended to be overridden in deriving classes!
        """
        pass

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be done
        (for each step) before solving the solution step.
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived\
            analysis classes")

    def Predict(self):
        """This function predicts the solution
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived\
            analysis classes")

    def SolveSolutionStep(self):
        """This function solves the current step
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived\
            analysis classes")

    def FinalizeSolutionStep(self):
        """This function Performs all the required operations that should be done
        (for each step) after solving the solution step.
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived\
            analysis classes")

    def _ReadModelPart(self):
        """This function reads the ModelPart, in case it is not provided to the AnalysisStage
        This function is NOT intended to be overridden in deriving classes!
        """
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            # Import model part from mdpa file.
            KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]::", "Reading model part from file: " + os.path.join(problem_path, input_filename) + ".mdpa")
            KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)
            KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]::", "Finished reading model part from mdpa file.")
        elif
            ...
        else:
            Error
        # import materials if applicable
        if self.ProjectParameters["material_import_settings"].Has("materials_filename"):
            materials_filename = self.ProjectParameters["material_import_settings"]["materials_filename"].GetString()
            if (materials_filename != ""):
                import read_materials_process
                # Create a dictionary of model parts.
                Model = KratosMultiphysics.Model()
                Model.AddModelPart(self.main_model_part)
                # Add constitutive laws and material properties from json file to model parts.
                read_materials_process.ReadMaterialsProcess(Model, self.settings["material_import_settings"])

    def _CreateSolver(self):
        """This function creates the solver of the AnalysisStage
        (typically by importing it using the python solver wrappers)
        This function has to be implemented in deriving classes!
        """
        raise NotImplementedError("This function has to be implemented by derived\
            analysis classes")