# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

# Import applications
import KratosMultiphysics.ShallowWaterApplication as SW

def CreateSolver(model, custom_settings):
    return ShallowWaterBaseSolver(model, custom_settings)

class ShallowWaterBaseSolver(PythonSolver):
    def __init__(self, model, settings):  # Constructor of the class
        super().__init__(model, settings)

        ## Set the element and condition names for the replace settings
        ## These should be defined in derived classes
        self.element_name = None
        self.condition_name = None
        self.min_buffer_size = 2

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        self.main_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, domain_size)

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(SW.HEIGHT)
        self.main_model_part.AddNodalSolutionStepVariable(KM.MOMENTUM)
        self.main_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(SW.FREE_SURFACE_ELEVATION)
        self.main_model_part.AddNodalSolutionStepVariable(KM.GRAVITY)
        self.main_model_part.AddNodalSolutionStepVariable(SW.BATHYMETRY)
        self.main_model_part.AddNodalSolutionStepVariable(SW.TOPOGRAPHY)
        self.main_model_part.AddNodalSolutionStepVariable(SW.MANNING)
        self.main_model_part.AddNodalSolutionStepVariable(SW.RAIN)
        self.main_model_part.AddNodalSolutionStepVariable(KM.NORMAL)

    def AddDofs(self):
        raise Exception("Calling the base class instead of the derived one")

    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        # Definition of the variables
        gravity = self.settings["gravity"].GetDouble()

        # Set ProcessInfo variables
        self.main_model_part.ProcessInfo.SetValue(KM.STEP, 0)
        self.main_model_part.ProcessInfo.SetValue(KM.GRAVITY_Z, gravity)

        if not self.main_model_part.ProcessInfo[KM.IS_RESTARTED]:
            ## Replace default elements and conditions
            self._ReplaceElementsAndConditions()
            ## Executes the check and prepare model process (Create computing_model_part)
            self._CheckAndPrepare()
            ## Set buffer size
            self.main_model_part.SetBufferSize(self.GetMinimumBufferSize())

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def GetComputingModelPart(self):
        return self.main_model_part

    def Initialize(self):
        self._GetSolutionStrategy().Check()
        self._GetSolutionStrategy().Initialize()
        KM.Logger.PrintInfo(self.__class__.__name__, "Initialization finished")

    def AdvanceInTime(self, current_time):
        current_time += self._GetEstimateDeltaTimeUtility().Execute()
        self.main_model_part.CloneTimeStep(current_time)
        self.main_model_part.ProcessInfo[KM.STEP] += 1
        return current_time

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self._GetSolutionStrategy().InitializeSolutionStep()

    def Predict(self):
        if self._TimeBufferIsInitialized():
            self._GetSolutionStrategy().Predict()

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            is_converged = self._GetSolutionStrategy().SolveSolutionStep()
            if not is_converged:
                KM.Logger.PrintInfo(self.__class__.__name__, "The solver did not converge")
            return is_converged
        else:
            return True

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self._GetSolutionStrategy().FinalizeSolutionStep()

    def Check(self):
        self._GetSolutionStrategy().Check()

    def Clear(self):
        self._GetSolutionStrategy().Clear()

    #### Specific internal functions ####

    def _TimeBufferIsInitialized(self):
        # We always have one extra old step (step 0, read from input)
        return self.main_model_part.ProcessInfo[KM.STEP] + 1 >= self.GetMinimumBufferSize()

    def _GetEstimateDeltaTimeUtility(self):
        if not hasattr(self, '_delta_time_utility'):
            self._delta_time_utility = self._CreateEstimateDeltaTimeUtility()
        return self._delta_time_utility

    def _CreateEstimateDeltaTimeUtility(self):
        # The c++ utility manages all the time step settings
        return SW.EstimateTimeStepUtility(self.GetComputingModelPart(), self.settings["time_stepping"])

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""
        {
            "solver_type"              : "shallow_water_base_solver",
            "model_part_name"          : "main_model_part",
            "domain_size"              : 2,
            "gravity"                  : 9.81,
            "model_import_settings"    : {
                "input_type"               : "mdpa",
                "input_filename"           : "unknown_name"
            },
            "echo_level"               : 0,
            "relative_tolerance"       : 1e-6,
            "absolute_tolerance"       : 1e-9,
            "maximum_iterations"       : 20,
            "compute_reactions"        : false,
            "reform_dofs_at_each_step" : false,
            "move_mesh_flag"           : false,
            "linear_solver_settings"   : {
                "solver_type"              : "amgcl"
            },
            "time_stepping"            : {
                "automatic_time_step"      : false,
                "time_step"                : 0.01
            }
        }""")
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def _ReplaceElementsAndConditions(self):
        ## Get number of nodes and domain size
        elem_num_nodes = self.__get_element_num_nodes()
        cond_num_nodes = self.__get_condition_num_nodes()
        domain_size = self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE]

        ## Complete the element name
        if (self.element_name is not None):
            new_elem_name = self.element_name + str(int(domain_size)) + "D" + str(int(elem_num_nodes)) + "N"
        else:
            raise Exception("There is no element name. Define the self.element_name string variable in your derived solver.")

        ## Complete the condition name
        if (self.condition_name is not None):
            new_cond_name = self.condition_name + str(int(domain_size)) + "D" + str(int(cond_num_nodes)) + "N"
        else:
            raise Exception("There is no condition name. Define the self.condition_name string variable in your derived solver.")

        ## Set the element and condition names in the Json parameters
        self.settings.AddValue("element_replace_settings", KM.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)

        ## Call the replace elements and conditions process
        KM.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

    def __get_element_num_nodes(self):
        if self.main_model_part.NumberOfElements() != 0:
            element_num_nodes = len(self.main_model_part.Elements.__iter__().__next__().GetNodes())
        else:
            element_num_nodes = 0
        element_num_nodes = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(element_num_nodes)
        return element_num_nodes

    def __get_condition_num_nodes(self):
        if self.main_model_part.NumberOfConditions() != 0:
            condition_num_nodes = len(self.main_model_part.Conditions.__iter__().__next__().GetNodes())
        else:
            condition_num_nodes = 2
        condition_num_nodes = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(condition_num_nodes)
        return condition_num_nodes

    def _CheckAndPrepare(self):
        pass

    def _GetLinearSolver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._CreateLinearSolver()
        return self._linear_solver

    def _GetBuilderAndSolver(self):
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._CreateBuilderAndSolver()
        return self._builder_and_solver

    def _GetConvergenceCriterion(self):
        if not hasattr(self, '_convergence_criterion'):
            self._convergence_criterion = self._CreateConvergenceCriterion()
        return self._convergence_criterion

    def _GetScheme(self):
        if not hasattr(self, '_scheme'):
            self._scheme = self._CreateScheme()
        return self._scheme

    def _GetSolutionStrategy(self):
        if not hasattr(self, '_solution_strategy'):
            self._solution_strategy = self._CreateSolutionStrategy()
        return self._solution_strategy

    def _CreateLinearSolver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        return linear_solver_factory.ConstructSolver(linear_solver_configuration)

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        builder_and_solver = KM.ResidualBasedBlockBuilderAndSolver(linear_solver)
        return builder_and_solver

    def _CreateConvergenceCriterion(self):
        convergence_criterion = KM.DisplacementCriteria(
            self.settings["relative_tolerance"].GetDouble(),
            self.settings["absolute_tolerance"].GetDouble())
        convergence_criterion.SetEchoLevel(self.echo_level)
        return convergence_criterion

    def _CreateScheme(self):
        time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        return time_scheme

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        scheme = self._GetScheme()
        convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        strategy = KM.ResidualBasedNewtonRaphsonStrategy(
            computing_model_part,
            scheme,
            convergence_criterion,
            builder_and_solver,
            self.settings["maximum_iterations"].GetInt(),
            self.settings["compute_reactions"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool(),
            self.settings["move_mesh_flag"].GetBool())
        strategy.SetEchoLevel(max(0, self.echo_level-1))
        return strategy
