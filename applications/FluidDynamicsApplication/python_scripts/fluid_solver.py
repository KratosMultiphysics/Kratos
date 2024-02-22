# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics import auxiliary_solver_utilities
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.FluidDynamicsApplication import check_and_prepare_model_process_fluid

def CreateSolver(model, custom_settings):
    return FluidSolver(model, custom_settings)

class FluidSolver(PythonSolver):
    """The base class for fluid dynamics solvers.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Depending on the formulation type, derived classes may require to
    override some (or all) the following functions:

    _CreateScheme
    _CreateConvergenceCriterion
    _CreateLinearSolver
    _CreateBuilderAndSolver
    _CreateSolutionStrategy

    The solution strategy, builder_and_solver, etc. should alway be retrieved
    using the getter functions _GetSolutionStrategy, _GetBuilderAndSolver,
    etc. from this base class.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the model containing the modelpart used to construct the solver.
    settings -- Kratos parameters containing solver settings.
    """
    def __init__(self, model, settings):

        super(FluidSolver,self).__init__(model, settings)

        ## Set the element and condition names for the replace settings
        ## These should be defined in derived classes
        self.element_name = None
        self.condition_name = None
        self.min_buffer_size = 3
        self._enforce_element_and_conditions_replacement = False #TODO: Remove once we remove the I/O from the solver

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size == -1:
            raise Exception('Please provide the domain size as the "domain_size" (int) parameter!')

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

    def AddVariables(self):
        raise Exception("Trying to call FluidSolver.AddVariables(). Implement the AddVariables() method in the specific derived solver.")

    def AddDofs(self):
        dofs_and_reactions_to_add = []
        dofs_and_reactions_to_add.append(["VELOCITY_X", "REACTION_X"])
        dofs_and_reactions_to_add.append(["VELOCITY_Y", "REACTION_Y"])
        dofs_and_reactions_to_add.append(["VELOCITY_Z", "REACTION_Z"])
        dofs_and_reactions_to_add.append(["PRESSURE", "REACTION_WATER_PRESSURE"])
        KratosMultiphysics.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Fluid solver DOFs added correctly.")

    def GetDofsList(self):
        """This function creates and returns a list with the DOFs defined in the conditions and elements specifications
        Note that this requires the main_model_part to be already set, that is to say to have already performed the element substitution (see PrepareModelPart).
        """
        return KratosMultiphysics.SpecificationsUtilities.GetDofsListFromSpecifications(self.main_model_part)

    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if not self.is_restarted():
            ## Set fluid properties from materials json file
            materials_imported = self._SetPhysicalProperties()
            if not materials_imported:
                KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, "Material properties have not been imported. Check \'material_import_settings\' in your ProjectParameters.json.")
            ## Replace default elements and conditions
            use_input_model_part = self.settings["model_import_settings"]["input_type"].GetString() == "use_input_model_part"
            if not (use_input_model_part and self._enforce_element_and_conditions_replacement):
                self._ReplaceElementsAndConditions()
            ## Set and fill buffer
            self._SetAndFillBuffer()

        ## Executes the check and prepare model process. Always executed as it also assigns neighbors which are not saved in a restart
        self._ExecuteCheckAndPrepare()

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Model reading finished.")

    def ExportModelPart(self):
        ## Model part writing
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Model export finished.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def Initialize(self):
        raise Exception("Calling FluidSolver.Initialize() base method. Please implement a custom Initialize() method for your solver.")

    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        self._GetSolutionStrategy().InitializeSolutionStep()

    def Predict(self):
        self._GetSolutionStrategy().Predict()

    def SolveSolutionStep(self):
        is_converged = self._GetSolutionStrategy().SolveSolutionStep()
        if not is_converged:
            msg  = "Fluid solver did not converge for step " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]) + "\n"
            msg += "corresponding to time " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]) + "\n"
            KratosMultiphysics.Logger.PrintWarning(self.__class__.__name__, msg)
        return is_converged

    def FinalizeSolutionStep(self):
        self._GetSolutionStrategy().FinalizeSolutionStep()

    def Check(self):
        self._GetSolutionStrategy().Check()

    def Clear(self):
        self._GetSolutionStrategy().Clear()

    def GetComputingModelPart(self):
        if not self.main_model_part.HasSubModelPart("fluid_computational_model_part"):
            raise Exception("The ComputingModelPart was not created yet!")
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")

    ## FluidSolver specific methods.
    def _ReplaceElementsAndConditions(self):
        ## Get number of nodes and domain size
        elem_num_nodes = self._GetElementNumNodes()
        cond_num_nodes = self._GetConditionNumNodes()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ## If there are no elements and/or conditions, default to triangles/tetra meshes to avoid breaking the ReplaceElementsAndConditionsProcess
        ## This only affects the input name (if there are no elements or conditions to replace, nothing is replaced).
        if elem_num_nodes == 0:
            elem_num_nodes = domain_size + 1
        if cond_num_nodes == 0:
            cond_num_nodes = domain_size

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
        #self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""{}""")
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)

        ## Call the replace elements and conditions process
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

    def _GetElementNumNodes(self):
        if self.main_model_part.NumberOfElements() != 0:
            element_num_nodes = len(self.main_model_part.Elements.__iter__().__next__().GetNodes())
        else:
            element_num_nodes = 0

        element_num_nodes = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(element_num_nodes)
        return element_num_nodes

    def _GetConditionNumNodes(self):
        if self.main_model_part.NumberOfConditions() != 0:
            condition_num_nodes = len(self.main_model_part.Conditions.__iter__().__next__().GetNodes())
        else:
            condition_num_nodes = 0

        condition_num_nodes = self.main_model_part.GetCommunicator().GetDataCommunicator().MaxAll(condition_num_nodes)
        return condition_num_nodes

    def _ExecuteCheckAndPrepare(self):
        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
        prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])
        prepare_model_part_settings.AddValue("assign_neighbour_elements_to_conditions",self.settings["assign_neighbour_elements_to_conditions"])

        check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()

    def _SetAndFillBuffer(self):
        init_dt = self._ComputeInitialDeltaTime()
        required_buffer_size = self.GetMinimumBufferSize()
        auxiliary_solver_utilities.SetAndFillBuffer(self.main_model_part, required_buffer_size, init_dt)

    def _ComputeDeltaTime(self):
        # Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            delta_time = self.GetEstimateDtUtility().EstimateDt()
        # User-defined delta time
        else:
            delta_time = self.settings["time_stepping"]["time_step"].GetDouble()

        return delta_time

    def _ComputeInitialDeltaTime(self):
        # Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            initial_delta_time = self.settings["time_stepping"]["minimum_delta_time"].GetDouble()
        # User-defined delta time
        else:
            initial_delta_time = self.settings["time_stepping"]["time_step"].GetDouble()

        return initial_delta_time

    def _SetPhysicalProperties(self):
        # Check if the fluid properties are provided using a .json file
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            materials_imported = True
        else:
            materials_imported = False

        # If the element uses nodal material properties, transfer them to the nodes
        if self.element_has_nodal_properties:
            self._SetNodalProperties()

        return materials_imported

    def _SetNodalProperties(self):
        err_msg = "Calling base FluidSolver \'_SetNodalProperties\' method.\n"
        err_msg += "This must be implemented in the derived solver in accordance with the element formulation."
        raise Exception(err_msg)

    # TODO: I THINK THIS SHOULD BE MOVED TO THE BASE PYTHON SOLVER
    def is_restarted(self):
        # this function avoids the long call to ProcessInfo and is also safer
        # in case the detection of a restart is changed later
        return self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]

    def GetEstimateDtUtility(self):
        if not hasattr(self, '_estimate_dt_utility'):
            self._estimate_dt_utility = self._CreateEstimateDtUtility()
        return self._estimate_dt_utility

    def _GetScheme(self):
        if not hasattr(self, '_scheme'):
            self._scheme = self._CreateScheme()
        return self._scheme

    def _GetConvergenceCriterion(self):
        if not hasattr(self, '_convergence_criterion'):
            self._convergence_criterion = self._CreateConvergenceCriterion()
        return self._convergence_criterion

    def _GetLinearSolver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._CreateLinearSolver()
        return self._linear_solver

    def _GetBuilderAndSolver(self):
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._CreateBuilderAndSolver()
        return self._builder_and_solver

    def _GetSolutionStrategy(self):
        if not hasattr(self, '_solution_strategy'):
            self._solution_strategy = self._CreateSolutionStrategy()
        return self._solution_strategy

    def _CreateEstimateDtUtility(self):
        estimate_dt_utility = KratosCFD.EstimateDtUtility(
                self.GetComputingModelPart(),
                self.settings["time_stepping"])

        return estimate_dt_utility

    def _CreateScheme(self):
        domain_size = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # Cases in which the element manages the time integration
        if self.element_integrates_in_time:
            # "Fake" scheme for those cases in where the element manages the time integration
            # It is required to perform the nodal update once the current time step is solved
            scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(
                domain_size,
                domain_size + 1)
            # In case the BDF2 scheme is used inside the element, the BDF time discretization utility is required to update the BDF coefficients
            if (self.settings["time_scheme"].GetString() == "bdf2"):
                time_order = 2
                self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)
            else:
                if  (self.settings["time_scheme"].GetString()!= "crank_nicolson"):
                    err_msg = "Requested elemental time scheme \"" + self.settings["time_scheme"].GetString()+ "\" is not available.\n"
                    err_msg += "Available options are: \"bdf2\" and \"crank_nicolson\""
                    raise Exception(err_msg)
        # Cases in which a time scheme manages the time integration
        else:
            # Bossak time integration scheme
            if self.settings["time_scheme"].GetString() == "bossak":
                if self.settings["consider_periodic_conditions"].GetBool() == True:
                    scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                        self.settings["alpha"].GetDouble(),
                        domain_size,
                        KratosCFD.PATCH_INDEX)
                else:
                    scheme = KratosCFD.ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
                        self.settings["alpha"].GetDouble(),
                        self.settings["move_mesh_strategy"].GetInt(),
                        domain_size)
            # BDF2 time integration scheme
            elif self.settings["time_scheme"].GetString() == "bdf2":
                scheme = KratosCFD.BDF2TurbulentScheme()
            # Time scheme for steady state fluid solver
            elif self.settings["time_scheme"].GetString() == "steady":
                scheme = KratosCFD.ResidualBasedSimpleSteadyScheme(
                        self.settings["velocity_relaxation"].GetDouble(),
                        self.settings["pressure_relaxation"].GetDouble(),
                        domain_size)
            else:
                if  (self.settings["time_scheme"].GetString()!= "crank_nicolson"):
                    err_msg = "Requested time scheme " + self.settings["time_scheme"].GetString() + " is not available.\n"
                    err_msg += "Available options are: \"bossak\", \"bdf2\" ,\"steady\" and \"crank_nicolson\""
                    raise Exception(err_msg)

        return scheme

    def _CreateLinearSolver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        return linear_solver_factory.ConstructSolver(linear_solver_configuration)

    def _CreateConvergenceCriterion(self):
        if self.settings["time_scheme"].GetString() == "steady":
            convergence_criterion = KratosMultiphysics.ResidualCriteria(
                self.settings["relative_velocity_tolerance"].GetDouble(),
                self.settings["absolute_velocity_tolerance"].GetDouble())
        else:
            convergence_criterion = KratosMultiphysics.MixedGenericCriteria(
                [(KratosMultiphysics.VELOCITY, self.settings["relative_velocity_tolerance"].GetDouble(), self.settings["absolute_velocity_tolerance"].GetDouble()),
                (KratosMultiphysics.PRESSURE, self.settings["relative_pressure_tolerance"].GetDouble(), self.settings["absolute_pressure_tolerance"].GetDouble())])
        convergence_criterion.SetEchoLevel(self.settings["echo_level"].GetInt())
        return convergence_criterion

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        if self.settings["consider_periodic_conditions"].GetBool():
            builder_and_solver = KratosCFD.ResidualBasedBlockBuilderAndSolverPeriodic(
                linear_solver,
                KratosCFD.PATCH_INDEX)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        return builder_and_solver

    def _CreateSolutionStrategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            solution_strategy = self._CreateLinearStrategy()
        elif analysis_type == "non_linear":
            solution_strategy = self._CreateNewtonRaphsonStrategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        return solution_strategy

    def _CreateLinearStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()
        calculate_norm_dx = False
        return KratosMultiphysics.ResidualBasedLinearStrategy(
            computing_model_part,
            time_scheme,
            builder_and_solver,
            self.settings["compute_reactions"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool(),
            calculate_norm_dx,
            self.settings["move_mesh_flag"].GetBool())

    def _CreateNewtonRaphsonStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        time_scheme = self._GetScheme()
        convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
            computing_model_part,
            time_scheme,
            convergence_criterion,
            builder_and_solver,
            self.settings["maximum_iterations"].GetInt(),
            self.settings["compute_reactions"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool(),
            self.settings["move_mesh_flag"].GetBool())
