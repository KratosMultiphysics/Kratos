from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(custom_settings, Model):
    return MonolithicSolver(Model, custom_settings)

#Base class to develop other solvers
class MonolithicSolver(object):
    """The base class for solid mechanics solvers

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Derived classes must override the function _create_mechanical_solver which
    constructs and returns a valid solving strategy. Depending on the type of
    solver, derived classes may also need to override the following functions:

    _create_solution_scheme
    _create_convergence_criterion
    _create_linear_solver
    _create_builder_and_solver
    _create_mechanical_solver

    The mechanical_solver, builder_and_solver, etc. should alway be retrieved
    using the getter functions _get_mechanical_solver, _get_builder_and_solver,
    etc. from this base class.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    settings -- Kratos parameters containing solver settings.
    model_part -- the model part used to construct the solver (computing_model_part).
    """
    def __init__(self, Model, custom_settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solving_model_part": "computing_domain",
            "dofs": [],
            "time_integration_settings":{
                "solution_type": "Dynamic",
  	        "analysis_type": "Non-linear",
                "time_integration": "Implicit",
                "integration_method": "Newmark",
                "time_integration_order": 1,
                "buffer_size": 2
            },
            "convergence_criterion_settings":{
                "convergence_criterion": "Residual_criterion",
                "variable_relative_tolerance": 1.0e-4,
                "variable_absolute_tolerance": 1.0e-9,
                "residual_relative_tolerance": 1.0e-4,
                "residual_absolute_tolerance": 1.0e-9,
                "separate_dofs": true
            },
            "solving_strategy_settings":{
                "builder_type": "block_builder",
                "line_search": false,
                "line_search_type": 0,
                "implex": false,
                "compute_reactions": true,
                "move_mesh_flag": true,
                "iterative_update": true,
                "clear_storage": false,
                "reform_dofs_at_each_step": false,
                "adaptive_solution": false,
                "max_iteration": 10
            },
            "linear_solver_settings":{
                "solver_type": "SuperLUSolver",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
            "processes": []
        }
        """)

        # Linear solver settings can have different number of settings
        if(custom_settings.Has("linear_solver_settings")):
            default_settings.RemoveValue("linear_solver_settings")
            default_settings.AddValue("linear_solver_settings",custom_settings["linear_solver_settings"])

        # Check and fix supplied settings compatibility (experimental)
        from json_settings_utility import JsonSettingsUtility
        JsonSettingsUtility.CheckAndFixNotMatchingSettings(custom_settings,default_settings)

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # Validate and assign other values
        self.settings["time_integration_settings"].ValidateAndAssignDefaults(default_settings["time_integration_settings"]) #default values in factory
        self.settings["solving_strategy_settings"].ValidateAndAssignDefaults(default_settings["solving_strategy_settings"])
        #self.settings["convergence_criterion_settings"].ValidateAndAssignDefaults(default_settings["convergence_criterion_settings"]) #default values in factory
        self.settings["linear_solver_settings"].ValidateAndAssignDefaults(default_settings["linear_solver_settings"])
        #print("Monolithic Solver Settings",self.settings.PrettyPrintJsonString())

        # Model
        self.model = Model

        # Echo level
        self.echo_level = 0


    def ExecuteInitialize(self):

        # Set model and info
        self._set_model_info()

        # Configure model and solver
        self._set_integration_parameters()

    def ExecuteBeforeSolutionLoop(self):

        # The mechanical solver is created here if it does not already exist.
        mechanical_solver = self._get_mechanical_solver()
        mechanical_solver.SetEchoLevel(self.echo_level)

        self.Clear()

        self._check_initialized()

        self.Check()

        print(self._class_prefix(),self.settings["time_integration_settings"]["integration_method"].GetString()+"_Scheme Ready")

    def GetVariables(self):

        import schemes_factory
        Schemes = schemes_factory.SolutionScheme(self.settings["time_integration_settings"],self.settings["dofs"])
        return Schemes.GetVariables()

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def SetEchoLevel(self, level):
        self.echo_level = level

    def Clear(self):
        if self.settings["solving_strategy_settings"]["clear_storage"].GetBool():
            self._get_mechanical_solver().Clear()
        self._check_reform_dofs()

    def Check(self):
        self._get_mechanical_solver().Check()


    #### Solve loop methods ####

    def Solve(self):
        self.Clear()
        return self._get_mechanical_solver().Solve()

    # step by step:

    def InitializeSolutionStep(self):
        self.Clear()
        self._get_mechanical_solver().InitializeSolutionStep()

    def SolveSolutionStep(self):
        return self._get_mechanical_solver().SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self._get_mechanical_solver().FinalizeSolutionStep()


    #### Solver internal methods ####

    def _check_adaptive_solution(self):
        return self.settings["solving_strategy_settings"]["adaptive_solution"].GetBool()

    def _check_reform_dofs(self):
        if self._domain_parts_updated():
            if not self._get_mechanical_solver().GetOptions().Is(KratosSolid.SolverLocalFlags.REFORM_DOFS):
                self._get_mechanical_solver().GetOptions().Set(KratosSolid.SolverLocalFlags.REFORM_DOFS, True)
                KratosMultiphysics.Logger.PrintInfo("", self._class_prefix()+" Set flag (REFORM_DOFS:true)")
        else:
            if self._get_mechanical_solver().GetOptions().Is(KratosSolid.SolverLocalFlags.REFORM_DOFS):
                self._get_mechanical_solver().GetOptions().Set(KratosSolid.SolverLocalFlags.REFORM_DOFS, False)
                KratosMultiphysics.Logger.PrintInfo("", self._class_prefix()+" Set flag (REFORM_DOFS:false)")

    def _domain_parts_updated(self):
        update_time = False
        if not self._is_not_restarted():
            if self.process_info.Has(KratosSolid.RESTART_STEP_TIME):
                update_time = self._check_current_time_step(self.process_info[KratosSolid.RESTART_STEP_TIME])

        if not update_time and self.process_info.Has(KratosSolid.MESHING_STEP_TIME):
            update_time = self._check_previous_time_step(self.process_info[KratosSolid.MESHING_STEP_TIME])

        if not update_time and self.process_info.Has(KratosSolid.CONTACT_STEP_TIME):
            update_time = self._check_previous_time_step(self.process_info[KratosSolid.CONTACT_STEP_TIME])

        return update_time

    def _check_current_time_step(self, step_time):
        current_time  = self.process_info[KratosMultiphysics.TIME]
        delta_time    = self.process_info[KratosMultiphysics.DELTA_TIME]
        #arithmetic floating point tolerance
        tolerance = delta_time * 0.001

        if( step_time > current_time-tolerance and step_time < current_time+tolerance ):
            return True
        else:
            return False

    def _check_previous_time_step(self, step_time):
        current_time  = self.process_info[KratosMultiphysics.TIME]
        delta_time    = self.process_info[KratosMultiphysics.DELTA_TIME]
        previous_time = current_time - delta_time

        #arithmetic floating point tolerance
        tolerance = delta_time * 0.001

        if( step_time > previous_time-tolerance and step_time < previous_time+tolerance ):
            return True
        else:
            return False

    def _check_initialized(self):
        if not self._is_not_restarted():
            mechanical_solver = self._get_mechanical_solver()
            if hasattr(mechanical_solver, 'SetInitializePerformedFlag'):
                mechanical_solver.SetInitializePerformedFlag(True)
            else:
                mechanical_solver.Set(KratosSolid.SolverLocalFlags.INITIALIZED, True)

    def _is_not_restarted(self):
        if self.process_info.Has(KratosMultiphysics.IS_RESTARTED):
            if self.process_info[KratosMultiphysics.IS_RESTARTED]:
                return False
            else:
                return True
        else:
            return True

    def _set_model_info(self):

        # Get solving model part
        self.model_part = self.model[self.settings["solving_model_part"].GetString()]

        # Main model part from computing model part
        self.main_model_part = self.model_part.GetRootModelPart()

        # Process information
        self.process_info = self.main_model_part.ProcessInfo

    def _set_integration_parameters(self):
        # Add dofs
        if( self._is_not_restarted() ):
            self._add_dofs()

        # Create integration information (needed in other processes)
        self._set_time_integration_methods()

        # Set buffer
        if( self._is_not_restarted() ):
            self._set_and_fill_buffer()

    def _get_solution_scheme(self):
        if not hasattr(self, '_solution_scheme'):
            self._solution_scheme = self._create_solution_scheme()
        return self._solution_scheme

    def _get_convergence_criterion(self):
        if not hasattr(self, '_convergence_criterion'):
            self._convergence_criterion = self._create_convergence_criterion()
        return self._convergence_criterion

    def _get_linear_solver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._create_linear_solver()
        return self._linear_solver

    def _get_builder_and_solver(self):
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._create_builder_and_solver()
        return self._builder_and_solver

    def _get_mechanical_solver(self):
        if not hasattr(self, '_mechanical_solver'):
            self._mechanical_solver = self._create_mechanical_solver()
        return self._mechanical_solver

    def _get_minimum_buffer_size(self):
        buffer_size = self.settings["time_integration_settings"]["buffer_size"].GetInt()
        time_integration_order = self.settings["time_integration_settings"]["time_integration_order"].GetInt()
        if( buffer_size <= time_integration_order ):
            buffer_size = time_integration_order + 1
        return buffer_size;

    def _set_and_fill_buffer(self):
        """Prepare nodal solution step data containers and time step information. """
        # Set the buffer size for the nodal solution steps data. Existing nodal
        # solution step data may be lost.
        buffer_size = self._get_minimum_buffer_size()
        self.main_model_part.SetBufferSize(buffer_size)
        # Cycle the buffer. This sets all historical nodal solution step data to
        # the current value and initializes the time stepping in the process info.
        delta_time = self.process_info[KratosMultiphysics.DELTA_TIME]
        time = self.process_info[KratosMultiphysics.TIME]
        step =-buffer_size
        time = time - delta_time * buffer_size
        self.process_info.SetValue(KratosMultiphysics.TIME, time)
        for i in range(0, buffer_size):
            step = step + 1
            time = time + delta_time
            self.process_info.SetValue(KratosMultiphysics.STEP, step)
            self.main_model_part.CloneTimeStep(time)

    def _create_solution_scheme(self):
        import schemes_factory
        Schemes = schemes_factory.SolutionScheme(self.settings["time_integration_settings"],self.settings["dofs"])
        solution_scheme = Schemes.GetSolutionScheme()
        solution_scheme.SetProcessVector(self._get_scheme_custom_processes())
        return solution_scheme

    def _create_convergence_criterion(self):
        criterion_parameters = self.settings["convergence_criterion_settings"]
        criterion_parameters.AddEmptyValue("echo_level").SetInt(self.echo_level)
        # Construction of the class convergence_criterion
        import convergence_criteria_factory
        convergence_criterion = convergence_criteria_factory.ConvergenceCriterion(criterion_parameters,self.settings["dofs"])
        return convergence_criterion.GetConvergenceCriterion()

    def _create_linear_solver(self):
        import linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        return linear_solver

    def _create_builder_and_solver(self):
        linear_solver = self._get_linear_solver()
        if(self.settings["solving_strategy_settings"]["builder_type"].GetString() == "block_builder"):
            # To keep matrix blocks in builder
            builder_and_solver = KratosSolid.BlockBuilderAndSolver(linear_solver)
        else:
            builder_and_solver = KratosSolid.ReductionBuilderAndSolver(linear_solver)

        return builder_and_solver

    def _create_mechanical_solver(self):
        raise Exception("please implement the Custom Choice of your Mechanical Solver (_create_mechanical_solver) in your solver")


    def _get_dofs(self):
        import schemes_factory
        Schemes = schemes_factory.SolutionScheme(self.settings["time_integration_settings"],self.settings["dofs"])
        return Schemes.GetDofsAndReactions()

    def _add_dofs(self):
        dof_variables, dof_reactions = self._get_dofs()
        AddDofsProcess = KratosSolid.AddDofsProcess(self.main_model_part, dof_variables, dof_reactions)
        AddDofsProcess.Execute()
        if( self.echo_level > 1 ):
            print(dof_variables + dof_reactions)
            print(self._class_prefix()+" DOF's ADDED")

    @classmethod
    def _get_scheme_custom_processes(self):
        process_list = []
        return process_list

    def _set_scheme_process_info_parameters(self):
        pass

    def _get_time_integration_methods(self):

        # set solution scheme
        import schemes_factory
        Schemes = schemes_factory.SolutionScheme(self.settings["time_integration_settings"],self.settings["dofs"])

        # set integration method dictionary for scalars
        scalar_integration_methods = Schemes.GetScalarIntegrationMethods()

        # set integration method dictionary for components
        component_integration_methods = Schemes.GetComponentIntegrationMethods()

        #print(scalar_integration_methods)
        #print(component_integration_methods)

        # set time order
        self.process_info[KratosSolid.TIME_INTEGRATION_ORDER] = Schemes.GetTimeIntegrationOrder()

        # set scheme parameters to process_info
        self._set_scheme_process_info_parameters()

        # first: calculate parameters (only once permitted for each monolithic dof set "components + scalar")
        dofs_list = self.settings["dofs"]

        dofs = []
        for i in range(0, dofs_list.size() ):
            dofs.append(dofs_list[i].GetString())

        # add default DISPLACEMENT dof
        if( len(dofs) == 0 or (len(dofs) == 1 and dofs[0] =="ROTATION") ):
            dofs.append('DISPLACEMENT')

        #print(" DOFS ",dofs)

        scalar_dof_method_set = False
        vector_dof_method_set = False

        for dof in dofs:
            kratos_variable = KratosMultiphysics.KratosGlobals.GetVariable(dof)
            if( isinstance(kratos_variable,KratosMultiphysics.DoubleVariable) and (not scalar_dof_method_set) ):
                scalar_integration_methods[dof].CalculateParameters(self.process_info)
                scalar_dof_method_set = True

                print("::[----Integration----]::",scalar_integration_methods[dof],"("+dof+")")
            if( isinstance(kratos_variable,KratosMultiphysics.Array1DVariable3) and (not vector_dof_method_set) ):
                component_integration_methods[dof+"_X"].CalculateParameters(self.process_info)
                vector_dof_method_set = True
                print("::[----Integration----]::",component_integration_methods[dof+"_X"],"("+dof+")")

        return scalar_integration_methods, component_integration_methods


    def _set_time_integration_methods(self):

        scalar_integration_methods, component_integration_methods = self._get_time_integration_methods();

        # second: for the same method the parameters (already calculated)
        scalar_integration_methods_container = KratosSolid.ScalarTimeIntegrationMethods()
        for dof, method in scalar_integration_methods.items():
            method.SetParameters(self.process_info) #set same parameters to all methods from process_info values
            scalar_integration_methods_container.Set(dof,method)


        component_integration_methods_container = KratosSolid.ComponentTimeIntegrationMethods()
        for dof, method in component_integration_methods.items():
            method.SetParameters(self.process_info) #set same parameters to all methods from process_info values
            component_integration_methods_container.Set(dof,method)

        # set time integration methods (for scalars and components) to process_info for processes access
        if( len(scalar_integration_methods) ):
            scalar_integration_methods_container.AddToProcessInfo(KratosSolid.SCALAR_TIME_INTEGRATION_METHODS, scalar_integration_methods_container, self.process_info)

        component_integration_methods_container.AddToProcessInfo(KratosSolid.COMPONENT_TIME_INTEGRATION_METHODS, component_integration_methods_container, self.process_info)

        #print(scalar_integration_methods)
        #print(component_integration_methods)

    #
    @classmethod
    def _class_prefix(self):
        header = "::[-Monolithic_Solver-]::"
        return header
