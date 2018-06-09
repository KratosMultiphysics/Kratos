from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(custom_settings):
    return MonolithicSolver(custom_settings)

#Base class to develop other solvers
class MonolithicSolver(object):
    """The base class for solid mechanics solvers.

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
    def __init__(self, custom_settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "dofs": [],
            "time_integration_settings":{
                "solution_type": "Dynamic",
  	        "analysis_type": "Non-linear",
                "time_integration": "Implicit",
                "integration_method": "Newmark",
                "time_integration_order": 1,
                "buffer_size": 2
            },
            "solving_strategy_settings":{
                "builder_type": "block_builder",
                "line_search": false,
                "implex": false,
                "compute_reactions": true,
                "move_mesh_flag": true,
                "clear_storage": false,
                "reform_dofs_at_each_step": false,
                "max_iteration": 10
            },
            "convergence_criterion_settings":{
                "convergence_criterion": "Residual_criterion",
                "variable_relative_tolerance": 1.0e-4,
                "variable_absolute_tolerance": 1.0e-9,
                "residual_relative_tolerance": 1.0e-4,
                "residual_absolute_tolerance": 1.0e-9,
                "separate_dofs": false
            },
            "linear_solver_settings":{
                "solver_type": "SuperLUSolver",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            }
        }
        """)

        # Check and fix supplied settings compatibility (experimental)
        from json_settings_utility import JsonSettingsUtility
        JsonSettingsUtility.CheckAndFixNotMatchingSettings(custom_settings,default_settings)

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # Validate and assign other values
        self.settings["time_integration_settings"].ValidateAndAssignDefaults(default_settings["time_integration_settings"])
        self.settings["solving_strategy_settings"].ValidateAndAssignDefaults(default_settings["solving_strategy_settings"])
        self.settings["convergence_criterion_settings"].ValidateAndAssignDefaults(default_settings["convergence_criterion_settings"])
        self.settings["linear_solver_settings"].ValidateAndAssignDefaults(default_settings["linear_solver_settings"])

        # Echo level
        self.echo_level = 0

    def GetMinimumBufferSize(self):
        buffer_size = self.settings["time_integration_settings"]["buffer_size"].GetInt()
        time_integration_order = self.settings["time_integration_settings"]["time_integration_order"].GetInt()
        if( buffer_size <= time_integration_order ):
            buffer_size = time_integration_order + 1
        return buffer_size;

    def SetComputingModelPart(self, computing_model_part):
        self.model_part = computing_model_part


    def ExecuteInitialize(self):

        # Main model part and computing model part
        self.main_model_part = self.model_part.GetRootModelPart()

        # Process information
        self.process_info = self.main_model_part.ProcessInfo

        # Add dofs
        if( self._is_not_restarted() ):
            self._add_dofs()

        # Create integration information (needed in other processes)
        self._set_time_integration_methods()

        # Set buffer
        if( self._is_not_restarted() ):
            self._set_and_fill_buffer()

    def ExecuteBeforeSolutionLoop(self):

        # The mechanical solver is created here if it does not already exist.
        mechanical_solver = self._get_mechanical_solver()
        mechanical_solver.SetEchoLevel(self.echo_level)

        self.Clear()

        self._check_initialized()

        self.Check()

        print("::[Solver]:: Ready")

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

    def Check(self):
        self._get_mechanical_solver().Check()


    #### Solve loop methods ####

    def Solve(self):
        self.Clear()
        self._get_mechanical_solver().Solve()

    # step by step:

    def InitializeSolutionStep(self):
        self.Clear()
        self._get_mechanical_solver().InitializeSolutionStep()

    def SolveSolutionStep(self):
        is_converged = self._get_mechanical_solver().SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self._get_mechanical_solver().FinalizeSolutionStep()


    #### Solver internal methods ####

    def _check_initialized(self):
        if( not self._is_not_restarted() ):
            self._get_solution_scheme().Initialize(self.main_model_part)
            if hasattr(mechanical_solver, 'SetInitializePerformedFlag'):
                self.get_mechanical_solver.SetInitializePerformedFlag(True)
            else:
                self.get_mechanical_solver.Set(KratosSolid.SolverLocalFlags.INITIALIZED, True)

    def _is_not_restarted(self):
        if( self.process_info.Has(KratosMultiphysics.IS_RESTARTED) ):
            if( self.process_info[KratosMultiphysics.IS_RESTARTED] == False ):
                return True
            else:
                return False
        else:
            return True

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


    def _set_and_fill_buffer(self):
        """Prepare nodal solution step data containers and time step information. """
        # Set the buffer size for the nodal solution steps data. Existing nodal
        # solution step data may be lost.
        buffer_size = self.settings["time_integration_settings"]["buffer_size"].GetInt()
        if buffer_size < self.GetMinimumBufferSize():
            buffer_size = self.GetMinimumBufferSize()
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
        self.process_info[KratosMultiphysics.IS_RESTARTED] = False

    def _create_solution_scheme(self):
        import schemes_factory
        Schemes = schemes_factory.SolutionScheme(self.settings["time_integration_settings"],self.settings["dofs"])
        mechanical_scheme = Schemes.GetSolutionScheme()
        return mechanical_scheme

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
            print("::[Solver]:: DOF's ADDED")


    def _set_scheme_parameters(self):
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
        self._set_scheme_parameters()

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
                print("::[Integration]:: ",scalar_integration_methods[dof],"(",dof,")")
            if( isinstance(kratos_variable,KratosMultiphysics.Array1DVariable3) and (not vector_dof_method_set) ):
                component_integration_methods[dof+"_X"].CalculateParameters(self.process_info)
                vector_dof_method_set = True
                print("::[Integration]:: ",component_integration_methods[dof+"_X"],"(",dof,")")

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
