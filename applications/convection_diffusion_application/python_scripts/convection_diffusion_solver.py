from __future__ import print_function, absolute_import, division

from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):
    return ConvectionDiffusionSolver(main_model_part, custom_settings)

class ConvectionDiffusionSolver(object):
    ''' Python class for the convection-diffusion solver.
    '''

    def __init__(self, model_part, custom_settings):

        ##settings string in json format
        default_settings = Parameters("""
        {
            "solver_type": "convection_diffusion_solver",
            "dynamic_tau": 0.0,
            "echo_level": 1,
            "calculate_reactions": false,
            "reform_dofs_at_each_step": false,
            "time_stepping"                : {
                "theta"      : 0.5,
                "time_step"  : 0.1
            }
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.solver_parameters = custom_settings
        self.solver_parameters.ValidateAndAssignDefaults(default_settings)

        self.main_model_part = model_part

        if not self.main_model_part.ProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS):
            raise Exception("the provided model_part does not have CONVECTION_DIFFUSION_SETTINGS defined.")

        self.linear_solver = None
        self.strategy = None

        self.echo_level = self.solver_parameters["echo_level"].GetInt()

        self.time_step = self.solver_parameters["time_stepping"]["time_step"].GetDouble()

    def AddVariables(self):
        ''' Add nodal solution step variables based on provided CONVECTION_DIFFUSION_SETTINGS
        '''
        if self.main_model_part.ProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS):
            settings = self.main_model_part.ProcessInfo[CONVECTION_DIFFUSION_SETTINGS]

            if settings.IsDefinedUnknownVariable():
                self.main_model_part.AddNodalSolutionStepVariable(settings.GetUnknownVariable())
            if settings.IsDefinedVelocityVariable():
                self.main_model_part.AddNodalSolutionStepVariable(settings.GetVelocityVariable())
            if settings.IsDefinedMeshVelocityVariable():
                self.main_model_part.AddNodalSolutionStepVariable(settings.GetMeshVelocityVariable())
            if settings.IsDefinedDiffusionVariable():
                self.main_model_part.AddNodalSolutionStepVariable(settings.GetDiffusionVariable())
            if settings.IsDefinedSpecificHeatVariable():
                self.main_model_part.AddNodalSolutionStepVariable(settings.GetSpecificHeatVariable())
            if settings.IsDefinedDensityVariable():
                self.main_model_part.AddNodalSolutionStepVariable(settings.GetDensityVariable())
            if settings.IsDefinedVolumeSourceVariable():
                self.main_model_part.AddNodalSolutionStepVariable(settings.GetVolumeSourceVariable())
            if settings.IsDefinedSurfaceSourceVariable():
                self.main_model_part.AddNodalSolutionStepVariable(settings.GetSurfaceSourceVariable())
            if settings.IsDefinedReactionVariable():
                self.main_model_part.AddNodalSolutionStepVariable(settings.GetReactionVariable())
        else:
            raise Exception("the provided model_part does not have CONVECTION_DIFFUSION_SETTINGS defined.")

    def AddDofs(self):
        ''' Add nodal degrees of freedom based on provided CONVECTION_DIFFUSION_SETTINGS
        '''
        if not self.main_model_part.ProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS):
            raise Exception("the provided model_part does not have CONVECTION_DIFFUSION_SETTINGS defined.")

        settings = self.main_model_part.ProcessInfo[CONVECTION_DIFFUSION_SETTINGS]

        if not settings.IsDefinedUnknownVariable():
            raise Exception("the provided CONVECTION_DIFFUSION_SETTINGS does not define an Unknown variable.")

        unknown_variable = settings.GetUnknownVariable()

        if settings.IsDefinedReactionVariable():
            reaction_variable = settings.GetReactionVariable()

            for node in self.main_model_part.Nodes:
                node.AddDof(unknown_variable,reaction_variable)
        else:
            for node in self.main_model_part.Nodes:
                node.AddDof(unknown_variable)

        Logger.PrintInfo("Convection-diffusion solver","DOFs for the convection diffusion solver added correctly")

    def ImportModelPart(self):
        pass # the model part is read by the fluid solver/passed from outside

    def PrepareModelPart(self):
        # Duplicate model part
        self.thermal_model_part = ModelPart("Thermal")
        if self.main_model_part.ProcessInfo[DOMAIN_SIZE] == 2:
            conv_diff_element = "EulerianConvDiff2D"
            conv_diff_condition = "Condition2D2N"
        elif self.main_model_part.ProcessInfo[DOMAIN_SIZE] == 3:
            conv_diff_element = "EulerianConvDiff3D"
            conv_diff_condition = "Condition3D3N"

        modeler = ConnectivityPreserveModeler()
        modeler.GenerateModelPart(self.main_model_part,self.thermal_model_part,conv_diff_element,conv_diff_condition)


    def Initialize(self):
        ''' Initialize the underlying C++ objects and validate input
        '''

        if self.solver_parameters.Has("linear_solver_settings"):
            import linear_solver_factory
            self.linear_solver = linear_solver_factory.ConstructSolver(self.solver_parameters["linear_solver_settings"])
        else:
            diagonal_preconditioner = DiagonalPreconditioner()
            linear_solver_tolerance = 1e-9
            linear_solver_maximum_iterations = 5000
            self.linear_solver = BICGSTABSolver(linear_solver_tolerance,
                                                linear_solver_maximum_iterations,
                                                diagonal_preconditioner)

        if self.strategy is None:

            scheme = ResidualBasedIncrementalUpdateStaticScheme()
            builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.linear_solver)

            calculate_norm_dx = False
            move_mesh = False

            self.strategy = ResidualBasedLinearStrategy(
                self.thermal_model_part,
                scheme,
                self.linear_solver,
                builder_and_solver,
                self.solver_parameters["calculate_reactions"].GetBool(),
                self.solver_parameters["reform_dofs_at_each_step"].GetBool(),
                calculate_norm_dx,
                move_mesh)

        self.strategy.SetEchoLevel(self.echo_level)
        self.strategy.Check()

        verbose = True
        self._ValidateInput(verbose)

        #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.thermal_model_part.ProcessInfo[THETA] = self.solver_parameters["time_stepping"]["theta"].GetDouble()
        self.thermal_model_part.ProcessInfo[DYNAMIC_TAU] = self.solver_parameters["dynamic_tau"].GetDouble()

    def Solve(self):
        ''' Solve an iteration of the convection-diffusion problem
        '''
        self.strategy.Solve()

    def AdvanceInTime(self, current_time):
        dt = self.time_step
        new_time = current_time + dt

        self.thermal_model_part.CloneTimeStep(new_time)
        self.thermal_model_part.ProcessInfo[STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self.strategy.InitializeSolutionStep()

    def Predict(self):
        if self._TimeBufferIsInitialized():
            self.strategy.Predict()

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            is_converged = self.strategy.SolveSolutionStep()
            if not is_converged and self._IsPrintingRank():
                msg  = "Not converge for iteration " + str(self.thermal_model_part.ProcessInfo[STEP]) + "\n"
                msg += "corresponding to time " + str(self.thermal_model_part.ProcessInfo[TIME]) + "\n"
                Logger.PrintWarning("Convection-diffusion solver",msg)

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            (self.strategy).FinalizeSolutionStep()

    def _TimeBufferIsInitialized(self):
        # This is a Backward Euler/Crank-Nicolson element, only one old step is required.
        # We always have one extra old step (step 0, read from input)
        return self.thermal_model_part.ProcessInfo[STEP] + 1 >= 2


    def _ValidateInput(self,verbose=False):
        ''' Verify that the convection-diffusion settings have the required variables
        '''

        if not self.thermal_model_part.ProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS):
            raise Exception("the provided model_part does not have CONVECTION_DIFFUSION_SETTINGS defined.")

        settings = self.thermal_model_part.ProcessInfo[CONVECTION_DIFFUSION_SETTINGS]

        # hackish way to get first node in model part
        ref_node = next( self.thermal_model_part.Nodes.__iter__() )

        # Density
        if settings.IsDefinedDensityVariable():
            if not ref_node.SolutionStepsDataHas( settings.GetDensityVariable() ):
                raise Exception("Density variable defined in CONVECTION_DIFFUSION_SETTINGS but not added to the nodes.")
        else:
            if verbose:
                print("Density variable not defined in CONVECTION_DIFFUSION_SETTINGS. Assuming unit density.")

        # Diffusion
        if settings.IsDefinedDiffusionVariable():
            if not ref_node.SolutionStepsDataHas( settings.GetDiffusionVariable() ):
                raise Exception("Diffusion variable defined in CONVECTION_DIFFUSION_SETTINGS but not added to the nodes.")
        else:
            if verbose:
                print("Diffusion variable not defined in CONVECTION_DIFFUSION_SETTINGS. Assuming pure convection.")

        # Unknown variable
        if settings.IsDefinedUnknownVariable():
            if not ref_node.SolutionStepsDataHas( settings.GetUnknownVariable() ):
                raise Exception("Unknown variable defined in CONVECTION_DIFFUSION_SETTINGS but not added to the nodes.")
        else:
            raise Exception("Unknown variable not defined in CONVECTION_DIFFUSION_SETTINGS.")

        # Volume source
        if settings.IsDefinedVolumeSourceVariable():
            if not ref_node.SolutionStepsDataHas( settings.GetVolumeSourceVariable() ):
                raise Exception("Volume source variable defined in CONVECTION_DIFFUSION_SETTINGS but not added to the nodes.")
        else:
            if verbose:
                print("Volume source variable not defined in CONVECTION_DIFFUSION_SETTINGS. Assuming zero source term.")

        # Surface source
        #if settings.IsDefinedSurfaceSourceVariable():
        #    if not ref_node.SolutionStepsDataHas( settings.GetSurfaceSourceVariable() ):
        #        raise Exception("Surface source variable defined in CONVECTION_DIFFUSION_SETTINGS but not added to the nodes.")
        #else:
        #    if verbose:
        #        print("Surface source variable not defined in CONVECTION_DIFFUSION_SETTINGS. Assuming zero source term.")

        # Velocity
        if settings.IsDefinedVelocityVariable():
            if not ref_node.SolutionStepsDataHas( settings.GetVelocityVariable() ):
                raise Exception("Velocity variable defined in CONVECTION_DIFFUSION_SETTINGS but not added to the nodes.")
        else:
            if verbose:
                print("Velocity variable not defined in CONVECTION_DIFFUSION_SETTINGS. Assuming zero velocity.")

        # Mesh Velocity
        if settings.IsDefinedMeshVelocityVariable():
            if not ref_node.SolutionStepsDataHas( settings.GetMeshVelocityVariable() ):
                raise Exception("Mesh velocity variable defined in CONVECTION_DIFFUSION_SETTINGS but not added to the nodes.")
        else:
            if verbose:
                print("Mesh velocity variable not defined in CONVECTION_DIFFUSION_SETTINGS. Assuming zero mesh velocity.")

        # Specific heat variable
        if settings.IsDefinedSpecificHeatVariable():
            if not ref_node.SolutionStepsDataHas( settings.GetSpecificHeatVariable() ):
                raise Exception("Specific heat variable defined in CONVECTION_DIFFUSION_SETTINGS but not added to the nodes.")
        else:
            if verbose:
                print("Specific heat variable not defined in CONVECTION_DIFFUSION_SETTINGS. Assuming unit specific heat.")
