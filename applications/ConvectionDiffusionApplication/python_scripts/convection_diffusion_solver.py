from __future__ import print_function, absolute_import, division

from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionApplication import *

def AddVariables(model_part):
    ''' Add nodal solution step variables based on provided CONVECTION_DIFFUSION_SETTINGS
    '''

    if model_part.ProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS):
        settings = model_part.ProcessInfo[CONVECTION_DIFFUSION_SETTINGS]

        if settings.IsDefinedUnknownVariable():
            model_part.AddNodalSolutionStepVariable(settings.GetUnknownVariable())
        if settings.IsDefinedVelocityVariable():
            model_part.AddNodalSolutionStepVariable(settings.GetVelocityVariable())
        if settings.IsDefinedMeshVelocityVariable():
            model_part.AddNodalSolutionStepVariable(settings.GetMeshVelocityVariable())
        if settings.IsDefinedDiffusionVariable():
            model_part.AddNodalSolutionStepVariable(settings.GetDiffusionVariable())
        if settings.IsDefinedSpecificHeatVariable():
            model_part.AddNodalSolutionStepVariable(settings.GetSpecificHeatVariable())
        if settings.IsDefinedDensityVariable():
            model_part.AddNodalSolutionStepVariable(settings.GetDensityVariable())
        if settings.IsDefinedVolumeSourceVariable():
            model_part.AddNodalSolutionStepVariable(settings.GetVolumeSourceVariable())
        if settings.IsDefinedSurfaceSourceVariable():
            model_part.AddNodalSolutionStepVariable(settings.GetSurfaceSourceVariable())
        if settings.IsDefinedReactionVariable():
            model_part.AddNodalSolutionStepVariable(settings.GetReactionVariable())
    else:
        raise Exception("the provided model_part does not have CONVECTION_DIFFUSION_SETTINGS defined.")

def AddDofs(model_part):
    ''' Add nodal degrees of freedom based on provided CONVECTION_DIFFUSION_SETTINGS
    '''
    if not model_part.ProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS):
        raise Exception("the provided model_part does not have CONVECTION_DIFFUSION_SETTINGS defined.")

    settings = model_part.ProcessInfo[CONVECTION_DIFFUSION_SETTINGS]

    if not settings.IsDefinedUnknownVariable():
        raise Exception("the provided CONVECTION_DIFFUSION_SETTINGS does not define an Unknown variable.")

    unknown_variable = settings.GetUnknownVariable()

    if settings.IsDefinedReactionVariable():
        reaction_variable = settings.GetReactionVariable()

        for node in model_part.Nodes:
            node.AddDof(unknown_variable,reaction_variable)
    else:
        for node in model_part.Nodes:
            node.AddDof(unknown_variable)

    print("DOFs for the convection diffusion solver added correctly")

class ConvectionDiffusionSolver(object):
    ''' Python class for the convection-diffusion solver.
    '''

    def __init__(self, model_part, domain_size):

        self.model_part = model_part
        self.domain_size = domain_size

        #Variable defining the temporal scheme (0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
        self.theta = 0.5
        self.dynamic_tau = 0.0

        if not self.model_part.ProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS):
            raise Exception("the provided model_part does not have CONVECTION_DIFFUSION_SETTINGS defined.")

        self.linear_solver = None
        self.strategy = None

        self.calculate_reactions = False
        self.reform_dofs_at_each_step = False

        self.echo_level = 1


    def Initialize(self):
        ''' Initialize the underlying C++ objects and validate input
        '''

        if self.linear_solver is None:
            diagonal_preconditioner = DiagonalPreconditioner()
            linear_solver_tolerance = 1e-9
            linear_solver_maximum_iterations = 5000
            self.linear_solver = BICGSTABSolver(linear_solver_tolerance,
                                                linear_solver_maximum_iterations,
                                                diagonal_preconditioner)

        if self.strategy is None:

            scheme = ResidualBasedIncrementalUpdateStaticScheme()
            builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.linear_solver)

            calculate_norm_dx = True
            move_mesh = False

            self.strategy = ResidualBasedLinearStrategy(
                self.model_part,
                scheme,
                builder_and_solver,
                self.calculate_reactions,
                self.reform_dofs_at_each_step,
                calculate_norm_dx,
                move_mesh)

        self.strategy.SetEchoLevel(self.echo_level)
        self.strategy.Check()

        verbose = True
        self.ValidateInput(verbose)

        self.model_part.ProcessInfo[TIME_INTEGRATION_THETA] = self.theta
        self.model_part.ProcessInfo[DYNAMIC_TAU] = self.dynamic_tau

    def Solve(self):
        ''' Solve an iteration of the convection-diffusion problem
        '''
        self.strategy.Solve()

    def ValidateInput(self,verbose=False):
        ''' Verify that the convection-diffusion settings have the required variables
        '''

        if not self.model_part.ProcessInfo.Has(CONVECTION_DIFFUSION_SETTINGS):
            raise Exception("the provided model_part does not have CONVECTION_DIFFUSION_SETTINGS defined.")

        settings = self.model_part.ProcessInfo[CONVECTION_DIFFUSION_SETTINGS]

        # hackish way to get first node in model part
        ref_node = next( self.model_part.Nodes.__iter__() )

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
