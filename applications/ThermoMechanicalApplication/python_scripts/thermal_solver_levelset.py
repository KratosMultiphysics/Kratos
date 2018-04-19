from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

from KratosMultiphysics.StructuralApplication import *


# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part, settings):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(settings.GetMeshVelocityVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetUnknownVariable());
    model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    model_part.AddNodalSolutionStepVariable(settings.GetVolumeSourceVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDensityVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDiffusionVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetSurfaceSourceVariable());
    # model_part.AddNodalSolutionStepVariable(CONVECTION_COEFFICIENT);
    model_part.AddNodalSolutionStepVariable(NORMAL);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(settings.GetTransferCoefficientVariable());
    model_part.AddNodalSolutionStepVariable(HTC);
    model_part.AddNodalSolutionStepVariable(ENTHALPY);
    model_part.AddNodalSolutionStepVariable(SOLIDFRACTION);
    model_part.AddNodalSolutionStepVariable(SOLIDFRACTION_RATE);
    model_part.AddNodalSolutionStepVariable(DISTANCE);
    model_part.AddNodalSolutionStepVariable(NODAL_PAUX);
    model_part.AddNodalSolutionStepVariable(NODAL_VOLUME);
    model_part.AddNodalSolutionStepVariable(SOLIDIF_TIME);
    model_part.AddNodalSolutionStepVariable(SOLIDIF_MODULUS);
    model_part.AddNodalSolutionStepVariable(IS_VISITED);
    model_part.AddNodalSolutionStepVariable(SOLIDFRACTION);
    print("variables for the THERMAL_SOLVER added correctly")


def AddDofs(model_part, settings):
    for node in model_part.Nodes:
        node.AddDof(settings.GetUnknownVariable());

    print("dofs for the THERMAL_SOLVER added correctly")


class Solver:
    #

    def __init__(self, model_part, domain_size, my_settings):

        self.model_part = model_part

        # self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
        self.time_scheme = ResidualBasedIncrementalUpdateStaticVariablePropertyScheme()
        self.settings = my_settings
        self.domain_size = domain_size
        # definition of the solvers
        # self.linear_solver =  SkylineLUFactorizationSolver()
# self.linear_solver =SuperLUSolver()

# pPrecond = DiagonalPreconditioner()
# pPrecond = ILU0Preconditioner()
        # self.linear_solver =  BICGSTABSolver(1e-6, 5000,pPrecond)
       # self.linear_solver = BICGSTABSolver(1e-4, 5000,pPrecond)

        # new solvers
        gmres_size = 50
        tol = 1e-3
        verbosity = 0
        self.linear_solver = AMGCLSolver(AMGCLSmoother.ILU0, AMGCLIterativeSolverType.GMRES, tol, 200, verbosity, gmres_size)

        self.dynamic_tau = 0.0

        self.echo_level = 0
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = False
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False

        self.neigh_finder = FindNodalNeighboursProcess(self.model_part, 9, 18)
        if (self.domain_size == 2):
            self.elem_neighbor_finder = FindElementalNeighboursProcess(self.model_part, 2, 20)
        else:
            self.elem_neighbor_finder = FindElementalNeighboursProcess(self.model_part, 3, 20)

        self.Nmax = len(model_part.Properties)
        self.contact_matix = Matrix()

        # calculate normals
        self.normal_tools = BodyNormalCalculationUtils()

        # self.conv_criteria = ResidualCriteria(1e-3,1e-4)
        self.conv_criteria = IncrementalDisplacementCriteria(1e-4, 1e-6)
        self.max_iter = 3

        BiphasicFillingUtilities().ComputeNodalVolume(self.model_part)

        # linesearch solver
        # definition of parameters
        self.MaxLineSearchIterations = 20
        self.tolls = 0.8           # energy tolerance factor on LineSearch (0.8 is ok)
        self.amp = 1.618         # maximum amplification factor
        self.etmxa = 3.0           # maximum allowed step length
        self.etmna = 0.1           # minimum allowed step length
        self.toler = 1.0E-9
        self.norm = 1.0E-6
        self.ApplyLineSearches = True

    #
    def Initialize(self):

        if (self.domain_size == 2):
            self.duplicate_and_create_conditions = DuplicateInterfaceNodesCreateConditionsProcess(self.model_part, "HeatContact2D", self.Nmax, self.contact_matix)
        else:
            print("duplicate_and_create_conditions is not implemented for 3D")

        (self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, self.settings)

        self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,self.MaxNewtonRapshonIterations,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)

        #self.solver = ResidualBasedNewtonRaphsonLineSearchesStrategy(self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria,
                                                                     #self.MaxNewtonRapshonIterations, self.MaxLineSearchIterations, self.tolls, self.amp, self.etmxa, self.etmna,
                                                                     #self.CalculateReactionFlag,
                                                                     #self.ReformDofSetAtEachStep,
                                                                     #self.MoveMeshFlag,
                                                                    #self.ApplyLineSearches)

        (self.solver).SetEchoLevel(self.echo_level)
        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);


    #
    def Solve(self):
        (self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, self.settings)
        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);
# self.ApplyFluidProperties()
        (self.solver).Solve()
# print "solving step monolithic solver finished"

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
