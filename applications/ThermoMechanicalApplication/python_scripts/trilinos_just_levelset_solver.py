from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MetisApplication import *


# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part, settings):
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(settings.GetMeshVelocityVariable())
    model_part.AddNodalSolutionStepVariable(settings.GetUnknownVariable())
    model_part.AddNodalSolutionStepVariable(settings.GetConvectionVariable())
    # model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    # model_part.AddNodalSolutionStepVariable(settings.GetVolumeSourceVariable());
    # model_part.AddNodalSolutionStepVariable(settings.GetDensityVariable());
    # model_part.AddNodalSolutionStepVariable(settings.GetDiffusionVariable());
    # model_part.AddNodalSolutionStepVariable(HTC);
    # model_part.AddNodalSolutionStepVariable(DISTANCE);
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
    print("variables for the LEVELSET_SOLVER added correctly")


def AddDofs(model_part, settings):
    for node in model_part.Nodes:
        node.AddDof(settings.GetUnknownVariable())

    print("dofs for the LEVELSET_SOLVER added correctly")


class Solver:
    #

    def __init__(self, model_part, domain_size, my_settings):

        self.model_part = model_part

        self.time_scheme = TrilinosResidualBasedIncrementalUpdateStaticScheme()
        # self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
        self.Comm = CreateCommunicator()

        self.settings = my_settings
        self.domain_size = domain_size
        # definition of the solvers
        # self.linear_solver =  SkylineLUFactorizationSolver()
# self.linear_solver =SuperLUSolver()

        # pPrecond = DiagonalPreconditioner()
        # self.linear_solver = BICGSTABSolver(1e-5, 5000,pPrecond)

        self.buildertype = "standard"

        # definition of the solvers
        aztec_parameters = ParameterList()
        aztec_parameters.set("AZ_solver", "AZ_gmres")
        aztec_parameters.set("AZ_kspace", 200)
        aztec_parameters.set("AZ_output", "AZ_none")
        aztec_parameters.set("AZ_output", 10);
        preconditioner_type = "ILU"
        preconditioner_parameters = ParameterList()
        preconditioner_parameters.set("fact: drop tolerance", 1e-9);
        preconditioner_parameters.set("fact: level-of-fill", 1);
        overlap_level = 0
        nit_max = 1000
        linear_tol = 1e-6
        self.linear_solver = AztecSolver(aztec_parameters, preconditioner_type, preconditioner_parameters, linear_tol, nit_max, overlap_level);
        self.guess_row_size = 18

        self.dynamic_tau = 0.0

        self.echo_level = 0
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = False
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False

        # self.neigh_finder = FindNodalNeighboursProcess(self.model_part,9,18)
        # if (self.domain_size == 2):
           # self.elem_neighbor_finder = FindElementalNeighboursProcess(self.model_part, 2, 20)
        # else:
           # self.elem_neighbor_finder = FindElementalNeighboursProcess(self.model_part, 3, 20)

        self.Nmax = len(model_part.Properties)
        self.contact_matix = Matrix()

        # calculate normals
        self.normal_tools = BodyNormalCalculationUtils()

        # self.conv_criteria = ResidualCriteria(1e-3,1e-4)
        self.conv_criteria = TrilinosDisplacementCriteria(1e-3, 1e-4, self.Comm)
        self.max_iter = 5

        # material settings
        # self.rho_mat = 100.0
        # self.rho_empty = 1.0

        # self.specific_heat_mat = 1006.0
        # self.specific_heat_empty = 1.0

        # self.conductivity_mat = 0.024
        # self.conductivity_empty = 1.0
    #
    def Initialize(self):

        #(self.neigh_finder).ClearNeighbours();
        #(self.neigh_finder).Execute();
        #(self.elem_neighbor_finder).ClearNeighbours()
        #(self.elem_neighbor_finder).Execute()
        print("INSIDE INITIALIZE")
        (self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, self.settings)

        # self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,self.max_iter,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)
        import trilinos_strategy_python
        self.solver = trilinos_strategy_python.SolvingStrategyPython(self.buildertype, self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.Comm, self.guess_row_size)
        self.solver.max_iter = self.max_iter
        mpi.world.barrier()
        (self.solver).SetEchoLevel(self.echo_level)

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);

        if (self.domain_size == 2):
            self.normal_tools.CalculateBodyNormals(self.model_part, 2);
        else:
            self.normal_tools.CalculateBodyNormals(self.model_part, 3);

# print "Initialization monolithic solver finished"
    #
    def Solve(self):
        # print "*****************entering solve?????????????"
        #(self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS,self.settings)
        # self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);
        # self.ApplyFluidProperties()
        mpi.world.barrier()
        (self.solver).Solve()
# print "solving step monolithic solver finished"

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
