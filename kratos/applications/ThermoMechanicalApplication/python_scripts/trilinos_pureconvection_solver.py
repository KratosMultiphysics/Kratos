from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.mpi import *
# try:
 # import boost.mpi as mpi
# except ImportError:
 # import mpi


def AddVariables(model_part, settings):
    model_part.AddNodalSolutionStepVariable(settings.GetConvectionVariable())
    model_part.AddNodalSolutionStepVariable(settings.GetMeshVelocityVariable())
    model_part.AddNodalSolutionStepVariable(settings.GetUnknownVariable())

    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)

    if(mpi.rank == 0):
        print("variables for the trilinos thermal solver added correctly")


def AddDofs(model_part, settings):
    for node in model_part.Nodes:
        node.AddDof(settings.GetUnknownVariable())

    if(mpi.rank == 0):
        print("dofs for the trilinos thermal solver solver added correctly")


class Solver:
    #

    def __init__(self, model_part, domain_size, my_settings):

        self.model_part = model_part
        self.domain_size = domain_size
#        self.settings = my_settings
        (self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, my_settings)

        self.time_scheme = TrilinosResidualBasedIncrementalUpdateStaticScheme()

        self.Comm = CreateCommunicator()

        self.buildertype = "standard"

        # definition of the solvers
        aztec_parameters = ParameterList()
        aztec_parameters.set("AZ_solver", "AZ_gmres")
        aztec_parameters.set("AZ_kspace", 200)
        aztec_parameters.set("AZ_output", "AZ_none")

        preconditioner_type = "ILU"
        preconditioner_parameters = ParameterList()
        preconditioner_parameters.set("fact: drop tolerance", 1e-9)
        preconditioner_parameters.set("fact: level-of-fill", 1);
        overlap_level = 0
        nit_max = 1000
        linear_tol = 1e-9
        self.linear_solver = AztecSolver(aztec_parameters, preconditioner_type, preconditioner_parameters, linear_tol, nit_max, overlap_level);

#        aztec_parameters = ParameterList()
#        aztec_parameters.set("AZ_solver","AZ_gmres");
#        aztec_parameters.set("AZ_kspace", 100);
# aztec_parameters.set("AZ_output","AZ_none");
#        aztec_parameters.set("AZ_output",10);
#        MLList = ParameterList()
#        default_settings = EpetraDefaultSetter()
#        default_settings.SetDefaults(MLList,"NSSA");
# MLList.set("ML output", 10);
#        MLList.set("PDE equations", 1);
#        MLList.setboolvalue("null space: add default vectors",True);
#        MLList.set("aggregation: type","Uncoupled");
#        nit_max = 1000
#        linear_tol = 1e-9
#        print MLList
#        self.linear_solver =  MultiLevelSolver(aztec_parameters,MLList,linear_tol,nit_max);

        # definition of the convergence criteria
        self.conv_criteria = TrilinosDisplacementCriteria(1e-9, 1e-12, self.Comm)

        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = False
        self.MoveMeshFlag = False
        self.calculate_norm_dx_flag = False
        self.max_iterations = 4

        if(domain_size == 2):
            self.guess_row_size = 20
        else:
            self.guess_row_size = 45

        self.dynamic_tau = 1.0

    #
    def Initialize(self):
# 	(self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS,self.settings)

        print(self.dynamic_tau)
        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);

        import trilinos_strategy_python
        self.solver = trilinos_strategy_python.SolvingStrategyPython(self.buildertype, self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.Comm, self.guess_row_size)
        mpi.world.barrier()
        self.solver.max_iter = self.max_iterations

    #
    def Solve(self):
        mpi.world.barrier()
        (self.solver).Solve()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
