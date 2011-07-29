#importing the Kratos Library
from Kratos import *
from KratosTrilinosApplication import *
import mpi
from trilinos_strategy_python import SolvingStrategyPython # Base class

class DecoupledUPStrategyPython(SolvingStrategyPython):

    def __init__(self,builder_and_solver_type,model_part,time_scheme,\
                 vel_linear_solver,press_linear_solver,\
                 convergence_criteria,max_iter,CalculateReactionsFlag,\
                 ReformDofSetAtEachStep,MoveMeshFlag,Comm,guess_row_size):

        #save the input parameters
        self.model_part = model_part
        self.scheme = time_scheme
        self.vel_linear_solver = vel_linear_solver
        self.press_linear_solver = press_linear_solver
        self.convergence_criteria = convergence_criteria
        self.CalculateReactionsFlag = CalculateReactionsFlag
        self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
        self.MoveMeshFlag = MoveMeshFlag

        self.space_utils = TrilinosSparseSpace()

        #default values for some variables
        self.max_iter = max_iter
        self.echo_level = 1
##        self.builder_and_solver_type = builder_and_solver_type
	if(builder_and_solver_type == "PressureSplitting"):
		VelocityCorrection = 2
		UseInexactNewton = False
		IN_MinTol = 0.001
		IN_MaxTol = 0.1
		IN_Gamma = 0.9
        	self.builder_and_solver = \
                    TrilinosPressureSplittingBuilderAndSolver(\
                        Comm,guess_row_size,self.vel_linear_solver,\
                        self.press_linear_solver,\
                        VelocityCorrection,UseInexactNewton,\
                        IN_MinTol,IN_MaxTol,IN_Gamma)
        else:
            raise "Unknown Builder And Solver Type"
        	
        self.SetEchoLevel(self.echo_level)

  

        #local matrices and vectors
        self.pA = self.space_utils.CreateEmptyMatrixPointer(Comm)
        self.pDx = self.space_utils.CreateEmptyVectorPointer(Comm)
        self.pb = self.space_utils.CreateEmptyVectorPointer(Comm)

        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()

        #initialize flags
        self.SolutionStepIsInitialized = False
        self.InitializeWasPerformed = False
        self.StiffnessMatrixIsBuilt = False

        #provide settings to the builder and solver
        (self.builder_and_solver).SetCalculateReactionsFlag(\
            self.CalculateReactionsFlag);
        (self.builder_and_solver).SetReshapeMatrixFlag(\
            self.ReformDofSetAtEachStep);
