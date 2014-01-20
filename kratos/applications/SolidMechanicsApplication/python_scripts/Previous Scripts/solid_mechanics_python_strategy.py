from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *

CheckForPreviousImport()


class ResidualStrategy:
    #

    def __init__(self, model_part, time_scheme, linear_solver, convergence_criteria, CalculateReactionsFlag, ReformDofSetAtEachStep, LineSearchActivation):
        # save the input parameters
        self.space_utils = UblasSparseSpace()
        self.model_part = model_part
        self.scheme = time_scheme
        self.linear_solver = linear_solver
        self.convergence_criteria = convergence_criteria
        self.CalculateReactionsFlag = CalculateReactionsFlag
        self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
        self.EquationSystemSize = 0

        # default values for some variables
        # for normal execution
        self.builder_and_solver = ResidualBasedBuilderAndSolver(self.linear_solver)
        # to conserve matrix blocks AMG solving
        # self.builder_and_solver = BlockResidualBasedBuilderAndSolver(self.linear_solver)

        self.max_iter = 15
        self.echo_level = 1

        # local matrices and vectors
        self.pA = self.space_utils.CreateEmptyMatrixPointer()
        self.pDx = self.space_utils.CreateEmptyVectorPointer()
        self.pb = self.space_utils.CreateEmptyVectorPointer()

        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()

        # initialize flags
        self.SolutionStepIsInitialized = False
        self.InitializeWasPerformed = False

        # provide settings to the builder and solver
        (self.builder_and_solver).SetCalculateReactionsFlag(self.CalculateReactionsFlag)
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)

    #
    def Initialize(self):
        if(self.scheme.SchemeIsInitialized() == False):
            self.scheme.Initialize(self.model_part)

        if(self.scheme.ElementsAreInitialized() == False):
            self.scheme.InitializeElements(self.model_part)

        if(self.scheme.ConditionsAreInitialized() == False):
            self.scheme.InitializeConditions(self.model_part)

    #
    #
    def Solve(self):
        # perform the operations to be performed ONCE and ensure they will not be repeated
        # elemental function "Initialize" is called here
        if(self.InitializeWasPerformed == False):
            print(" INITIALIZE ")
            self.Initialize()
            self.InitializeWasPerformed = True

        # perform initializations for the current step
        # this operation implies:
        # identifying the set of DOFs that will be solved during this step
        # organizing the DOFs so to identify the dirichlet conditions
        # resizing the matrix preallocating the "structure"
        if(self.SolutionStepIsInitialized == False):
            self.InitializeSolutionStep()
            self.ReformStepDofs()

        # initialize iteration number
        self.model_part.ProcessInfo[NL_ITERATION_NUMBER] = 1

        # perform prediction
        print(" Update ")
        self.Predict()

        # execute iteration - first iteration is ALWAYS executed
        calculate_norm = False
        normDx = self.ExecuteIteration(self.echo_level, calculate_norm)
        it = 1

        print(" ")
        print("Initial iteration DONE")
        print(" ")

        # Start the loop of solution
        print(" ## --START THE NON LINEAR LOOP OF SOLUTION-- ## ")
        print(" ")

        converged = False
        while(it < self.max_iter and converged == False):

            # verify convergence
            print(" ")
            print("Iteration Number =", it)

            # Call Initialize iteration on Elements and Conditions
            self.scheme.InitializeNonLinIteration(self.model_part, self.A, self.Dx, self.b)
            # self.model_part.ProcessInfo[NL_ITERATION_NUMBER] = it
            converged = self.convergence_criteria.PreCriteria(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)

            # print " [ ELEMENTS: ",self.model_part.NumberOfElements(), " CONDICIONS :",self.model_part.NumberOfConditions()," ]"

            # calculate iteration
            # - system is built and solved
            # - database is updated depending on the solution
            # - nodal coordinates are updated if required
            normDx = self.ExecuteIteration(self.echo_level, calculate_norm)

            # verify convergence
            print(" ## --Convergence-- ## ")
            converged = self.convergence_criteria.PostCriteria(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)

            # perform update
            # Displacement, Velocity and Acceleration Update
            print(" Update ")
            self.scheme.Update(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)

            # Call Finalize iteration on Elements and Conditions
            self.scheme.FinalizeNonLinIteration(self.model_part, self.A, self.Dx, self.b)

            # update iteration count
            it = it + 1

            if(it >= self.max_iter):
                print(" ******************************************** ")
                print(" Warning: No Solution For This Step Was Found ")
                print(" ******************************************** ")

        # compute reactions
        if(self.CalculateReactionsFlag):
            self.ComputeReactions()

        # finalize the solution step
        self.FinalizeSolutionStep()

        # clean step
        self.CleanStep()

        print(" Step Finished ")

        # clear if needed - deallocates memory
        if(self.ReformDofSetAtEachStep):
            self.ClearMemory()

    #
    #
    #
    def Predict(self):
        self.scheme.Predict(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)

    #
    def InitializeSolutionStep(self):
        # Initialize of the solution step
        print(" INITIALIZE SOLUTION STEP ")
        self.builder_and_solver.InitializeSolutionStep(self.model_part, self.A, self.Dx, self.b)
        self.scheme.InitializeSolutionStep(self.model_part, self.A, self.Dx, self.b)

        self.SolutionStepIsInitialized = True

    #
    def ReformStepDofs(self):

        if(self.builder_and_solver.GetDofSetIsInitializedFlag() == False or self.ReformDofSetAtEachStep == True):
            print(" REFORM DOF ")
            # initialize the list of degrees of freedom to be used
            self.builder_and_solver.SetUpDofSet(self.scheme, self.model_part)
            # reorder the list of degrees of freedom to identify fixity and system size
            self.builder_and_solver.SetUpSystem(self.model_part)
            # allocate memory for the system and preallocate the structure of the matrix
            self.builder_and_solver.ResizeAndInitializeVectors(self.pA, self.pDx, self.pb, self.model_part.Elements, self.model_part.Conditions, self.model_part.ProcessInfo)

            self.EquationSystemSize = self.builder_and_solver.GetEquationSystemSize()

            # updating references
            self.A = (self.pA).GetReference()
            self.Dx = (self.pDx).GetReference()
            self.b = (self.pb).GetReference()

    #
    def FinalizeSolutionStep(self):
        # Finalization of the solution step
        print(" FINALIZE SOLUTION STEP ")
        self.builder_and_solver.FinalizeSolutionStep(self.model_part, self.A, self.Dx, self.b)
        self.scheme.FinalizeSolutionStep(self.model_part, self.A, self.Dx, self.b)

        self.scheme.Clean()
        # reset flags for the next step
        self.SolutionStepIsInitialized = False

    #
    def ExecuteIteration(self, echo_level, CalculateNormDxFlag):

        # reset system matrices and vectors prior to rebuild
        self.space_utils.SetToZeroMatrix(self.A)
        self.space_utils.SetToZeroVector(self.Dx)
        self.space_utils.SetToZeroVector(self.b)

        print(" Start iteration ")
        print(" ---->>><<<<---- ")

        # build and solve the problem
        self.builder_and_solver.BuildAndSolve(self.scheme, self.model_part, self.A, self.Dx, self.b)

        # full output if needed
        if(echo_level == 2):
            print(" ")
            print(" //---GLOBAL LINEAR SYSTEM SOLUTION---// ")
            # print " "
            # print "LHS = ", self.A
            print(" ")
            print("RHS = ", self.b)
            print(" ")
            print("Solution = ", self.Dx)
            print(" ")

        if(echo_level == 3):
            print(" ")
            print(" //---GLOBAL LINEAR SYSTEM SOLUTION---// ")
            print(" ")
            print("RHS = ", self.b)
            print(" ")
            print("Solution = ", self.Dx)
            print(" ")
            for node in self.model_part.Nodes:
                print("Displacement = ", node.GetSolutionStepValue(DISPLACEMENT, 0))
                print("Pressure = ", node.GetSolutionStepValue(PRESSURE, 0))
                print("Node = (", node.X, " , ", node.Y, " , ", node.Z, " )")

        print(" End iteration ")

        # calculate the norm of the "correction" Dx (not in use)
        if(CalculateNormDxFlag):
            normDx = self.space_utils.TwoNorm(self.Dx)
        else:
            normDx = 0.0

        return normDx

    #
    def ComputeReactions(self):
        print("***********************************************************")
        if(self.builder_and_solver.GetDofSetIsInitializedFlag() == False):
            self.ReformStepDofs()
        print(" -- CALCULATE REACTIONS START -- ")
        self.builder_and_solver.CalculateReactions(self.scheme, self.model_part, self.A, self.Dx, self.b)
        print(" -- CALCULATE REACTIONS END -- ")

    #
    def CleanStep(self):
        self.scheme.Clean()
        # reset flags for the next step
        self.SolutionStepIsInitialized = False

    #
    def ClearMemory(self):
        self.space_utils.ClearMatrix(self.pA)
        self.space_utils.ResizeMatrix(self.A, 0, 0)

        self.space_utils.ClearVector(self.pDx)
        self.space_utils.ResizeVector(self.Dx, 0)

        self.space_utils.ClearVector(self.pb)
        self.space_utils.ResizeVector(self.b, 0)

        # updating references
        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()

        self.builder_and_solver.SetDofSetIsInitializedFlag(False)
        self.builder_and_solver.Clear()

    #
    def SetEchoLevel(self, level):
        self.echo_level = level
        self.builder_and_solver.SetEchoLevel(level)

    #
    def SetInitialized(self):
        self.InitializeWasPerformed = True

    #
    def SetMaxIters(self, iters):
        self.max_iter = iters

    #
    def Check(self):
        self.builder_and_solver.Check(self.model_part)
        self.scheme.Check(self.model_part)
        self.convergence_criteria.Check(self.model_part)
