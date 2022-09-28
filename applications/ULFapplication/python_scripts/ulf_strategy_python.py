 #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ULFApplication import *
# import cProfile
import time


class ULFStrategyPython:
    #

    def __init__(self, model_part, time_scheme, linear_solver, convergence_criteria, CalculateReactionsFlag, ReformDofSetAtEachStep, MoveMeshFlag, domain_size):
        # save the input parameters
        self.model_part = model_part
        self.scheme = time_scheme
        self.linear_solver = linear_solver
        self.convergence_criteria = convergence_criteria
        self.CalculateReactionsFlag = CalculateReactionsFlag
        self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
        self.MoveMeshFlag = MoveMeshFlag
        self.domain_size = domain_size

        self.space_utils = UblasSparseSpace()

        # temporary ... i need it to calculate the nodal area
        self.UlfUtils = UlfUtils()

        # default values for some variables
        self.max_iter = 30
        self.rebuild_level = 1  # rebuild at each solution step
        self.echo_level = 1
        self.builder_and_solver = ResidualBasedEliminationBuilderAndSolver(self.linear_solver)
        self.builder_and_solver.SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)

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
        self.StiffnessMatrixIsBuilt = False

        # provide settings to the builder and solver
        (self.builder_and_solver).SetCalculateReactionsFlag(self.CalculateReactionsFlag)
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)

        (self.VariableUtils) = VariableUtils()

    #
    def SetEchoLevel(self, level):
        self.echo_level = level

    #
    def Initialize(self):
        if(self.scheme.SchemeIsInitialized() == False):
            self.scheme.Initialize(self.model_part)

        if (self.scheme.ElementsAreInitialized() == False):
            self.scheme.InitializeElements(self.model_part)

    #
    #
    def Solve(self, domain_size, UlfUtils):
        # perform the operations to be performed ONCE and ensure they will not be repeated
        # elemental function "Initialize" is called here
        if(self.InitializeWasPerformed == False):
            self.Initialize()
            self.InitializeWasPerformed = True

        print("in the solve")
        print(self.model_part)

        # perform initializations for the current step
        # this operation implies:
        # identifying the set of DOFs that will be solved during this step
        # organizing the DOFs so to identify the dirichlet conditions
        # resizing the matrix preallocating the "structure"
        import time
        step_initialization_start = time.clock()
        reform_dofs = True
        self.InitializeSolutionStep(reform_dofs)
        print("step initizalization time =", time.clock() - step_initialization_start)

        # perform prediction
        self.Predict()

        if(self.MoveMeshFlag):
            self.scheme.MoveMesh(self.model_part.Nodes)

        # check for inverted elements
        inverted_elements = False
        print("volume calculation")
        volume = (UlfUtils).CalculateVolume(self.model_part, self.domain_size)
        print("finished volume calculation")
        if(volume <= 0.00):
            inverted_elements = True
            print("INVERTED ELEMENT FOUND - right after prediction")

        print("before iteration")
        # execute iteration - first iteration is ALWAYS executed
        calculate_norm = False
        if(inverted_elements == False):
            normDx = self.ExecuteIteration(self.echo_level, self.MoveMeshFlag, calculate_norm)
        it = 1

        print("after iteration")
# print self.A
# print self.b

        # check for inverted elements
        volume = (UlfUtils).CalculateVolume(self.model_part, self.domain_size)
        if(volume <= 0.00):
            inverted_elements = True
            print("INVERTED ELEMENT FOUND - just after first iteration")

        # non linear loop
        print("non linear loop")
        converged = False
        while(it < self.max_iter and converged == False and inverted_elements == False):
            # verify convergence
            converged = self.convergence_criteria.PreCriteria(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)

            # calculate iteration
            # - system is built and solved
            # - database is updated depending on the solution
            # - nodal coordinates are updated if required
            normDx = self.ExecuteIteration(self.echo_level, self.MoveMeshFlag, calculate_norm)

            # verify convergence
            converged = self.convergence_criteria.PostCriteria(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)

            # check for inverted elements
            volume = (UlfUtils).CalculateVolume(self.model_part, self.domain_size)
            if(volume <= 0.00):
                inverted_elements = True
                print("INVERTED ELEMENT FOUND")

            # update iteration count
            it = it + 1

        # moving lonely nodes
        (self.UlfUtils).MoveLonelyNodes(self.model_part)

        # finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)

        # clear if needed - deallocates memory
        self.Clear()

        return inverted_elements

    #
    #
    #
    def Predict(self):
        self.scheme.Predict(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)

    #
    def InitializeSolutionStep(self, reform_dofs):
        if(reform_dofs):
            print("reforming dofs")
            print("Inside the Initialize Solution Step --> ", self.model_part)
            # initialize the list of degrees of freedom to be used
            self.builder_and_solver.SetUpDofSet(self.scheme, self.model_part);

            # reorder the list of degrees of freedom to identify fixity and system size
            self.builder_and_solver.SetUpSystem(self.model_part)

            # allocate memory for the system and preallocate the structure of the matrix
            self.builder_and_solver.ResizeAndInitializeVectors(self.scheme, self.pA, self.pDx, self.pb, self.model_part);

            # updating references
            self.A = (self.pA).GetReference()
            self.Dx = (self.pDx).GetReference()
            self.b = (self.pb).GetReference()

# print "A size1 ",self.space_utils.Size1(self.A)
# print "b size ",self.space_utils.Size(self.b)
# print "Dx size ",self.space_utils.Size(self.Dx)
        self.builder_and_solver.InitializeSolutionStep(self.model_part, self.A, self.Dx, self.b)
        self.scheme.InitializeSolutionStep(self.model_part, self.A, self.Dx, self.b)

    #
    def ExecuteIteration(self, echo_level, MoveMeshFlag, CalculateNormDxFlag):
# print "A size ",len(self.A)
# print "b size ",self.space_utils.Size(self.b)
# print "Dx size ",self.space_utils.Size(self.Dx)
        # reset system matrices and vectors prior to rebuild
        self.space_utils.SetToZeroMatrix(self.A)
        self.space_utils.SetToZeroVector(self.Dx)
        self.space_utils.SetToZeroVector(self.b)

# print "ln188"
# print "b=",self.b
        print("size of problem ", self.space_utils.Size(self.b))

        self.scheme.InitializeNonLinIteration(self.model_part, self.A, self.Dx, self.b)
        print("ln191")
        # build and solve the problem
        self.builder_and_solver.BuildAndSolve(self.scheme, self.model_part, self.A, self.Dx, self.b)

        # full output if needed
        if(echo_level >= 3):
            print("SystemMatrix = ", self.A)
            print("solution obtained = ", self.Dx)
            print("RHS = ", self.b)
        print("ln200")
        # perform update
        self.scheme.Update(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b);

        # move the mesh as needed
        if(MoveMeshFlag):
            self.scheme.MoveMesh(self.model_part.Nodes);

        self.scheme.FinalizeNonLinIteration(self.model_part, self.A, self.Dx, self.b)

        # calculate the norm of the "correction" Dx
        if(CalculateNormDxFlag):
            normDx = self.space_utils.TwoNorm(self.Dx)
        else:
            normDx = 0.0

        return normDx

    #
    def FinalizeSolutionStep(self, CalculateReactionsFlag):
        if(CalculateReactionsFlag):
            # in the other solvers the reactions are reset within the function CalculateReactions of the quasi inc builder and solver
            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(REACTION_X, 0, 0.0)
                node.SetSolutionStepValue(REACTION_Y, 0, 0.0)
                node.SetSolutionStepValue(REACTION_Z, 0, 0.0)

            self.builder_and_solver.CalculateReactions(self.scheme, self.model_part, self.A, self.Dx, self.b)

        # Finalisation of the solution step,
        self.scheme.FinalizeSolutionStep(self.model_part, self.A, self.Dx, self.b)
        self.builder_and_solver.FinalizeSolutionStep(self.model_part, self.A, self.Dx, self.b)

        self.scheme.Clean()

        # reset flags for the next step
        self.mSolutionStepIsInitialized = False

    #
    def PredictionStep(self, domain_size, UlfUtils):
        reform_dofs = True
        self.InitializeSolutionStep(reform_dofs);

        # perform prediction
        self.Predict()

        if(self.MoveMeshFlag):
            self.scheme.MoveMesh(self.model_part.Nodes);

    #
    def MoveMesh(self):
        self.scheme.MoveMesh(self.model_part.Nodes);

    #
    def Clear(self):
        # self.space_utils.ClearMatrix(self.A)
        self.space_utils.ResizeMatrix(self.A, 0, 0)

        # self.space_utils.ClearVector(self.Dx)
        self.space_utils.ResizeVector(self.Dx, 0)

        # self.space_utils.ClearVector(self.b)
        self.space_utils.ResizeVector(self.b, 0)

        self.builder_and_solver.SetDofSetIsInitializedFlag(False)

        self.builder_and_solver.Clear()
        print("clear completed")
