import KratosMultiphysics

class ResidualBasedNewtonRaphsonStrategyPython():


    def __init__(self, main_model_part, scheme, convergence_criterion, builder_and_solver, max_iterations, compute_reactions, reform_dofs_at_each_step, move_mesh_flag):
        self.main_model_part = main_model_part
        self.scheme = scheme
        self.convergence_criteria = convergence_criterion
        self.builder_and_solver = builder_and_solver
        self.max_iterations = max_iterations
        self.CalculateReactionsFlag = compute_reactions
        self.ReformDofSetAtEachStep = reform_dofs_at_each_step
        self.MoveMeshFlag = move_mesh_flag

        self.builder_and_solver.SetCalculateReactionsFlag(self.CalculateReactionsFlag)
        self.builder_and_solver.SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)

        self.rebuild_level = 2
        self.echo_level = 1

        self.space_utils = KratosMultiphysics.UblasSparseSpace()

        self.pA = self.space_utils.CreateEmptyMatrixPointer()
        self.pDx = self.space_utils.CreateEmptyVectorPointer()
        self.pb = self.space_utils.CreateEmptyVectorPointer()

        self.A = (self.pA)
        self.Dx = (self.pDx)
        self.b = (self.pb)

        #initialize flags
        self.SolutionStepIsInitialized = False
        self.InitializeWasPerformed = False
        self.StiffnessMatrixIsBuilt = False
        self.UseOldStiffnessInFirstIteration = False

    def Initialize(self):
        if (self.InitializeWasPerformed == False):
            if(self.scheme.SchemeIsInitialized() == False):
                self.scheme.Initialize(self.main_model_part)

            if (self.scheme.ElementsAreInitialized() == False):
                self.scheme.InitializeElements(self.main_model_part)

            if (self.scheme.ConditionsAreInitialized() == False):
                self.scheme.InitializeConditions(self.main_model_part)

            #if (self.convergence_criteria.IsInitialized() == False):
            self.convergence_criteria.Initialize(self.main_model_part)

            self.InitializeWasPerformed == True

    def Predict(self):
        self.scheme.Predict(self.main_model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

    def InitializeSolutionStep(self):
        if (self.builder_and_solver.GetDofSetIsInitializedFlag() == False or self.ReformDofSetAtEachStep == True):
            #initialize the list of degrees of freedom to be used
            self.builder_and_solver.SetUpDofSet(self.scheme,self.main_model_part)
            #reorder the list of degrees of freedom to identify fixity and system size
            self.builder_and_solver.SetUpSystem(self.main_model_part)
            #allocate memory for the system and preallocate the structure of the matrix
            self.builder_and_solver.ResizeAndInitializeVectors(self.scheme, self.pA,self.pDx,self.pb,self.main_model_part)

            self.A = (self.pA)
            self.Dx = (self.pDx)
            self.b = (self.pb)

        self.builder_and_solver.InitializeSolutionStep(self.main_model_part,self.A,self.Dx,self.b)
        self.scheme.InitializeSolutionStep(self.main_model_part,self.A,self.Dx,self.b)

    def SolveSolutionStep(self, projection_module):
        self.projection_module = projection_module

        max_iterations = self.max_iterations
        dof_set = self.builder_and_solver.GetDofSet()
        iteration_number = 1

        self.main_model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] = iteration_number

        self.scheme.InitializeNonLinIteration(self.main_model_part,self.A,self.Dx,self.b)
        self.convergence_criteria.InitializeNonLinearIteration(self.main_model_part,dof_set,self.A,self.Dx,self.b)
        is_converged = self.convergence_criteria.PreCriteria(self.main_model_part,dof_set,self.A,self.Dx,self.b)

        if (self.rebuild_level > 0 or self.StiffnessMatrixIsBuilt == False):
            self.space_utils.SetToZeroVector(self.Dx)
            self.space_utils.SetToZeroVector(self.b)
            self.space_utils.SetToZeroMatrix(self.A)
            if (self.UseOldStiffnessInFirstIteration == True):
                self.builder_and_solver.BuildAndSolveLinearizedOnPreviousIteration(self.scheme,self.main_model_part,self.A,self.Dx,self.b, self.MoveMeshFlag)
            else:
                self.builder_and_solver.BuildAndSolve(self.scheme,self.main_model_part,self.A,self.Dx,self.b)
        else:
            self.space_utils.SetToZeroVector(self.Dx)
            self.space_utils.SetToZeroVector(self.b)

            self.builder_and_solver.BuildRHSAndSolve(self.scheme,self.main_model_part,self.A,self.Dx,self.b)

        self.scheme.Update(self.main_model_part,dof_set,self.A,self.Dx,self.b)

        self.scheme.FinalizeNonLinIteration(self.main_model_part,self.A,self.Dx,self.b)
        self.convergence_criteria.FinalizeNonLinearIteration(self.main_model_part,dof_set,self.A,self.Dx,self.b)

        # self.projection_module.ApplyForwardCoupling(1)
        # self.projection_module.UpdateHydrodynamicForces()
        # self.projection_module.ProjectFromParticles(False)

        if (is_converged):
            if (self.convergence_criteria.GetActualizeRHSflag()):
                self.space_utils.SetToZeroVector(self.b)
                self.builder_and_solver.BuildRHS(self.scheme,self.main_model_part,self.b)

            is_converged = self.convergence_criteria.PostCriteria(self.main_model_part,dof_set,self.A,self.Dx,self.b)

        while (is_converged == False and iteration_number < max_iterations):

            iteration_number += 1

            self.main_model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] = iteration_number

            self.scheme.InitializeNonLinIteration(self.main_model_part,self.A,self.Dx,self.b)
            self.convergence_criteria.InitializeNonLinearIteration(self.main_model_part,dof_set,self.A,self.Dx,self.b)

            is_converged = self.convergence_criteria.PreCriteria(self.main_model_part,dof_set,self.A,self.Dx,self.b)

            if (self.space_utils.Size(self.Dx) != 0):
                if (self.rebuild_level > 1 or self.StiffnessMatrixIsBuilt == False):
                    self.space_utils.SetToZeroVector(self.Dx)
                    self.space_utils.SetToZeroVector(self.b)
                    self.space_utils.SetToZeroMatrix(self.A)
                    self.builder_and_solver.BuildAndSolve(self.scheme,self.main_model_part,self.A,self.Dx,self.b)
                else:
                    self.space_utils.SetToZeroVector(self.Dx)
                    self.space_utils.SetToZeroVector(self.b)
                    self.builder_and_solver.BuildRHSAndSolve(self.scheme,self.main_model_part,self.A,self.Dx,self.b)

            self.scheme.Update(self.main_model_part,dof_set,self.A,self.Dx,self.b)

            self.scheme.FinalizeNonLinIteration(self.main_model_part,self.A,self.Dx,self.b)
            self.convergence_criteria.FinalizeNonLinearIteration(self.main_model_part,dof_set,self.A,self.Dx,self.b)

            # self.projection_module.ApplyForwardCoupling(1.0)
            # self.projection_module.UpdateHydrodynamicForces()
            # self.projection_module.ProjectFromParticles(False)

            residual_is_updated = False

            if (is_converged):
                if (self.convergence_criteria.GetActualizeRHSflag()):
                    self.space_utils.SetToZeroVector(self.b)
                    self.builder_and_solver.BuildRHS(self.scheme,self.main_model_part,self.b)
                    residual_is_updated = True
                is_converged = self.convergence_criteria.PostCriteria(self.main_model_part,dof_set,self.A,self.Dx,self.b)

        if (iteration_number >= max_iterations):
            print('NOT CONVERGENCE!!!!!')
        else:
            print("Convergence achieved after " + str(iteration_number) + " / " + str(max_iterations) + " iterations")

        return is_converged

    def Check(self):
        self.main_model_part.Check()
        self.builder_and_solver.Check(self.main_model_part)
        self.scheme.Check(self.main_model_part)
        self.convergence_criteria.Check(self.main_model_part)

    def FinalizeSolutionStep(self):
        dof_set = self.builder_and_solver.GetDofSet()
        self.scheme.FinalizeSolutionStep(self.main_model_part,self.A,self.Dx,self.b)
        self.builder_and_solver.FinalizeSolutionStep(self.main_model_part,self.A,self.Dx,self.b)
        self.convergence_criteria.FinalizeSolutionStep(self.main_model_part,dof_set,self.A,self.Dx,self.b)
        self.scheme.Clean()
        if (self.ReformDofSetAtEachStep == True):
            self.Clear()

    def Clear(self):
        self.builder_and_solver.SetDofSetIsInitializedFlag(False)
        self.builder_and_solver.Clear()
        self.space_utils.Clear(self.Dx)
        self.space_utils.Clear(self.b)
        self.space_utils.Clear(self.A)
        self.scheme.Clear()
        self.InitializeWasPerformed = False

    def SetEchoLevel(self, echo_level):
        self.builder_and_solver.SetEchoLevel(echo_level)