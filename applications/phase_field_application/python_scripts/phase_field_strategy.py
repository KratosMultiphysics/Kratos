from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import math

from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.PhaseFieldApplication import *
from KratosMultiphysics.mpi import *
CheckForPreviousImport()

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(INTERNAL_FORCE)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_DT)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL_DT)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS_DT)
    model_part.AddNodalSolutionStepVariable(ACCELERATION_NULL)
    model_part.AddNodalSolutionStepVariable(ACCELERATION_EINS)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(ELASTIC_LEFT_CAUCHY_GREEN_OLD)
    model_part.AddNodalSolutionStepVariable(INSITU_STRESS)
    model_part.AddNodalSolutionStepVariable(PRESTRESS)
    model_part.AddNodalSolutionStepVariable(STRESSES)
    model_part.AddNodalSolutionStepVariable(STRAIN)
    model_part.AddNodalSolutionStepVariable(FACE_LOAD)
    model_part.AddNodalSolutionStepVariable(ROTATION)
    model_part.AddNodalSolutionStepVariable(PHASE_FIELD)
    model_part.AddNodalSolutionStepVariable(PHASE_FIELD_DUAL_VARIABLE)
    model_part.AddNodalSolutionStepVariable(PRESCRIBED_DELTA_DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
    print("variables for the phase field problem added correctly")

def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(DISPLACEMENT_X, REACTION_X)
        node.AddDof(DISPLACEMENT_Y, REACTION_Y)
        node.AddDof(DISPLACEMENT_Z, REACTION_Y)
        node.AddDof(PHASE_FIELD, PHASE_FIELD_DUAL_VARIABLE)
    print("dofs for the phase field problem added correctly")

class SampleSolver():

    def __init__(self, model_part, abs_tol, rel_tol):
        self.echo_level = 0
        self.toll = rel_tol
        self.absolute_tol = abs_tol
        self.structure_linear_solver = SkylineLUFactorizationSolver()
        self.conv_criteria = DisplacementCriteria(rel_tol, abs_tol)
        self.CalculateReactionFlag = False

    def Initialize(self):
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
#        self.conv_criteria = MultiPhaseFlowCriteria(self.toll, self.absolute_tol)
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag = False
        self.MaxNewtonRapshonIterations = 30
        self.parallel_space = UblasSparseSpace()
        self.solver = PhaseFieldStaggeredSolver(self.parallel_space, None, self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.builder_and_solver, self.MaxNewtonRapshonIterations, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag)

class PhaseFieldStaggeredSolver:
    def __init__(self, parallel_space, comm, model_part, time_scheme, linear_solver, conv_criteria, builder_and_solver, MaxNewtonRapshonIterations, CalculateReactionsFlag, ReformDofSetAtEachStep, MoveMeshFlag, params):
        self.comm = comm
        self.parallel_space = parallel_space
        self.model_part = model_part
        self.time_scheme = time_scheme
        self.linear_solver = linear_solver
        self.conv_criteria = conv_criteria
        self.builder_and_solver = builder_and_solver
        self.max_iter = MaxNewtonRapshonIterations
        self.CalculateReactionsFlag = CalculateReactionsFlag
        self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
        self.MoveMeshFlag = MoveMeshFlag
        self.params = params

        if self.params == None: # default values 
            self.echo_level = 1
            self.toll = 1.0e-6
        else:
            self.echo_level = self.params['Echo Level']
            self.toll = self.params['Tolerance']

        # local matrices and vectors
        if(self.comm != None):
            self.pA = self.parallel_space.CreateEmptyMatrixPointer(self.comm)
            self.pDx = self.parallel_space.CreateEmptyVectorPointer(self.comm)
            self.pb = self.parallel_space.CreateEmptyVectorPointer(self.comm)
        else:
            self.pA = self.parallel_space.CreateEmptyMatrixPointer()
            self.pDx = self.parallel_space.CreateEmptyVectorPointer()
            self.pb = self.parallel_space.CreateEmptyVectorPointer()

        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()

        # provide settings to the builder and solver
        (self.builder_and_solver).SetCalculateReactionsFlag(self.CalculateReactionsFlag)
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)

        self.iterations_last_solution = 0

    #
    def Initialize(self):
        if(self.time_scheme.SchemeIsInitialized() == False):
            self.time_scheme.Initialize(self.model_part)

        if (self.time_scheme.ElementsAreInitialized() == False):
            self.time_scheme.InitializeElements(self.model_part)

        if (self.time_scheme.ConditionsAreInitialized() == False):
            self.time_scheme.InitializeConditions(self.model_part)

    def Solve(self):
        return self.SolveStaggeredOneStep()
#        self.SolveStaggeredMultiStep()

    def SolveStaggeredOneStep(self):
        # freeze the PHASE_FIELD dofs
        free_node_list = []
        for node in self.model_part.Nodes:
            if not node.IsFixed(PHASE_FIELD):
                free_node_list.append(node)
                node.Fix(PHASE_FIELD)

        # initialize flags
        self.SolutionStepIsInitialized = False
        self.InitializeWasPerformed = False
        self.StiffnessMatrixIsBuilt = False

        # set the block size (useful for AMG preconditioner)
        if("BlockBuilderAndSolver" in self.builder_and_solver.__class__.__name__):
            self.builder_and_solver.SetCustomBlockSize(3)

        # solve for displacement
        self.model_part.ProcessInfo[FRACTIONAL_STEP] = 0
        converged = self.SolveOneStep()
        for node in self.model_part.Nodes:
            if node.IsFixed(DISPLACEMENT_X) and node.SolutionStepsDataHas(PRESCRIBED_DELTA_DISPLACEMENT_X):
                ux = node.GetSolutionStepValue(DISPLACEMENT_X)
                dux = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X)
                node.SetSolutionStepValue(DISPLACEMENT_X, ux + dux)
            if node.IsFixed(DISPLACEMENT_Y) and node.SolutionStepsDataHas(PRESCRIBED_DELTA_DISPLACEMENT_Y):
                uy = node.GetSolutionStepValue(DISPLACEMENT_Y)
                duy = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y)
                node.SetSolutionStepValue(DISPLACEMENT_Y, uy + duy)
            if node.IsFixed(DISPLACEMENT_Z) and node.SolutionStepsDataHas(PRESCRIBED_DELTA_DISPLACEMENT_Z):
                uz = node.GetSolutionStepValue(DISPLACEMENT_Z)
                duz = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z)
                node.SetSolutionStepValue(DISPLACEMENT_Z, uz + duz)
        print("Solve for displacements completed")

        # unfreeze the PHASE_FIELD dofs
        for node in free_node_list:
            node.Free(PHASE_FIELD)

        # if we did not get convergence for displacement, then we stop
        if not converged:
            return converged

        ################################################

        # freeze the DISPLACEMENT dofs
        free_node_list_x = []
        free_node_list_y = []
        free_node_list_z = []
        for node in self.model_part.Nodes:
            if not node.IsFixed(DISPLACEMENT_X):
                free_node_list_x.append(node)
                node.Fix(DISPLACEMENT_X)
            if not node.IsFixed(DISPLACEMENT_Y):
                free_node_list_y.append(node)
                node.Fix(DISPLACEMENT_Y)
            if not node.IsFixed(DISPLACEMENT_Z):
                free_node_list_z.append(node)
                node.Fix(DISPLACEMENT_Z)

        # initialize flags
        self.SolutionStepIsInitialized = False
        self.InitializeWasPerformed = False
        self.StiffnessMatrixIsBuilt = False

        # set the block size (useful for AMG preconditioner)
        if("BlockBuilderAndSolver" in self.builder_and_solver.__class__.__name__):
            self.builder_and_solver.SetCustomBlockSize(1)

        # solve for phase field
        self.model_part.ProcessInfo[FRACTIONAL_STEP] = 1
        converged = self.SolveOneStep()
        print("Solve for phase field completed")

        # unfreeze the displacement dofs
        for node in free_node_list_x:
            node.Free(DISPLACEMENT_X)
        for node in free_node_list_y:
            node.Free(DISPLACEMENT_Y)
        for node in free_node_list_z:
            node.Free(DISPLACEMENT_Z)

        # if we did not get the convergence for phase field, then we stop
        if not converged:
            return converged

    def SolveStaggeredMultiStep(self):
        it = 0
        converged = False
        staggered_max_iter = self['Staggered Max Iteration Number'] # 30
        while (it < staggered_max_iter) and (converged == False):
            # freeze the PHASE_FIELD dofs
            free_node_list = []
            for node in self.model_part.Nodes:
                if not node.IsFixed(PHASE_FIELD):
                    free_node_list.append(node)
                    node.Fix(PHASE_FIELD)

            # solve for displacement
            self.model_part.ProcessInfo[FRACTIONAL_STEP] = 0
            self.SolveOneStep()
            for node in self.model_part.Nodes:
                if node.IsFixed(DISPLACEMENT_X) and node.SolutionStepsDataHas(PRESCRIBED_DELTA_DISPLACEMENT_X):
                    ux = node.GetSolutionStepValue(DISPLACEMENT_X)
                    dux = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X)
                    node.SetSolutionStepValue(DISPLACEMENT_X, ux + dux)
                if node.IsFixed(DISPLACEMENT_Y) and node.SolutionStepsDataHas(PRESCRIBED_DELTA_DISPLACEMENT_Y):
                    uy = node.GetSolutionStepValue(DISPLACEMENT_Y)
                    duy = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y)
                    node.SetSolutionStepValue(DISPLACEMENT_Y, uy + duy)
                if node.IsFixed(DISPLACEMENT_Z) and node.SolutionStepsDataHas(PRESCRIBED_DELTA_DISPLACEMENT_Z):
                    uz = node.GetSolutionStepValue(DISPLACEMENT_Z)
                    duz = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z)
                    node.SetSolutionStepValue(DISPLACEMENT_Z, uz + duz)
            print("Solve for displacements completed")

            # unfreeze the PHASE_FIELD dofs
            for node in free_node_list:
                node.Free(PHASE_FIELD)

            ################################################

            # freeze the DISPLACEMENT dofs
            free_node_list_x = []
            free_node_list_y = []
            free_node_list_z = []
            for node in self.model_part.Nodes:
                if not node.IsFixed(DISPLACEMENT_X):
                    free_node_list_x.append(node)
                    node.Fix(DISPLACEMENT_X)
                if not node.IsFixed(DISPLACEMENT_Y):
                    free_node_list_y.append(node)
                    node.Fix(DISPLACEMENT_Y)
                if not node.IsFixed(DISPLACEMENT_Z):
                    free_node_list_z.append(node)
                    node.Fix(DISPLACEMENT_Z)
    #        for node in self.model_part.Nodes:
    #            node.Fix(DISPLACEMENT_X)
    #            node.Fix(DISPLACEMENT_Y)
    #            node.Fix(DISPLACEMENT_Z)

            # solve for phase field
            self.model_part.ProcessInfo[FRACTIONAL_STEP] = 1
            self.SolveOneStep()
            print("Solve for phase field completed")

            # unfreeze the displacement dofs
            for node in free_node_list_x:
                node.Free(DISPLACEMENT_X)
            for node in free_node_list_y:
                node.Free(DISPLACEMENT_Y)
            for node in free_node_list_z:
                node.Free(DISPLACEMENT_Z)

            ################################################
            it = it + 1

            # TODO compute the energy functional

#            print("################PHASE FIELD CONVERGENCE CRITERIA##################")
#            print("displacement: norm_disp_changed = " + str(norm_disp_changed) + "\t, norm_disp = " + str(norm_disp))
#            print("phase field : norm_phase_field_changed = " + str(norm_phase_field_changed) + "\t, norm_phase_field = " + str(norm_phase_field))
#            print("current total phase field problem convergence criteria = " + str(current_tol) + ", tol = " + str(self.toll))
#            print("##################################################################")

#        if(it == staggered_max_iter or converged == False):
#            print("Staggered Iteration does not converge in " + str(self.max_iter) + " steps")
#            print("Solve for staggered solution failed")
#            sys.exit(0)

    #
    def SolveOneStep(self):
        # perform the operations to be performed ONCE and ensure they will not be repeated
        # elemental function "Initialize" is called here
        if(self.InitializeWasPerformed == False):
            self.Initialize()
            self.InitializeWasPerformed = True

        # perform initializations for the current step
        # this operation implies:
        # identifying the set of DOFs that will be solved during this step
        # organizing the DOFs so to identify the dirichlet conditions
        # resizing the matrix preallocating the "structure"
        if (self.SolutionStepIsInitialized == False):
            if(self.builder_and_solver.GetDofSetIsInitializedFlag() == False or self.ReformDofSetAtEachStep == True):
                reform_dofs = True
            else:
                reform_dofs = False
            self.InitializeSolutionStep(reform_dofs)
            self.SolutionStepIsInitialized = True

        # perform prediction
        self.Predict()

        # execute iteration - first iteration is ALWAYS executed
        calculate_norm = False
        normDx = self.ExecuteIteration(self.echo_level, self.MoveMeshFlag, calculate_norm)
        it = 1

        # non linear loop
        converged = False

        while(it < self.max_iter and converged == False):
            # verify convergence
            converged = self.conv_criteria.PreCriteria(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)

            # calculate iteration
            # - system is built and solved
            # - database is updated depending on the solution
            # - nodal coordinates are updated if required
            normDx = self.ExecuteIteration(self.echo_level, self.MoveMeshFlag, calculate_norm)

            # verify convergence
            converged = self.conv_criteria.PostCriteria(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)

            # update iteration count
            it = it + 1

        # finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)
        self.SolutionStepIsInitialized = False

        self.iterations_last_solution = it - 1

        # clear if needed - deallocates memory
        mpi.world.barrier()
        if(self.ReformDofSetAtEachStep):
            self.Clear()
        if(mpi.rank == 0):
            print("SolveOneStep is Finished")

        if(it == self.max_iter or converged == False):
            print("Iteration does not converge in " + str(self.max_iter) + " steps")
            if(self.model_part.ProcessInfo[FRACTIONAL_STEP] == 0):
                print("Solve for displacements failed at time " + str(self.model_part.ProcessInfo[TIME]))
            elif(self.model_part.ProcessInfo[FRACTIONAL_STEP] == 1):
                print("Solve for phase field failed at time " + str(self.model_part.ProcessInfo[TIME]))
#            sys.exit(0)
        return converged

    #
    def Predict(self):
        self.time_scheme.Predict(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b)

    #
    def InitializeSolutionStep(self, reform_dofs):
        if(reform_dofs):
            # initialize the list of degrees of freedom to be used
            self.builder_and_solver.SetUpDofSet(self.time_scheme, self.model_part)
            # reorder the list of degrees of freedom to identify fixity and system size
            self.builder_and_solver.SetUpSystem(self.model_part)
            # allocate memory for the system and preallocate the structure of the matrix
            self.builder_and_solver.ResizeAndInitializeVectors(self.pA, self.pDx, self.pb, self.model_part.Elements, self.model_part.Conditions, self.model_part.ProcessInfo)

            # updating references
            self.A = (self.pA).GetReference()
            self.Dx = (self.pDx).GetReference()
            self.b = (self.pb).GetReference()

            # clear scheme so dof update map is recomputed for the new Dof set
            self.time_scheme.Clear()
            
        self.builder_and_solver.InitializeSolutionStep(self.model_part, self.A, self.Dx, self.b)
        self.time_scheme.InitializeSolutionStep(self.model_part, self.A, self.Dx, self.b)

    #
    def ExecuteIteration(self, echo_level, MoveMeshFlag, CalculateNormDxFlag):
        # reset system matrices and vectors prior to rebuild
        self.parallel_space.SetToZeroMatrix(self.A)
        self.parallel_space.SetToZeroVector(self.Dx)
        self.parallel_space.SetToZeroVector(self.b)

        self.time_scheme.InitializeNonLinIteration(self.model_part, self.A, self.Dx, self.b)

        # build and solve the problem
        self.builder_and_solver.Build(self.time_scheme, self.model_part, self.A, self.b)
        self.builder_and_solver.ApplyDirichletConditions(self.time_scheme, self.model_part, self.A, self.Dx, self.b)
        if self.linear_solver.AdditionalPhysicalDataIsNeeded():
            self.linear_solver.ProvideAdditionalData(self.A, self.Dx, self.b, self.builder_and_solver.GetDofSet(), self.model_part)
        self.builder_and_solver.SystemSolve(self.A, self.Dx, self.b)
#        self.builder_and_solver.BuildAndSolve(self.time_scheme, self.model_part, self.A, self.Dx, self.b)


        # perform update
        self.time_scheme.Update(self.model_part, self.builder_and_solver.GetDofSet(), self.A, self.Dx, self.b);

        # move the mesh as needed
        if(MoveMeshFlag):
            self.time_scheme.MoveMesh(self.model_part.Nodes);

        self.time_scheme.FinalizeNonLinIteration(self.model_part, self.A, self.Dx, self.b)

        # calculate the norm of the "correction" Dx
        if(CalculateNormDxFlag):
            normDx = self.parallel_space.TwoNorm(self.Dx)
        else:
            normDx = 0.0

        return normDx

    #
    def FinalizeSolutionStep(self, CalculateReactionsFlag):
        if(CalculateReactionsFlag):
            self.builder_and_solver.CalculateReactions(self.time_scheme, self.model_part, self.A, self.Dx, self.b)
        # Finalisation of the solution step,
        self.time_scheme.FinalizeSolutionStep(self.model_part, self.A, self.Dx, self.b)
        self.builder_and_solver.FinalizeSolutionStep(self.model_part, self.A, self.Dx, self.b)
        self.time_scheme.Clean()
        # reset flags for the next step
        self.mSolutionStepIsInitialized = False

    #
    def Clear(self):
        mpi.world.barrier()
        if(mpi.rank == 0):
            print("Entered in Clear")
        self.parallel_space.ClearMatrix(self.pA)
        self.parallel_space.ClearVector(self.pDx)
        self.parallel_space.ClearVector(self.pb)

        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()

        self.builder_and_solver.SetDofSetIsInitializedFlag(False)

        self.builder_and_solver.Clear()

        self.time_scheme.Clear()

        if(mpi.rank == 0):
            print("Clear is completed")

    #
    def SetEchoLevel(self, level):
        self.echo_level = level
        self.builder_and_solver.SetEchoLevel(level)

