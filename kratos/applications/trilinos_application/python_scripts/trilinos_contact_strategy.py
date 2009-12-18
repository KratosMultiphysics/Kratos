#importing the Kratos Library
from Kratos import *
from KratosStructuralApplication import *
from KratosTrilinosApplication import *
import os,time

import mpi

## Parameters contain:
# perform_contact_analysis_flag
# penalty value for normal contact
# maximum number of uzawa iterations
# friction coefficient
# penalty value for frictional contact
# contact_double_check_flag
# contact_ramp_penalties_flag
# maximum penalty value for normal contact
# ramp criterion for normal contact
# ramp factor for normal contact
# maximum penalty value for frictional contact
# ramp criterion for frictional contact
# ramp factor for frictional contact

class SolvingStrategyPython:
    #######################################################################
    def __init__( self, builder_and_solver_type, model_part, time_scheme, linear_solver, convergence_criteria, CalculateReactionsFlag, ReformDofSetAtEachStep, MoveMeshFlag, Parameters, space_utils, Comm, guess_row_size ):
        #save the input parameters
        self.model_part = model_part
        self.scheme = time_scheme
        self.linear_solver = linear_solver
        self.convergence_criteria = convergence_criteria
        self.CalculateReactionsFlag = CalculateReactionsFlag
        self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
        self.MoveMeshFlag = MoveMeshFlag
        self.Parameters = Parameters
        self.PerformContactAnalysis = self.Parameters[0]
        self.PrintSparsity = self.Parameters[13]
        self.space_utils = space_utils
        #contact utility
        self.cu = ContactUtility( 3 )
        #default values for some variables
        self.max_iter = 30
        self.echo_level = 1
        if(builder_and_solver_type == "standard"):
            self.builder_and_solver = TrilinosResidualBasedBuilderAndSolver(Comm,guess_row_size,self.linear_solver)
        elif(builder_and_solver_type == "ML2D"):
            self.builder_and_solver = TrilinosBuilderAndSolverML2D(Comm,guess_row_size,2,self.linear_solver)
        elif(builder_and_solver_type == "ML3D"):
            self.builder_and_solver = TrilinosBuilderAndSolverML2D(Comm,guess_row_size,3,self.linear_solver)
        elif(builder_and_solver_type == "ML2Dpress"):
            self.builder_and_solver = TrilinosBuilderAndSolverMLmixed(Comm,guess_row_size,2,self.linear_solver)
        elif(builder_and_solver_type == "ML3Dpress"):
            self.builder_and_solver = TrilinosBuilderAndSolverMLmixed(Comm,guess_row_size,3,self.linear_solver)
        elif(builder_and_solver_type == "superludist"):
            self.builder_and_solver = TrilinosResidualBasedBuilderAndSolver(Comm,guess_row_size,self.linear_solver)
        elif(builder_and_solver_type == "MLdeactivation"):
            self.builder_and_solver = TrilinosBuilderAndSolverMLDeactivation2D(Comm,guess_row_size,3,self.linear_solver)
        elif(builder_and_solver_type == "superludist_deactivation"):
            self.builder_and_solver = TrilinosResidualBasedBuilderAndSolverDeactivation(Comm,guess_row_size,self.linear_solver)


        
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
        (self.builder_and_solver).SetCalculateReactionsFlag(self.CalculateReactionsFlag);
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep);
        
    #######################################################################
    def Initialize(self):
        if(self.scheme.SchemeIsInitialized() == False):
            self.scheme.Initialize(self.model_part)
        if (self.scheme.ElementsAreInitialized() == False): 
            self.scheme.InitializeElements(self.model_part)

            
    #######################################################################
    def Solve(self):
        #print self.model_part
        ## - storing original condition size before adding virtual conditions.
        ## - performing contact search
        ## - creating virtual link conditions for the assembling
        if( self.PerformContactAnalysis == False ):
            self.PerformNewtonRaphsonIteration()
            #finalize the solution step
            self.FinalizeSolutionStep(self.CalculateReactionsFlag)
            #clear if needed - deallocates memory 
            if(self.ReformDofSetAtEachStep == True):
                self.Clear();
            return
        print "setting up contact conditions"
        originalPosition =  self.cu.SetUpContactConditions(self.model_part, self.Parameters[1], self.Parameters[4], self.Parameters[5] )
        uzawaConverged = False
        ##  First step: reform DOF set and check if uzawa iteration is necessary
        self.PerformNewtonRaphsonIteration()
        self.cu.Update( self.model_part, originalPosition, self.Parameters[3], self.Parameters[6], self.Parameters[8], self.Parameters[11], self.Parameters[9], self.Parameters[12], self.Parameters[7], self.Parameters[10]  )
        ## CHECK FOR CONTACT CONVERGENCE AMONG ALL PROCESSES
        contact_converged = self.cu.IsConverged( self.model_part, 0, originalPosition, self.Parameters[3] )
        print("rank "+str(mpi.rank)+": contact converged = "+str(contact_converged) )
        all_contact_converged = mpi.all_gather( mpi.world, contact_converged )
        print("in rank "+str(mpi.rank)+": "+str(all_contact_converged) )
        for flag in all_contact_converged:
            if( flag == False ):
                contact_converged = False
        mpi.world.barrier()
        if( contact_converged ):
            uzawaConverged = True
            (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)
            self.cu.Clean( self.model_part, originalPosition );
            #finalize the solution step
            self.FinalizeSolutionStep(self.CalculateReactionsFlag)
            #clear if needed - deallocates memory 
            if(self.ReformDofSetAtEachStep == True):
                self.Clear();
            return

        ## beginning of UZAWA loop
        print("##################")
        print("UZAWA LOOP STARTED")
        print("##################")
        (self.builder_and_solver).SetReshapeMatrixFlag(False)
        #props = self.model_part.Properties[1]
        for uzawaStep in range(1, self.Parameters[2] ):
            print "I am inside the uzawa loop, iteration no. " + str(uzawaStep)
            ## solving the standard newton-raphson iteration 
            self.PerformNewtonRaphsonIteration()
            ## updating the lagrange multipliers
            self.cu.Update( self.model_part, originalPosition, self.Parameters[3], self.Parameters[6], self.Parameters[8], self.Parameters[11], self.Parameters[9], self.Parameters[12], self.Parameters[7], self.Parameters[10]  )
            ## checking convergence
            uzawaConverged = self.cu.IsConverged(self.model_part, uzawaStep, originalPosition, self.Parameters[3])
            all_uzawa_converged = mpi.all_gather( mpi.world, uzawaConverged )
            for flag in all_uzawa_converged:
                if( flag == False ):
                    uzawaConverged = False
            mpi.world.barrier()
            if( uzawaConverged ):
                break
        if( uzawaConverged == False ):
            print "uzawa algorithm failes to converge within maximum number of iterations"
        ## end of UZAWA loop
        ## cleaning up the conditions
        self.cu.Clean( self.model_part, originalPosition )
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)
        #finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)
        #clear if needed - deallocates memory 
        if(self.ReformDofSetAtEachStep == True):
            self.Clear();
        #print "Solve is Finished for rank : ", mpi.rank


    #######################################################################
    #######################################################################
    def PerformNewtonRaphsonIteration( self ):
        #perform initializations for the current step
        #this operation implies:
        #identifying the set of DOFs that will be solved during this step
        #organizing the DOFs so to identify the dirichlet conditions
        #resizing the matrix preallocating the "structure"
        if (self.SolutionStepIsInitialized == False):
            if(self.builder_and_solver.GetDofSetIsInitializedFlag() == False or self.ReformDofSetAtEachStep == True):
                reform_dofs = True
            else:
                reform_dofs = False
	    self.InitializeSolutionStep(reform_dofs);
            self.SolutionStepIsInitialized = True;

        #perform prediction 
        self.Predict()

        #execute iteration - first iteration is ALWAYS executed
        calculate_norm = False
        normDx = self.ExecuteIteration(self.echo_level,self.MoveMeshFlag,calculate_norm)
        it = 1

        #non linear loop
        converged = False
        
        while(it < self.max_iter and converged == False):
            #verify convergence
            converged = self.convergence_criteria.PreCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)
           
            #calculate iteration
            # - system is built and solved
            # - database is updated depending on the solution
            # - nodal coordinates are updated if required
            normDx = self.ExecuteIteration(self.echo_level,self.MoveMeshFlag,calculate_norm)

            #verify convergence
            converged = self.convergence_criteria.PostCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

           
            #update iteration count
            it = it + 1
            if( it == self.max_iter ):
                print("Iteration did not converge")

    #######################################################################
    def Predict(self):
        self.scheme.Predict(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b);

    #######################################################################
    def InitializeSolutionStep(self,reform_dofs):
        if(self.builder_and_solver.GetDofSetIsInitializedFlag() == False or self.ReformDofSetAtEachStep == True):
            #initialize the list of degrees of freedom to be used
            self.builder_and_solver.SetUpDofSet(self.scheme,self.model_part);
            #reorder the list of degrees of freedom to identify fixity and system size	  			
            self.builder_and_solver.SetUpSystem(self.model_part)
            #allocate memory for the system and preallocate the structure of the matrix
            self.builder_and_solver.ResizeAndInitializeVectors(self.pA,self.pDx,self.pb,self.model_part.Elements,self.model_part.Conditions,self.model_part.ProcessInfo);

            #updating references
            self.A = (self.pA).GetReference()
            self.Dx = (self.pDx).GetReference()
            self.b = (self.pb).GetReference()


            
        self.builder_and_solver.InitializeSolutionStep(self.model_part,self.A,self.Dx,self.b)
        self.scheme.InitializeSolutionStep(self.model_part,self.A,self.Dx,self.b)


    #######################################################################
    def ExecuteIteration(self,echo_level,MoveMeshFlag,CalculateNormDxFlag):
        #reset system matrices and vectors prior to rebuild
        self.space_utils.SetToZeroMatrix(self.A)
        self.space_utils.SetToZeroVector(self.Dx)			
        self.space_utils.SetToZeroVector(self.b)
        self.scheme.InitializeNonLinIteration(self.model_part,self.A,self.Dx,self.b)
        #build and solve the problem
        
        self.builder_and_solver.BuildAndSolve(self.scheme,self.model_part,self.A,self.Dx,self.b)
        #full output if needed
        if(echo_level >= 3):
            print "SystemMatrix = ", self.A 
            print "solution obtained = ", self.Dx 
            print "RHS = ", self.b
            
        #perform update
        self.scheme.Update(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b);
        #move the mesh as needed
        if(MoveMeshFlag == True):
            self.scheme.MoveMesh(self.model_part.Nodes);

        self.scheme.FinalizeNonLinIteration(self.model_part,self.A,self.Dx,self.b)
    
        
        #calculate the norm of the "correction" Dx
        if(CalculateNormDxFlag == True):
            normDx = self.space_utils.TwoNorm(self.Dx)
        else:
            normDx = 0.0
            
        return normDx
      
    #######################################################################
    def FinalizeSolutionStep(self,CalculateReactionsFlag):
        if(CalculateReactionsFlag == True):
            self.builder_and_solver.CalculateReactions(self.scheme,self.model_part,self.A,self.Dx,self.b)
        #Finalisation of the solution step, 
        self.scheme.FinalizeSolutionStep(self.model_part,self.A,self.Dx,self.b)
        self.builder_and_solver.FinalizeSolutionStep(self.model_part,self.A,self.Dx,self.b)
        self.scheme.Clean()
        #reset flags for the next step
        self.SolutionStepIsInitialized = False
  
    #######################################################################
    def Clear(self):
        self.space_utils.ClearMatrix(self.pA)
        self.space_utils.ClearVector(self.pDx)
        self.space_utils.ClearVector(self.pb)
        
        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()
    
        self.builder_and_solver.SetDofSetIsInitializedFlag(False)

        self.builder_and_solver.Clear()

    #######################################################################   
    def SetEchoLevel(self,level):
        self.echo_level = level
        self.builder_and_solver.SetEchoLevel(level)

    #######################################################################   

