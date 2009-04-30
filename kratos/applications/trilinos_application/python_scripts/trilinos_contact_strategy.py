#importing the Kratos Library
from Kratos import *
from KratosTrilinosApplication import *
from KratosStructuralApplication import *
import os,time
import mpi

import trilinos_strategy_python

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

class SolvingStrategyPython(trilinos_strategy_python.SolvingStrategyPython):
    #######################################################################
    def __init__(self,builder_and_solver_type,model_part,time_scheme,linear_solver,convergence_criteria,CalculateReactionsFlag,ReformDofSetAtEachStep,MoveMeshFlag,Comm,guess_row_size,Parameters):
        trilinos_strategy_python.SolvingStrategyPython.__init__( self, builder_and_solver_type, model_part, time_scheme, linear_solver, convergence_criteria, CalculateReactionsFlag, ReformDofSetAtEachStep, MoveMeshFlag, Comm, guess_row_size )
        #save the input parameters
        self.Parameters = Parameters
        self.PerformContactAnalysis = self.Parameters[0]
        #contact utility
        self.cu = ContactUtility( 3 )        
            
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
        if( self.cu.IsConverged( self.model_part, 0,  originalPosition, self.Parameters[3] ) == True ):
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
        (self.builder_and_solver).SetReshapeMatrixFlag(False)
        #props = self.model_part.Properties[1]
        for uzawaStep in range(1, self.Parameters[2] ):
            print "I am inside the uzawa loop, iteration no. " + str(uzawaStep)
            ## solving the standard newton-raphson iteration 
            self.PerformNewtonRaphsonIteration()
            ## updating the lagrange multipliers
            self.cu.Update( self.model_part, originalPosition, self.Parameters[3], self.Parameters[6], self.Parameters[8], self.Parameters[11], self.Parameters[9], self.Parameters[12], self.Parameters[7], self.Parameters[10]  )
            ## checking convergence
            if( self.cu.IsConverged( self.model_part, uzawaStep, originalPosition, self.Parameters[3] ) == True ):
                uzawaConverged = True
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
        
    def PerformNewtonRaphsonIteration( self ):
        #perform the operations to be performed ONCE and ensure they will not be repeated
        # elemental function "Initialize" is called here
        if(self.InitializeWasPerformed == False):
            self.Initialize()
            self.InitializeWasPerformed = True
        #perform initializations for the current step
        #this operation implies:
        #identifying the set of DOFs that will be solved during this step
        #organizing the DOFs so to identify the dirichlet conditions
        #resizing the matrix preallocating the "structure"
        if (self.SolutionStepIsInitialized == False):
            self.InitializeSolutionStep( self.ReformDofSetAtEachStep )
            self.SolutionStepIsInitialized = True
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
    #######################################################################

