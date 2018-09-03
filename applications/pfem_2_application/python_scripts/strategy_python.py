#importing the Kratos Library
from KratosMultiphysics import *
import time
        
class SolvingStrategyPython:
    def __init__(self,model_part,time_scheme,linear_solver,convergence_criteria,CalculateReactionsFlag,ReformDofSetAtEachStep,MoveMeshFlag):
        #save the input parameters
        self.model_part = model_part
        self.scheme = time_scheme
        self.linear_solver = linear_solver
        self.convergence_criteria = convergence_criteria
        self.CalculateReactionsFlag = CalculateReactionsFlag
        self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
        self.MoveMeshFlag = MoveMeshFlag

        self.space_utils = UblasSparseSpace()

        #default values for some variables
        self.rebuild_level = 1 #rebuild at each solution step
        self.echo_level = 1        
        
        self.builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        self.builder_and_solver.SetEchoLevel(0)

        #local matrices and vectors
        self.pA = self.space_utils.CreateEmptyMatrixPointer()
        self.pDx = self.space_utils.CreateEmptyVectorPointer()
        self.pb = self.space_utils.CreateEmptyVectorPointer()

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

        #by default is linear. so only one iteration will be performed.
        self.max_iter = 1
        self.computepressureprojection=False;
        
    #######################################################################
    def Initialize(self):
        if(self.scheme.SchemeIsInitialized() == False):
            self.scheme.Initialize(self.model_part)
                        
        if (self.scheme.ElementsAreInitialized() == False): 
            self.scheme.InitializeElements(self.model_part)


    #######################################################################
    #######################################################################
    def Solve(self):
        #perform the operations to be performed ONCE and ensure they will not be repeated
        # elemental function "Initialize" is called here
        if(self.InitializeWasPerformed == False):
            self.Initialize()
            self.InitializeWasPerformed = True
        t1=time.time()
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
        #print("assembly time = ",time.time()-t1)
        #perform prediction 
        self.Predict()

        #execute iteration - first iteration is ALWAYS executed
        calculate_norm = False
        print ("about to perform first iteration")
        normDx = self.ExecuteIteration(self.echo_level,self.MoveMeshFlag,calculate_norm)
        it = 1
        
        #print ("fist iteration is done")
        if(self.computepressureprojection==True):
              self.pressureprojector()
        #non linear loop
        converged = False
        
        while(it < self.max_iter and converged == False): #by default self.max_iter=1, so we would not enter in the non linear loop
            #verify convergence
            converged = self.convergence_criteria.PreCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)
           
            #calculate iteration
            # - system is built and solved
            # - database is updated depending on the solution
            # - nodal coordinates are updated if required
            normDx = self.ExecuteIteration(self.echo_level,self.MoveMeshFlag,calculate_norm)

            #verify convergence
            converged = self.convergence_criteria.PostCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)
            if(self.computepressureprojection==True):
                  self.pressureprojector()
            #update iteration count
            it = it + 1

        #finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)
        self.SolutionStepIsInitialized = False
        
        if(self.max_iter>1):
             if(it < self.max_iter):
                 print ("converged in ", it, " iterations")       
             else:
                 print ("not converged after ",it," iterations")
        #clear if needed - deallocates memory 
        if(self.ReformDofSetAtEachStep == True):
            self.Clear();

            
    #######################################################################
    #######################################################################

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
            self.builder_and_solver.ResizeAndInitializeVectors(self.scheme, self.pA,self.pDx,self.pb,self.model_part);

            #updating references
            self.A = (self.pA).GetReference()
            self.Dx = (self.pDx).GetReference()
            self.b = (self.pb).GetReference()


            
        self.builder_and_solver.InitializeSolutionStep(self.model_part,self.A,self.Dx,self.b)
        self.scheme.InitializeSolutionStep(self.model_part,self.A,self.Dx,self.b)

    #######################################################################
    def ExecuteIteration(self,echo_level,MoveMeshFlag,CalculateNormDxFlag):
        #reset system matrices and vectors prior to rebuild

        

        self.scheme.InitializeNonLinIteration(self.model_part,self.A,self.Dx,self.b)

        self.space_utils.SetToZeroVector(self.Dx)                        
        self.space_utils.SetToZeroVector(self.b)
        #build and solve the problem
        self.space_utils.SetToZeroMatrix(self.A)
                       
        self.builder_and_solver.BuildAndSolve(self.scheme,self.model_part,self.A,self.Dx,self.b)
        
        #full output if needed
        #echo_level=4
        if(echo_level == 3):
            print("SystemMatrix = ", self.A)
            print("solution obtained = ", self.Dx)
            print("RHS = ", self.b)
        elif(echo_level == 4):
             filename = str("A_t_") + str(self.model_part.ProcessInfo[TIME]) + str(".mm");
             WriteMatrixMarketMatrix(filename, self.A, False);
             bname = str("b_t_") + str(self.model_part.ProcessInfo[TIME]) + str(".mm");
             #print "solution obtained = ", self.Dx 
             #print "RHS = ", self.b
        WriteMatrixMarketMatrix     
        #perform update
        self.scheme.Update(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b);

        #move the mesh as needed
        if(MoveMeshFlag == True):
            self.scheme.MoveMesh(self.model_part.Nodes);

        self.scheme.FinalizeNonLinIteration(self.model_part,self.A,self.Dx,self.b)

        if(echo_level == 3):
            print("after finalize iteration")
            print("SystemMatrix = ", self.A)
            print("solution obtained = ", self.Dx) 
            print("RHS = ", self.b)
    
        
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
        self.mSolutionStepIsInitialized = False

    #######################################################################
    def Clear(self):
        self.space_utils.ClearMatrix(self.pA)
        self.space_utils.ResizeMatrix(self.A,0,0)
        
        self.space_utils.ClearVector(self.pDx)
        self.space_utils.ResizeVector(self.Dx,0)

        self.space_utils.ClearVector(self.pb)
        self.space_utils.ResizeVector(self.b,0)

        #updating references
        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()
        
        self.builder_and_solver.SetDofSetIsInitializedFlag(False)
        self.builder_and_solver.Clear()

    #######################################################################   
    def SetEchoLevel(self,level):
        self.echo_level = level
        self.builder_and_solver.SetEchoLevel(level)

    def SetMaximumIterations(self,max_iter):
        self.max_iter=max_iter

    def ActivatePressureProjectionCalculation(self, CalculatePressureProjection ):
        self.pressureprojector =  CalculatePressureProjection
        self.computepressureprojection=True;
