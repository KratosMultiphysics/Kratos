#importing the Kratos Library
from __future__ import print_function
from KratosMultiphysics import *
from KratosMultiphysics.PFEM2Application import *


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
  	#self.builder_and_solver.SetEchoLevel(1)

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
        
    #######################################################################
    def Initialize(self):
	if(self.scheme.SchemeIsInitialized() == False):
	    self.scheme.Initialize(self.model_part)
			
	if (self.scheme.ElementsAreInitialized() == False): 
	    self.scheme.InitializeElements(self.model_part)


    #######################################################################
    #######################################################################
    def Solve(self,model_part,moveparticles,VariableUtils,CalculatePressureProjection,max_iter):
        #perform the operations to be performed ONCE and ensure they will not be repeated
        # elemental function "Initialize" is called here
        if(self.InitializeWasPerformed == False):
            self.Initialize()
            self.InitializeWasPerformed = True

        #perform initializations for the current step
        #this operation implie
        #identifying the set of DOFs that will be solved during this step
        #organizing the DOFs so to identify the dirichlet conditions
        #resizing the matrix preallocating the "structure"
        model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 10)
        if (self.SolutionStepIsInitialized == False):
            if(self.builder_and_solver.GetDofSetIsInitializedFlag() == False or self.ReformDofSetAtEachStep == True):
                reform_dofs = True
            else:
                reform_dofs = False
            self.InitializeSolutionStep(reform_dofs);
            self.SolutionStepIsInitialized = True;

        #perform prediction 
        self.Predict()

        #execute iteration - only one iteration is done since it is a linear system.
        calculate_norm = True
        non_linear_iteration_number = 1
        model_part.ProcessInfo.SetValue(NL_ITERATION_NUMBER, 1)
        normDx = self.ExecuteIteration(self.echo_level,self.MoveMeshFlag,calculate_norm,model_part,moveparticles,VariableUtils,CalculatePressureProjection)
        
        
        print ("fist iteration is done")
        #non linear loop
        converged = False

        while(non_linear_iteration_number < max_iter and converged == False):
            non_linear_iteration_number =  non_linear_iteration_number +1 
            model_part.ProcessInfo.SetValue(NL_ITERATION_NUMBER, non_linear_iteration_number)
            converged = self.convergence_criteria.PreCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)
            normDx = self.ExecuteIteration(self.echo_level,self.MoveMeshFlag,calculate_norm,model_part,moveparticles,VariableUtils,CalculatePressureProjection)
            converged = self.convergence_criteria.PostCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

        print("converged in ",non_linear_iteration_number," iterations")
        #finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)

        self.SolutionStepIsInitialized = False
        
        #clear if needed - deallocates memory 
        if(self.ReformDofSetAtEachStep == True):
                self.Clear();

            
    #######################################################################
    #######################################################################
    def Solve(self,model_part,moveparticles,VariableUtils,CalculatePressureProjection,CalculateVolumetricStrain,max_iter): #CalculatePressure,maxiter
        #perform the operations to be performed ONCE and ensure they will not be repeated
        # elemental function "Initialize" is called here
        if(self.InitializeWasPerformed == False):
            self.Initialize()
            self.InitializeWasPerformed = True
        
        #perform initializations for the current step
        #this operation implie
        #identifying the set of DOFs that will be solved during this step
        #organizing the DOFs so to identify the dirichlet conditions
        #resizing the matrix preallocating the "structure"
        model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 10)
        if (self.SolutionStepIsInitialized == False):
            if(self.builder_and_solver.GetDofSetIsInitializedFlag() == False or self.ReformDofSetAtEachStep == True):
                reform_dofs = True
            else:
                reform_dofs = False
            self.InitializeSolutionStep(reform_dofs);
            self.SolutionStepIsInitialized = True;

        #perform prediction 
        self.Predict()

        #execute iteration - only one iteration is done since it is a linear system.
        calculate_norm = True
        non_linear_iteration_number = 1
        model_part.ProcessInfo.SetValue(NL_ITERATION_NUMBER, 1)
        normDx = self.ExecuteIteration(self.echo_level,self.MoveMeshFlag,calculate_norm,model_part,moveparticles,VariableUtils,CalculatePressureProjection,CalculateVolumetricStrain)#,CalculatePressure)
        
        
        print ("fist iteration is done")
        #non linear loop
        converged = False

        while(non_linear_iteration_number < max_iter and converged == False):
            non_linear_iteration_number =  non_linear_iteration_number +1 
            model_part.ProcessInfo.SetValue(NL_ITERATION_NUMBER, non_linear_iteration_number)
            converged = self.convergence_criteria.PreCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)
            normDx = self.ExecuteIteration(self.echo_level,self.MoveMeshFlag,calculate_norm,model_part,moveparticles,VariableUtils,CalculatePressureProjection,CalculateVolumetricStrain)#,CalculatePressure)
            converged = self.convergence_criteria.PostCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

        print("converged in ",non_linear_iteration_number," iterations")
        #finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)

        self.SolutionStepIsInitialized = False
        
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
    def ExecuteIteration(self,echo_level,MoveMeshFlag,CalculateNormDxFlag,model_part,moveparticles,VariableUtils,CalculatePressureProjection):

        #implicit everything
        #self.CalculatePressureProjection()
        model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 20)
        (moveparticles).ComputeDeltaVelocityForNonLinearIteration();	 #saving velocity from previous iteration
        
        #########################
        model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 10)
        self.scheme.InitializeNonLinIteration(self.model_part,self.A,self.Dx,self.b)
        self.space_utils.SetToZeroVector(self.Dx)			
        self.space_utils.SetToZeroVector(self.b)
        self.space_utils.SetToZeroMatrix(self.A)
        
        self.builder_and_solver.BuildAndSolve(self.scheme,self.model_part,self.A,self.Dx,self.b)
        self.scheme.Update(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b);
        #move the mesh as needed
        if(MoveMeshFlag == True):
            self.scheme.MoveMesh(self.model_part.Nodes);
        self.scheme.FinalizeNonLinIteration(self.model_part,self.A,self.Dx,self.b)
        #########################
        CalculatePressureProjection()
        #bulk_modulus = 207000000.0 / (3.0 * (1.0-2.0* 0.0))
        #for node in self.model_part.Nodes:
        #   pressure = node.GetSolutionStepValue(PRESSURE,1) - node.GetSolutionStepValue(VOLUMETRIC_STRAIN)*bulk_modulus
        #   node.SetSolutionStepValue(PRESSURE,pressure)
        full_reset=True;
        (moveparticles).ResetBoundaryConditions(full_reset) 
        #(moveparticles).UpdateParticleStresses()	
        #finally we save the last pressure for future iterations. No need to do it for the velocity since it is saved in delta_velocity
        (VariableUtils).CopyScalarVar(PRESSURE,PRESSUREAUX,self.model_part.Nodes) #we save the pressure of this iteration in order to have the delta_pressure in the next one


        ########################
        #calculate the norm of the "correction" Dx
        if(CalculateNormDxFlag == True):
            normDx = self.space_utils.TwoNorm(self.Dx)
        else:
            normDx = 0.0
            
        return normDx
        
    #######################################################################
    def ExecuteIteration(self,echo_level,MoveMeshFlag,CalculateNormDxFlag,model_part,moveparticles,VariableUtils,CalculatePressureProjection,CalculateVolumetricStrain):#,CalculatePressure):

        #implicit everything
        #self.CalculatePressureProjection()
        model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 20)
        (moveparticles).ComputeDeltaVelocityForNonLinearIteration();	 #saving velocity from previous iteration
        
        #########################
        model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 10)
        self.scheme.InitializeNonLinIteration(self.model_part,self.A,self.Dx,self.b)
        self.space_utils.SetToZeroVector(self.Dx)			
        self.space_utils.SetToZeroVector(self.b)
        self.space_utils.SetToZeroMatrix(self.A)
        
        self.builder_and_solver.BuildAndSolve(self.scheme,self.model_part,self.A,self.Dx,self.b)
        self.scheme.Update(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b);
        #move the mesh as needed
        if(MoveMeshFlag == True):
            self.scheme.MoveMesh(self.model_part.Nodes);
        self.scheme.FinalizeNonLinIteration(self.model_part,self.A,self.Dx,self.b)
        #########################
        CalculatePressureProjection()
        #CalculateVolumetricStrain()
        #bulk_modulus = 207000000.0 / (3.0 * (1.0-2.0* 0.3))
        #for node in self.model_part.Nodes:
        #   pressure = node.GetSolutionStepValue(PRESSURE,1) - node.GetSolutionStepValue(VOLUMETRIC_STRAIN)*bulk_modulus
        #   node.SetSolutionStepValue(PRESSURE,pressure)
        full_reset=True;
        (moveparticles).ResetBoundaryConditions(full_reset) 
        (moveparticles).UpdateParticleStresses()
        #CalculatePressure()	
        #finally we save the last pressure for future iterations. No need to do it for the velocity since it is saved in delta_velocity
        (VariableUtils).CopyScalarVar(PRESSURE,PRESSUREAUX,self.model_part.Nodes) #we save the pressure of this iteration in order to have the delta_pressure in the next one


        ########################
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
