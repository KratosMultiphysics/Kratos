#importing the Kratos Library
from Kratos import *
from KratosULFApplication import *
#import cProfile

class ULFFracStrategyPython:
    #######################################################################
    def __init__(self,fluid_only_model_part, combined_model_part, fluid_model_part, disp_time_scheme, pres_time_scheme, pres_linear_solver, convergence_criteria,CalculateReactionsFlag,ReformDofSetAtEachStep,MoveMeshFlag,domain_size, bulk_modulus, density):
        outstring2 = "convergence_info.txt"
        out_file = open(outstring2, 'w')
        
        self.out_file=out_file
        self.out_file.write("Bulk" + str(bulk_modulus)+"\n")
        #save the input parameters
        #the model_part is the combined model part
        self.model_part = combined_model_part
        self.fluid_model_part = fluid_model_part

        self.fluid_only_model_part = fluid_only_model_part
        
        self.bulk_modulus = bulk_modulus
        self.density = density
        self.K = self.density*self.bulk_modulus
        self.scheme = disp_time_scheme
        self.dummy_solver = SkylineLUFactorizationSolver()

        self.pres_time_scheme=pres_time_scheme
        self.pres_linear_solver=pres_linear_solver
        
        self.convergence_criteria = convergence_criteria
        self.CalculateReactionsFlag = CalculateReactionsFlag
        self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
        self.MoveMeshFlag = MoveMeshFlag
        self.domain_size = domain_size

        self.CalculateMass=MassCalculateProcess(self.model_part)
        
        self.space_utils = UblasSparseSpace()

        #default values for some variables
        self.max_iter = 40

        self.rebuild_level = 1 #rebuild at each solution step
        self.echo_level = 1
        if (domain_size==2):
            self.builder_and_solver = ResidualBasedIncompressibleBuilder2D(self.dummy_solver)
        if (domain_size==3):
            self.builder_and_solver = ResidualBasedEliminationQuasiIncompressibleBuilderAndSolver3D(self.dummy_solver)
            
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)

        #local matrices and vectors
        self.A = CompressedMatrix()
        self.Dx = Vector()
        self.b = Vector()

        #divergence matrix D, MPinv (pressure mass matrix, lumped)
        self.WorkMatrix = CompressedMatrix()
        
        self.D = CompressedMatrix()
        self.MPconsistent = CompressedMatrix()
        self.MPinv = Vector()

        #this one is for fluid only
        self.D_fluid = CompressedMatrix()
        self.WorkMatrix=CompressedMatrix()

        #initialize flags
        self.SolutionStepIsInitialized = False
        self.InitializeWasPerformed = False
        self.StiffnessMatrixIsBuilt = False

        #provide settings to the builder and solver
        (self.builder_and_solver).SetCalculateReactionsFlag(self.CalculateReactionsFlag);
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep);

        ##########################################
        ## strategy for solution for pressure
        
        ReformDofSet=True
        #maybe should be just "fluid_model_part" below, and not the combined
        #self.PressureLinStrat=ResidualBasedLinearStrategy(self.model_part, self.pres_time_scheme, self.pres_linear_solver, False, ReformDofSet, False, False)
        self.PressureLinStrat=LapModifiedLinearStrategy2D(self.fluid_only_model_part, self.pres_time_scheme, self.pres_linear_solver, False, ReformDofSet, False, False)
        #self.PressureLinStrat=ResidualBasedLinearStrategy(self.fluid_only_model_part, self.pres_time_scheme, self.pres_linear_solver, False, ReformDofSet, False, False)
        self.PressureLinStrat.SetEchoLevel(1)
        #############################################
        (self.VariableUtils) = VariableUtils()

        
    #######################################################################
    def SetEchoLevel(self,level):
        self.echo_level = level

    #######################################################################
    def Initialize(self):
	if(self.scheme.SchemeIsInitialized() == False):
	    self.scheme.Initialize(self.model_part)
			
	if (self.scheme.ElementsAreInitialized() == False): 
	    self.scheme.InitializeElements(self.model_part)
    
    #######################################################################
    #######################################################################
    def Solve(self,domain_size,UlfUtils):
        #perform the operations to be performed ONCE and ensure they will not be repeated
        # elemental function "Initialize" is called here
        if(self.InitializeWasPerformed == False):
            self.Initialize()
            self.InitializeWasPerformed = True

        
        #perform initializations for the current step
        import time
        step_initialization_start = time.clock()
        reform_dofs = True
        self.InitializeSolutionStep(reform_dofs);
        print "step initizalization time =",time.clock()-step_initialization_start
        

        #perform prediction 
        self.Predict()

        #check for inverted elements
        inverted_elements = False

        volume = (UlfUtils).CalculateVolume(self.model_part,self.domain_size)

        if(volume <= 0.00):
            inverted_elements = True
            print "INVERTED ELEMENT FOUND - right after prediction"

        #execute iteration - first iteration is ALWAYS executed
        calculate_norm = False
        if(inverted_elements == False):
            print "INTEGRATION"
            normDx = self.ExecuteIteration(1, self.echo_level,self.MoveMeshFlag,calculate_norm, UlfUtils)
        it = 1
        
        print "size of problem!!!!!!!!!!!!!!! = ",len(self.b)
        
        #check for inverted elements      
        volume = (UlfUtils).CalculateVolume(self.model_part,self.domain_size)
        if(volume <= 0.00):
            inverted_elements = True
            print "INVERTED ELEMENT FOUND - just after first iteration"
     
        #non linear loop
        converged = False

        self.max_iter=40;
        global_iteration_number=1
        while(it < self.max_iter and converged == False and inverted_elements == False):
            #verify convergence
            converged = self.convergence_criteria.PreCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

            #calculate iteration
            # - system is built and solved
            # - database is updated depending on the solution
            # - nodal coordinates are updated if required
            
            #UlfUtils.CalculateNodalArea(self.fluid_model_part,self.domain_size);
            
            normDx = self.ExecuteIteration(global_iteration_number, self.echo_level,self.MoveMeshFlag,calculate_norm, UlfUtils)

            #verify convergence
            converged = self.convergence_criteria.PostCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

            #check for inverted elements      
            volume = (UlfUtils).CalculateVolume(self.model_part,self.domain_size)
            if(volume <= 0.00):
                inverted_elements = True
                print "INVERTED ELEMENT FOUND"
            
            #update iteration count
            it = it + 1
            
        (self.builder_and_solver).SavePressureIteration(self.model_part);
        for node in self.model_part.Nodes:
            if (node.GetSolutionStepValue(IS_STRUCTURE)==0):
                if (node.GetSolutionStepValue(IS_FREE_SURFACE)==1 and node.GetSolutionStepValue(IS_LAGRANGIAN_INLET)!=1):
		  # and node.X>0.5):
                  #node.SetSolutionStepValue(PRESSURE,0,0.0)
                  #print "FIXING PRESSURE AT THE OUTLET!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                  node.Fix(PRESSURE)
            if (node.GetSolutionStepValue(IS_INTERFACE)==1):
                node.Fix(PRESSURE)

        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 2);
        self.fluid_only_model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 2);          

        self.PressureLinStrat.Solve()        
        self.fluid_only_model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 1);          
        
        converged=False
        inverted_elements=False
        while(converged == False and inverted_elements == False):
            normDx = self.ExecuteIteration(3, self.echo_level,self.MoveMeshFlag,calculate_norm, UlfUtils)
            converged = self.convergence_criteria.PostCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)
            volume = (UlfUtils).CalculateVolume(self.model_part,self.domain_size)
            if(volume <= 0.00):
                inverted_elements = True
                print "INVERTED ELEMENT FOUND"
                

        

       # UlfUtils.CalculateNodalArea(self.model_part,self.domain_size);
        
        ###############################################################################################
        Atott=UlfUtils.CalculateVolume(self.model_part, self.domain_size)
        self.out_file.write(str(Atott)+" "+"\n")
        
        #moving lonely nodes
        UlfUtils.MoveLonelyNodes(self.model_part);

        #finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)
        
        #clear if needed - deallocates memory 
        self.Clear();

        return inverted_elements

            
    #######################################################################
    #######################################################################

    #######################################################################
    def Predict(self):
        self.scheme.Predict(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b);

    #######################################################################
    def InitializeSolutionStep(self,reform_dofs):
        if(reform_dofs == True):
            self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 1);
            print "reforming dofs"
            #initialize the list of degrees of freedom to be used 
            self.builder_and_solver.SetUpDofSet(self.scheme,self.model_part);
            #reorder the list of degrees of freedom to identify fixity and system size	  			
	    self.builder_and_solver.SetUpSystem(self.model_part)
	    #allocate memory for the system and preallocate the structure of the matrix
            self.builder_and_solver.ResizeAndInitializeVectors(self.A,self.D,self.Dx,self.b, self.MPconsistent, self.MPinv,self.model_part.Elements,self.model_part.Conditions,self.model_part.ProcessInfo);
            #auxiliary matrices (structure), needed SPECIFICALLY for the Quasi-Incompressible fluid
      	    
            self.builder_and_solver.ConstructMatrixStructure(self.A, self.model_part)
            self.builder_and_solver.ConstructMatrixStructure_DivergenceMatrixD(self.D,  self.model_part)
            self.builder_and_solver.ConstructMatrixStructure_Mconsistent(self.MPconsistent, self.model_part)
            
           
        self.builder_and_solver.InitializeSolutionStep(self.model_part,self.A,self.Dx,self.b)
        self.scheme.InitializeSolutionStep(self.model_part,self.A,self.Dx,self.b)
       
    #######################################################################
    def ExecuteIteration(self, global_iteration_number, echo_level,MoveMeshFlag,CalculateNormDxFlag, UlfUtils):
        #reset system matrices and vectors prior to rebuild
        self.space_utils.SetToZeroMatrix(self.A)
        self.space_utils.SetToZeroMatrix(self.MPconsistent)
        self.space_utils.SetToZeroMatrix(self.D)
        self.space_utils.SetToZeroVector(self.Dx)			
        self.space_utils.SetToZeroVector(self.b)
        self.space_utils.SetToZeroVector(self.MPinv)
        ############
        self.space_utils.SetToZeroMatrix(self.D_fluid)
        
        self.scheme.InitializeNonLinIteration(self.model_part,self.A,self.Dx,self.b)
        
    
        #build and solve the problem
        #build the standard part - here the GMPinvD is still not included!!!!
        #NOTE THAT BEFORE WE USED TO call the ConstructMatrixStructure in the InitAndResizeVecs - which we dont call anymore
        #


        #NOW WE SET THE FRACTIONAL STEP NUMBER TO 1 - to solve for displacement
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 1);
                                                         
        self.builder_and_solver.Build(self.scheme,self.model_part,self.A,self.b)
        #construct other necessary matrices: MPinv, G
        
        
        #assembling additional matrices
        self.builder_and_solver.AssembleMassMatrices(self.MPconsistent, self.MPinv, self.model_part)
        
        #now build the D = GT matrix
        self.builder_and_solver.BuildAuxiliaries(self.D, self.model_part)
        
        #print self.D
        #penalizing fixed DOF
        self.builder_and_solver.ModifyForDirichlet(self.A, self.b)
        
        #now we will use the CG algorith written here to solve the modified system: A+GMinvD=b
        self.prec_CG_Solve(20000, 1e-9, global_iteration_number)
        
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
        #updating pressures
	#parameter is the bulk modulus, and the second one is density
        
        UlfUtils.CalculateNodalArea(self.fluid_model_part,self.domain_size);

        if (global_iteration_number==1):
            print "Updating quasi-inc pressures"
            self.builder_and_solver.UpdatePressuresNew(self.MPconsistent, self.MPinv, self.model_part, self.bulk_modulus, self.density)
##                        
        self.scheme.FinalizeNonLinIteration(self.model_part,self.A,self.Dx,self.b)

        ##print "Found the U and now ready to do the Laplacian smoothing"

        print "WE FINISHED COMPUTING PRESSURES "
        
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
        #self.space_utils.ClearMatrix(self.A)
        self.space_utils.ResizeMatrix(self.A,0,0)

        #self.space_utils.ClearMatrix(self.D)
        self.space_utils.ResizeMatrix(self.D,0,0)

        #self.space_utils.ClearMatrix(self.MPconsistent)
        self.space_utils.ResizeMatrix(self.MPconsistent,0,0)
        
        
        #self.space_utils.ClearVector(self.Dx)
        self.space_utils.ResizeVector(self.Dx,0)

        #self.space_utils.ClearVector(self.MPinv)
        self.space_utils.ResizeVector(self.MPinv,0)

        
        #self.space_utils.ClearVector(self.b)
        self.space_utils.ResizeVector(self.b,0)
        
	self.builder_and_solver.SetDofSetIsInitializedFlag(False)
	self.builder_and_solver.Clear()
	print "clear completed"
    
        
    #####################################################################
    # due to the condensation of pressures we need to change the RHS - add the Gpn term
    # !!!!!!!!! NO NEED TO DO IT - RHS is already modified inside the element       
    ######################################################
    def calc_alpha(self, prod_ri_ri, prod_LHS_d, di):
        alpha=prod_ri_ri/self.space_utils.Dot(di, prod_LHS_d);
        return alpha;
    
    ######################################################
    def prec_CG_Solve(self, max_iterations, conv_criterion, global_iteration_number):

        print "STARTING CONJUGATE GRADIENT procedure"
        
        #first we compute the initial values of search direction d0 and residual r0
        #WE HAVE TO PASS THE INITIAL GUESS of "X0"
        #LHS_d is the variable to store the product of the LHS and the displacement vector (we just need the product. we
        #do not explicitely compute the LHS matrix (thats the advantage of our matrix-free method)
        LHS_d = Vector(self.space_utils.Size(self.b))

        di = Vector(self.space_utils.Size(self.b))
        prec_ri = Vector(self.space_utils.Size(self.b))
        preconditioner = Vector(self.space_utils.Size(self.b))
        destination = Vector(self.space_utils.Size(self.b))

        #claculating the preconditioner
        self.builder_and_solver.CalculatePreconditionerDiagonalMatrix(self.D, self.MPinv, self.A, preconditioner)
                
        #for the add_GMinvD_prod function 
        WorkArray = Vector(self.space_utils.Size(self.b)/self.domain_size)
        
        #we are assuming that the initial guess for the unknown X is 0 vecotr, =>r(0) = b - Ax = b
        #x0 = 0. mult by b just to get teh right size
          
        #here ri is r0 - before we enter the loop. we basically say: initialize ri=r0=d0
        #same is valid for xi=x0, and di = d0 =r0 in the beginning
        ri = Vector(self.b)
        
        #d0 = Prec_x_r0
        self.builder_and_solver.calc_prod_precond_vec(ri, preconditioner, di)

        xi = Vector(self.space_utils.Size(self.b))
        self.space_utils.SetToZeroVector(xi)
      
        counter = 0

        #kGM-1D should be added to the left hand side only at the first non-linear iteration!
        if (global_iteration_number<=2):
            print "ADDING GM-1D:........::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
            print "ADDING GM-1D:........::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
            print "ADDING GM-1D:........::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
            for i in range(1, max_iterations):
                #Calculate: PrecInv_ri
                self.builder_and_solver.calc_prod_precond_vec(ri, preconditioner, prec_ri)
                prod_ri_Prec_ri=self.space_utils.Dot(ri, prec_ri)
                #computing first (standard part, that is A_d) of LHS_d
                self.space_utils.Mult(self.A, di, LHS_d)
                #here we calculate the GMinvD_d term
                self.builder_and_solver.calc_GMinvD_prod(self.D, self.MPinv, di, WorkArray, destination)
                #print self.D
                #self.builder_and_solver.add_GMinvD_prod(self.D, self.MPinv, di, WorkArray, LHS_d, self.K)
                #and modifying the LHS by adding GM-1D_d term. NOTE, that result_aux IS REUSED for a different purpose below
                self.space_utils.UnaliasedAdd(LHS_d, -self.K, destination)

                #now alpha
                alpha = self.calc_alpha(prod_ri_Prec_ri, LHS_d,di);
                #updating/resetting xi,new residual ri, new search direction di... #xi=xi1;
                self.space_utils.UnaliasedAdd(xi, alpha, di)
                #updating residual... #ri=ri1
                #print "update r"
                self.space_utils.UnaliasedAdd(ri, (-1.0)*alpha, LHS_d)
                #here ri is already ri1
                self.builder_and_solver.calc_prod_precond_vec(ri, preconditioner, prec_ri)
                prod_ri1_Prec_ri1=self.space_utils.Dot(ri, prec_ri)
                #note, that since ri is updated, prec_ri means prec_ri1
                #di = self.update_search_direction(prec_ri, di, prod_ri1_Prec_ri1, prod_ri_Prec_ri);
                #updating the searcg direction
                #print "update search dir"
                self.space_utils.ScaleAndAdd(1.0, prec_ri,( prod_ri1_Prec_ri1/prod_ri_Prec_ri) , di)
                counter+=1

                if (self.builder_and_solver.ConvergenceCheck(ri, self.b, conv_criterion, counter, max_iterations)):            
                    print "CONVERGENCE ACHIEVED TO THE REQUIRED PRECISION"
                    print "counter", counter
                    #self.out_file.write(str(counter)+"\n")
                    self.builder_and_solver.ReturnDx(self.Dx, xi)
                    break
            #in the rest of non-linear iterations we dont need to add the "compressible" term to the tangent
        else:
            print "NOT ADDING GM-1D:........::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
            print "NOT ADDING GM-1D:........::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
            print "NOT ADDING GM-1D:........::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
                     
            for i in range(1, max_iterations):
                #Calculate: PrecInv_ri
                self.builder_and_solver.calc_prod_precond_vec(ri, preconditioner, prec_ri)
                prod_ri_Prec_ri=self.space_utils.Dot(ri, prec_ri)
                #computing first (standard part, that is A_d) of LHS_d
                self.space_utils.Mult(self.A, di, LHS_d)

                #now alpha
                alpha = self.calc_alpha(prod_ri_Prec_ri, LHS_d,di);
                #updating/resetting xi,new residual ri, new search direction di... #xi=xi1;
                self.space_utils.UnaliasedAdd(xi, alpha, di)
                #updating residual... #ri=ri1
                #print "update r"
                self.space_utils.UnaliasedAdd(ri, (-1.0)*alpha, LHS_d)
                #here ri is already ri1
                self.builder_and_solver.calc_prod_precond_vec(ri, preconditioner, prec_ri)
                prod_ri1_Prec_ri1=self.space_utils.Dot(ri, prec_ri)
                #note, that since ri is updated, prec_ri means prec_ri1
                #di = self.update_search_direction(prec_ri, di, prod_ri1_Prec_ri1, prod_ri_Prec_ri);
                #updating the searcg direction
                #print "update search dir"
                self.space_utils.ScaleAndAdd(1.0, prec_ri,( prod_ri1_Prec_ri1/prod_ri_Prec_ri) , di)
                counter+=1

                if (self.builder_and_solver.ConvergenceCheck(ri, self.b, conv_criterion, counter, max_iterations)):            
                    print "CONVERGENCE ACHIEVED TO THE REQUIRED PRECISION"
                    print "counter", counter
                    #self.out_file.write(str(counter)+"\n")
                    self.builder_and_solver.ReturnDx(self.Dx, xi)
                    break
                
                
    
    #######################################################################
    def PredictionStep(self,domain_size,UlfUtils):
        reform_dofs = True
        self.InitializeSolutionStep(reform_dofs);

        #perform prediction 
        self.Predict()

        if(self.MoveMeshFlag == True):
            self.scheme.MoveMesh(self.model_part.Nodes);
      
    #######################################################################
    def MoveMesh(self):
        self.scheme.MoveMesh(self.model_part.Nodes);
