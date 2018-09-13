#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()
import os,time

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
    def __init__( self, model_part, time_scheme, linear_solver, convergence_criteria, CalculateReactionsFlag, ReformDofSetAtEachStep, MoveMeshFlag, Parameters, space_utils, builder_and_solver ):
        #save the input parameters
        self.model_part = model_part
        self.scheme = time_scheme
        self.linear_solver = linear_solver
        self.convergence_criteria = convergence_criteria
        self.CalculateReactionsFlag = CalculateReactionsFlag
        self.ReformDofSetAtEachStep = ReformDofSetAtEachStep
        self.MoveMeshFlag = MoveMeshFlag
        self.Parameters = Parameters
        self.PerformContactAnalysis = self.Parameters['perform_contact_analysis_flag']
        self.PrintSparsity = self.Parameters['print_sparsity_info_flag']
        self.space_utils = space_utils
        #contact utility
        self.cu = ContactUtility( 3 )
        #default values for some variables
        self.max_iter = 30
        self.echo_level = 1
        self.builder_and_solver = builder_and_solver
        self.dof_util = DofUtility()

        #local matrices and vectors
        self.pA = self.space_utils.CreateEmptyMatrixPointer()
        self.pDx = self.space_utils.CreateEmptyVectorPointer()
        self.pb = self.space_utils.CreateEmptyVectorPointer()

        self.A = (self.pA).GetReference()
        self.Dx = (self.pDx).GetReference()
        self.b = (self.pb).GetReference()
        ##local matrices and vectors
        #self.A = CompressedMatrix()
        #self.Dx = Vector()
        #self.b = Vector()
        
        #initialize flags
        self.SolutionStepIsInitialized = False		
        self.InitializeWasPerformed = False
        self.StiffnessMatrixIsBuilt = False
        #provide settings to the builder and solver
        (self.builder_and_solver).SetCalculateReactionsFlag(self.CalculateReactionsFlag);
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep);
        
        self.solveCounter = 0; #hbui added this variable
        
        self.system_reorderer = Process()
#        self.system_reorderer = SystemRCMReordererProcess(self.model_part)
#        self.system_reorderer = SystemBoostRCMReordererProcess(self.model_part)
#        self.system_reorderer = SystemMetisReordererProcess(self.model_part)
#        self.system_reorderer = SystemAMDReordererProcess(self.model_part)

        self.attached_processes = []

    #######################################################################
    def Initialize(self):
        if(self.scheme.SchemeIsInitialized() == False):
            self.scheme.Initialize(self.model_part)
        if (self.scheme.ElementsAreInitialized() == False): 
            self.scheme.InitializeElements(self.model_part)
        for proc in self.attached_processes:
            proc.ExecuteInitialize()

    #######################################################################
    def SolveLagrange( self ):
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
        last_real_node = len(self.model_part.Nodes)
        originalPosition = self.cu.SetUpContactConditionsLagrangeTying(self.model_part )
        self.PerformNewtonRaphsonIteration()
        self.cu.CleanLagrangeTying( self.model_part, originalPosition, last_real_node )
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)
        #finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)
        #clear if needed - deallocates memory 
        if(self.ReformDofSetAtEachStep == True):
            self.Clear();

    #######################################################################
    def Solve(self):
        #print self.model_part
        ## - storing original condition size before adding virtual conditions.
        ## - performing contact search
        ## - creating virtual link conditions for the assembling
        self.solveCounter = self.solveCounter + 1
        if( self.PerformContactAnalysis == False ):
            self.PerformNewtonRaphsonIteration()
            #finalize the solution step
            self.FinalizeSolutionStep(self.CalculateReactionsFlag)
            #clear if needed - deallocates memory 
            if(self.ReformDofSetAtEachStep == True):
                self.Clear();
            return
        print "setting up contact conditions"
        originalPosition =  self.cu.SetUpContactConditions(self.model_part, self.Parameters['penalty'], self.Parameters['frictionpenalty'], self.Parameters['contact_double_check_flag'] )
        uzawaConverged = False
        ##  First step: reform DOF set and check if uzawa iteration is necessary
        self.PerformNewtonRaphsonIteration()
        self.cu.Update( self.model_part, originalPosition, self.Parameters['friction'], self.Parameters['contact_ramp_penalties_flag'], self.Parameters['rampcriterion'], self.Parameters['fricrampcriterion'], self.Parameters['rampfactor'], self.Parameters['fricrampfactor'], self.Parameters['maxpenalty'], self.Parameters['fricmaxpenalty']  )
        if( self.cu.IsConverged( self.model_part, 0,  originalPosition, self.Parameters['friction'] ) == True ):
            uzawaConverged = True
            (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)
            self.cu.Clean( self.model_part, originalPosition );
            #finalize the solution step
            self.FinalizeSolutionStep(self.CalculateReactionsFlag)
            #clear if needed - deallocates memory 
            if(self.ReformDofSetAtEachStep == True):
                self.Clear()
            return
        ## beginning of UZAWA loop
        (self.builder_and_solver).SetReshapeMatrixFlag(False)
        #props = self.model_part.Properties[1]
        for uzawaStep in range(1, self.Parameters['maxuzawa'] ):
            print "I am inside the uzawa loop, iteration no. " + str(uzawaStep)
            ## solving the standard newton-raphson iteration 
            self.PerformNewtonRaphsonIteration()
            ## updating the lagrange multipliers
            self.cu.Update( self.model_part, originalPosition, self.Parameters['friction'], self.Parameters['contact_ramp_penalties_flag'], self.Parameters['rampcriterion'], self.Parameters['fricrampcriterion'], self.Parameters['rampfactor'], self.Parameters['fricrampfactor'], self.Parameters['maxpenalty'], self.Parameters['fricmaxpenalty']  )
            ## checking convergence
            if( self.cu.IsConverged( self.model_part, uzawaStep, originalPosition, self.Parameters['friction'] ) == True ):
                uzawaConverged = True
                break
        if(self.Parameters['maxuzawa'] == 1):
            print "Congratulations. Newton-Raphson loop has converged."
        else:
            if( uzawaConverged == False ):
                if('stop_Uzawa_if_not_converge' in self.Parameters):
                    if(self.Parameters['stop_Uzawa_if_not_converge'] == True):
                        sys.exit("Stop. Uzawa algorithm failed to converge at time step " + str(self.model_part.ProcessInfo[TIME]) + ", uzawaStep = " + str(uzawaStep) + ", maxuzawa = " + str(self.Parameters['maxuzawa']))
                    else:
                        print('Uzawa algorithm failed to converge. However, the iteration will still be proceeded' + ", uzawaStep = " + str(uzawa) + ", maxuzawa = " + str(self.Parameters['maxuzawa']))
                else:
                    print("Stop. Uzawa algorithm failed to converge at time step " + str(self.model_part.ProcessInfo[TIME]) + ", uzawaStep = " + str(uzawaStep) + ", maxuzawa = " + str(self.Parameters['maxuzawa']))
                    print("However, I still want to proceed")
            else:
                print 'Congratulations. Uzawa loop has converged.'
        ### end of UZAWA loop
        ### cleaning up the conditions
        self.cu.Clean( self.model_part, originalPosition )
        (self.builder_and_solver).SetReshapeMatrixFlag(self.ReformDofSetAtEachStep)
        #finalize the solution step
        self.FinalizeSolutionStep(self.CalculateReactionsFlag)
        #clear if needed - deallocates memory 
        if(self.ReformDofSetAtEachStep == True):
            self.Clear();
        
    def PerformNewtonRaphsonIteration( self ):
        print("time = " + str(self.model_part.ProcessInfo[TIME]))
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
            self.InitializeSolutionStep()
            self.SolutionStepIsInitialized = True
        #perform prediction 
        self.Predict()

        #execute iteration - first iteration is ALWAYS executed
        calculate_norm = False
        self.iterationCounter = 0 #hbui added this variable
        self.iterationCounter = self.iterationCounter + 1
        normDx = self.ExecuteIteration(self.echo_level,self.MoveMeshFlag,calculate_norm)

        original_penalty = 0.0
        if( self.PerformContactAnalysis == True ):
            for cond in self.model_part.Conditions:
                if( cond.GetValue( IS_CONTACT_SLAVE ) ):
                    original_penalty = cond.GetValue( PENALTY )[0]
                    break

        #non linear loop
        converged = False
        it = 0
        while(it < self.max_iter and converged == False):
            #increase penalty...
            if( self.PerformContactAnalysis == True ):
                for cond in self.model_part.Conditions:
                    if( cond.GetValue( IS_CONTACT_SLAVE ) ):
                        penalty = cond.GetValue( PENALTY )
                        for i in range(0,len(penalty)):
                            penalty[i] = 1.0*penalty[i]
                        cond.SetValue( PENALTY, penalty )
            #end of increase penalty
            #verify convergence
            converged = self.convergence_criteria.PreCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

            #calculate iteration
            # - system is built and solved
            # - database is updated depending on the solution
            # - nodal coordinates are updated if required
            self.iterationCounter = self.iterationCounter + 1
            normDx = self.ExecuteIteration(self.echo_level,self.MoveMeshFlag,calculate_norm)

            #verify convergence
            converged = self.convergence_criteria.PostCriteria(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b)

            #update iteration count
            it = it + 1
        
        if( it == self.max_iter and converged == False):
            print("Iteration did not converge at time step " + str(self.model_part.ProcessInfo[TIME]))
            if('stop_Newton_Raphson_if_not_converge' in self.Parameters):
                if(self.Parameters['stop_Newton_Raphson_if_not_converge'] == True):
                    sys.exit("Sorry, my boss does not allow me to continue. The time step did not converge at time step " + str(self.model_part.ProcessInfo[TIME]) + ", it = " + str(it) + ", max_iter = " + str(self.max_iter))
                else:
                    print('However, the iteration will still be proceeded' + ", it = " + str(it) + ", max_iter = " + str(self.max_iter))
            else:
                sys.exit("Sorry, my boss does not allow me to continue. The time step did not converge at time step " + str(self.model_part.ProcessInfo[TIME]) + ", it = " + str(it) + ", max_iter = " + str(self.max_iter))
        print("PerformNewtonRaphsonIteration converged after " + str(it) + " steps")
        if( self.PerformContactAnalysis == True ):
            for cond in self.model_part.Conditions:
                if( cond.GetValue( IS_CONTACT_SLAVE ) ):
                    penalty = cond.GetValue( PENALTY )
                    for i in range(0,len(penalty)):
                        penalty[i] = original_penalty
                    cond.SetValue( PENALTY, penalty )

    #######################################################################
    #######################################################################

    #######################################################################
    def Predict(self):
        self.scheme.Predict(self.model_part,self.builder_and_solver.GetDofSet(),self.A,self.Dx,self.b);

    #######################################################################
    def InitializeSolutionStep(self):
        if(self.builder_and_solver.GetDofSetIsInitializedFlag() == False or self.ReformDofSetAtEachStep == True):
            #initialize the list of degrees of freedom to be used 
            self.builder_and_solver.SetUpDofSet(self.scheme,self.model_part);
            #reorder the list of degrees of freedom to identify fixity and system size	  			
            self.builder_and_solver.SetUpSystem(self.model_part)
            #reorder the system dof id
            self.system_reorderer.Execute()
            #allocate memory for the system and preallocate the structure of the matrix
            self.builder_and_solver.ResizeAndInitializeVectors(self.scheme, self.pA,self.pDx,self.pb,self.model_part)
            #updating references
            self.A = (self.pA).GetReference()
            self.Dx = (self.pDx).GetReference()
            self.b = (self.pb).GetReference()
        if(self.SolutionStepIsInitialized == False):
            self.builder_and_solver.InitializeSolutionStep(self.model_part,self.A,self.Dx,self.b)
            self.scheme.InitializeSolutionStep(self.model_part,self.A,self.Dx,self.b)
        for proc in self.attached_processes:
            proc.ExecuteInitializeSolutionStep()

    #######################################################################
    def ExecuteIteration(self,echo_level,MoveMeshFlag,CalculateNormDxFlag):
        #reset system matrices and vectors prior to rebuild
        self.space_utils.SetToZeroMatrix(self.A)
        self.space_utils.SetToZeroVector(self.Dx)			
        self.space_utils.SetToZeroVector(self.b)

        self.scheme.InitializeNonLinIteration(self.model_part,self.A,self.Dx,self.b)

        #provide data for the preconditioner and linear solver
#        self.linear_solver.ProvideAdditionalData(self.A,self.Dx,self.b,self.builder_and_solver.GetDofSet(),self.model_part)

        #build and solve the problem
        if(self.Parameters['decouple_build_and_solve'] == False):
            self.builder_and_solver.BuildAndSolve(self.scheme,self.model_part,self.A,self.Dx,self.b)
            self.dof_util.ListDofs(self.builder_and_solver.GetDofSet(),self.builder_and_solver.GetEquationSystemSize())
        else:
            self.builder_and_solver.Build(self.scheme,self.model_part,self.A,self.b)
            self.dof_util.ListDofs(self.builder_and_solver.GetDofSet(),self.builder_and_solver.GetEquationSystemSize())
            self.linear_solver.Solve(self.A,self.Dx,self.b)

#        diagAstr = ""
#        for i in range(0, self.A.Size1()):
#            diagAstr = diagAstr + ", " + str(self.A[(i, i)])
#        print("diagonal A:" + diagAstr)

        #full output if needed
        if( self.PrintSparsity ):
            #hbui edited
            #self.PlotSparsityScheme( self.A )
#            wr = UblasMatrixIO()
#            wr.WriteHB(self.A, self.b, "matrix" + str(self.solveCounter) + "." + str(self.iterationCounter) + ".hb.dat")
#            self.space_utils.WriteMatrixMarketMatrix("matrix" + str(self.solveCounter) + "." + str(self.iterationCounter) + ".mm",self.A,False)
            petsc_utils.DumpUblasCompressedMatrixVector("tempAb", self.A, self.b, False)
            
        if(echo_level >= 3):
            print "SystemMatrix = ", self.A 
        #printA = []
        #printdx = []
        #printb = []
        #for i in range(0,len(self.Dx)):
        #    if( abs(self.Dx[i]) < 1.0e-10 ):
         #       printdx.append(0.0)
         #   else:
         #       printdx.append(self.Dx[i])
         #   if( abs(self.b[i]) < 1.0e-6 ):
         #       printb.append(0.0)
         #   else:
         #       printb.append(self.b[i])
         #   row = []
         #   for j in range(0,len(self.Dx)):
         #       if( abs(self.A[(i,j)]) < 1.0 ):
         #           row.append( 0.0 )
         #       else:
         #           row.append(self.A[(i,j)])
         #   printA.append(row)
            print "solution obtained = ", self.Dx
            #formatted_printdx = [ '%.6f' % elem for elem in printdx ]
            #print formatted_printdx
            #formatted_printb = [ '%.4f' % elem for elem in printb ]
            print "RHS = ", self.b
        #print formatted_printb
        #print "Matrix: "
        #for i in range(0,len(self.Dx)):
        #    formatted_printA = [ '%.1f' % elem for elem in printA[i] ]
        #    print(formatted_printA)
        self.AnalyseSystemMatrix(self.A)
            
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

        for proc in self.attached_processes:
            proc.ExecuteFinalizeSolutionStep()

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
        self.scheme.Clear()

        for proc in self.attached_processes:
            proc.ExecuteFinalize()

    #######################################################################   
    def SetEchoLevel(self,level):
        self.echo_level = level
        self.builder_and_solver.SetEchoLevel(level)

    #######################################################################   
    def AnalyseSystemMatrix(self,  A):
        max = 0.0
        for i in range(0,  A.Size1()):
           if( abs(A[(i, i)]) > max ):
               max = A[(i, i)]
               
#        nonzeros = 0
#        for i in range(0,  A.Size1()):
#            for j in range(0,  A.Size2()):
#                if( abs(A[(i, j)]) > 1e-16 ):
#                    nonzeros = nonzeros + 1
                    
        print("#############################")
        print("Number of rows: " +str(A.Size1()) )
        print("Number of columns: " +str(A.Size2()) )
#        print("Number of entries: " +str(nonzeros) )
        print("Max in Diagonal: " +str(max) )
        print("#############################")

    #######################################################################   
    def PlotSparsityScheme(self, A):
        try:
            import Gnuplot, Gnuplot.PlotItems, Gnuplot.funcutils
        except ImportError:
            # kludge in case Gnuplot hasn't been installed as a module yet:
            import __init__
            Gnuplot = __init__
            import PlotItems
            Gnuplot.PlotItems = PlotItems
            import funcutils
            Gnuplot.funcutils = funcutils
        print("gnuplot-python imported")
        g = Gnuplot.Gnuplot(debug=1)
        g.clear()
        #g.plot(Gnuplot.Func('sin(x)'))
        #self.wait('hit Return to continue')
        file = open("matrix.dat",'w')
        for i in range(0, A.Size1()):
            for j in range(0, A.Size2()):
                tmp = A[(i,j)]
                if( (tmp > 1.0e-9) or (tmp < -1.0e-9) ):
                   #file.write( str(tmp) +"\t" )
                   file.write( "1.0 " )
                else:
                   file.write("0.0 ")
            file.write("\n")
        file.close()
        g("set term postscript")
        g("set size square")
        g("set output 'matrix.ps'")
        g("set zrange [0.5:1.5]")
        g("set pm3d map")
        g("splot 'matrix.dat' matrix with dots")
        

    def wait(self,str=None, prompt='Press return to show results...\n'):
        if str is not None:
            print(str)
        raw_input(prompt)

        
