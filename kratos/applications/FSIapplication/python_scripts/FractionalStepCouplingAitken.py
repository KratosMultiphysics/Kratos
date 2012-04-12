from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FSIApplication import *
from KratosMultiphysics.ALEApplication import *
CheckForPreviousImport()

import math

def AddVariables(fluid_model_part,structure_model_part):
    import fractional_step_solver
    fractional_step_solver.AddVariables(fluid_model_part);
    fluid_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    structure_model_part.AddNodalSolutionStepVariable(RELAXED_DISP);
    print "variables for FractionalStepCoupling added correctly"

def AddDofs(fluid_model_part,structure_model_part):
    import fractional_step_solver
    fractional_step_solver.AddDofs(fluid_model_part);

class FractionalStepCoupling:
    
    def __init__(self,fluid_model_part,structure_model_part,structural_solver,mesh_solver,mapper,domain_size):

        #saving solvers
        self.structural_solver = structural_solver
        self.mesh_solver = mesh_solver
        self.mapper = mapper
        self.domain_size = domain_size

        #saving model parts
        self.fluid_model_part = fluid_model_part
        self.structure_model_part = structure_model_part

        #creating utilities needed
        self.utilities = VariableUtils()
        self.fsi_utils = FSIUtils()
        self.aitken_utils = AitkenUtils()

        #initialize the list of interface nodes
        self.interface_structure_nodes = (self.utilities).SelectNodeList(IS_INTERFACE,1.0,structure_model_part.Nodes)
        self.interface_fluid_nodes = (self.utilities).SelectNodeList(IS_INTERFACE,1.0,fluid_model_part.Nodes)

        #settings for the fractional step fluid solver
        self.vel_toll = 0.001;
        self.press_toll = 0.001;
        self.max_vel_its = 4;
        self.max_press_its = 3;
        self.time_order = 2;
        self.CalculateReactions = False;
        self.ReformDofAtEachIteration = True; 
        self.CalculateNormDxFlag = True;
        self.laplacian_form = 1; #1 = laplacian, 2 = Discrete Laplacian
        self.predictor_corrector = False;
        self.echo_level = 0
        self.oss_stabilization = True

        #iterative coupling
        self.max_coupled_its = 20
        self.force_prediction_order = 2 #if 0 no pressure prediction ... using the value available
        self.incremental_structural_solution = True #setting to true the structural prediction is done just once and the solution is corrected between iterations
        self.complete_mesh_move_at_iterations = False
        self.switch_off_accelerator = False

        #setting fsi_convergence_toll
        self.fsi_convergence_toll = 0.001;
        self.fsi_absolute_toll = 1e-6;

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
##        pILUPrecond = ILU0Preconditioner()
        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
        self.pressure_linear_solver =  BICGSTABSolver(1e-5, 5000,pDiagPrecond)

        #definition of the neighbour search strategy
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(fluid_model_part,number_of_avg_elems,number_of_avg_nodes)

        self.step = 0

        self.max_mu = 0.0

        self.projections_are_initialized = False

        self.use_dt_in_stabilization = True
        self.fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, 1);
        
        self.TolFactor = 1e-3

    #if this is not used the beta is taken as 1
##    def CalculateFsiBeta(self,hs,hf,young,poisson):
##        beta = ...
##        (self.fluid_solver).SetFSIBeta(beta)
        
        
    
        
    def Initialize(self):
        (self.neighbour_search).Execute()
        
##        self.fluid_solver = ResidualBasedFluidStrategy(self.fluid_model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)   
        solver_configuration = FractionalStepConfiguration(self.fluid_model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.domain_size,self.laplacian_form )
        self.fluid_solver = FractionalStepStrategy( self.fluid_model_part, solver_configuration, self.ReformDofAtEachIteration, self.vel_toll, self.press_toll, self.max_vel_its, self.max_press_its, self.time_order, self.domain_size, self.predictor_corrector)

        (self.fluid_solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"



    def Solve_only_structure(self):
        (self.structural_solver).Solve()

    def Solve_structure_and_mesh(self):
        (self.structural_solver).Solve()
        (self.mapper).StructureToFluid_VectorMap(DISPLACEMENT,DISPLACEMENT)
        (self.mesh_solver).Solve()
        
    def SolutionStep1(self):
        normDx = Array3(); normDx[0] = 0.00; normDx[1] = 0.00; normDx[2] = 0.00;
        is_converged = False
        iteration = 0

        while(	is_converged == False and iteration < self.max_vel_its  ): 
	    normDx = (self.fluid_solver).FractionalVelocityIteration();
            is_converged = (self.fluid_solver).ConvergenceCheck(normDx,self.vel_toll);
	    #print iteration,normDx
            iteration = iteration + 1

    def Solve_only_fluid(self):
##        (self.fluid_solver).Solve()

        ##initialization of the fluid solution step
        (self.fluid_solver).InitializeFractionalStep(self.step, self.time_order);
	(self.fluid_solver).InitializeProjections(self.step,self.projections_are_initialized);
        self.projections_are_initialized = True
	(self.fluid_solver).AssignInitialStepValues();       

        ##solving the fractional step - taking in account the non linearity in the convection
        self.SolutionStep1()

        for i in range(0,4):
            dp = (self.fluid_solver).SolveStep2();
            
            (self.fluid_solver).SolveStep3();

            (self.fluid_solver).SolveStep4();       
                
            ## prepare the pressure for the next iteration (overwrites the pressure at the old step)
            (self.utilities).CopyVectorVar(VELOCITY,FRACT_VEL,(self.fluid_model_part).Nodes);
            (self.fluid_solver).SavePressureIteration()


    def SavePressure(self,temp_vect):
        i = 0;
        for node in (self.fluid_model_part).Nodes:
            temp_vect[i] = node.GetSolutionStepValue(PRESSURE)
            i = i + 1

        return temp_vect

    def WritePressure(self,relaxed_pressure):
        i = 0;
        for node in self.fluid_model_part.Nodes:
            node.SetSolutionStepValue(PRESSURE,0,relaxed_pressure[i])
            i = i + 1
     

    def Relax(self,pnew,pold,resnew,resold,w):
        num = 0.0
        denom = 0.0
        for i in range(0,len(pnew)):
            resnew[i] = pnew[i] - pold[i] #up to this point the prelaxed is still unrelaxed
            pnew[i] = w * pnew[i] + (1.0-w)*pold[i]
            num = num + resold[i]*(resnew[i]-resold[i])
            denom = denom + (resnew[i]-resold[i])*(resnew[i]-resold[i])

            #prepare for next step
            pold[i] = pnew[i];
            resold[i] = resnew[i]

        print "num=",num
        print "denom=",denom
        wnew = -w * num / denom
        print "wnew",wnew
        return [pnew,pold,resnew,resold,wnew]
        
            
    def Solve(self):

        ###################################################################################
        ############################# viscous computation  ################################
        ###################################################################################

        #predicting the pressure on the structural side
        if( self.force_prediction_order != 0):
            structure_buffer_size = self.structure_model_part.GetBufferSize()
            (self.fsi_utils).StructuralPressurePrediction(POSITIVE_FACE_PRESSURE,self.interface_structure_nodes, structure_buffer_size ,self.force_prediction_order)
        else:
            (self.fluid_solver).Solve()
            (self.mapper).FluidToStructure_ScalarMap(PRESSURE,POSITIVE_FACE_PRESSURE)
        
        #solve the structure (prediction)
        print "performing Structural Prediction"
        (self.structural_solver).solver.perform_prediction = True
        (self.structural_solver).Solve()
        print "Structural Prediction: OK"

        #saving the relaxed and unrelaxed vars at this iteration -> makes a copy 
        #from the VariableListDataValueContainer to the DataValueContainer
        #note that the first iteration is unrelaxed
##        (self.utilities).SaveVectorVar(DISPLACEMENT,DISPLACEMENT,self.interface_structure_nodes )
##
##        (self.utilities).CopyVectorVar(DISPLACEMENT,RELAXED_DISP,self.interface_structure_nodes);
##        (self.utilities).SaveVectorVar(RELAXED_DISP,RELAXED_DISP,self.interface_structure_nodes )

        
        ##map displacements to the structure
        (self.mapper).StructureToFluid_VectorMap(DISPLACEMENT,DISPLACEMENT)
##        (self.mapper).StructureToFluid_VectorMap(RELAXED_DISP,DISPLACEMENT)

        ##move the mesh
        (self.mesh_solver).Solve()
        
        ##set the fluid velocity at the interface to be equal to the corresponding mesh velocity
        (self.utilities).CopyVectorVar(MESH_VELOCITY,VELOCITY,self.interface_fluid_nodes);
        print "Displacement Map: OK"

        print "Mesh Movement: OK"

        ##initialization of the fluid solution step
        (self.fluid_solver).InitializeFractionalStep(self.step, self.time_order);
	(self.fluid_solver).InitializeProjections(self.step,self.projections_are_initialized);
	self.projections_are_initialized = True
	(self.fluid_solver).AssignInitialStepValues();       

        ##solving the fractional step - taking in account the non linearity in the convection
        self.SolutionStep1()

        ###################################################################################
        ################### coupled pressure-structure computation ########################
        ###################################################################################
        fsi_is_converged = False
        iteration = 0

        self.mu = 0.0; #first step unaccelerated
        #self.mu = self.max_mu

        initial_residual = 1.0

##        pold = Vector(len(self.fluid_model_part.Nodes))
##        pnew = Vector(len(self.fluid_model_part.Nodes))
##        resnew = Vector(len(self.fluid_model_part.Nodes))
##        resold = Vector(len(self.fluid_model_part.Nodes))
##
##        pold = self.SavePressure(pold)
##        w = 1.0
##        wmax = 2.0
##        wmin = 0.01
##
##        for i in range(0,len(pnew)):
##            resold[i] = 0.0

        full_tolerance = 1e-3;
        max_tolerance = 0.1;
        gamma = 0.9
        MaxDecreaseFactor = 0.1
        
        old_residual = 1.0;
        #self.pressure_linear_solver.SetTolerance(full_tolerance)
        self.TolFactor = full_tolerance
        print self.pressure_linear_solver

        while(	fsi_is_converged == False and iteration < self.max_coupled_its  ): 
            #pressure solution
            dp = (self.fluid_solver).SolveStep2();
            residual = (self.fluid_solver).GetStageResidualNorm(4);
            print "residual ============================================ ",residual
            if(iteration == 0):
                initial_residual = residual

            #linear tolerance reduction
            if(self.switch_off_accelerator == False):
                if(iteration > 0):
                    CandidateFactor = gamma*(residual**2)/(old_residual**2);
                    CandidateFactor_LimitedDecrease = gamma*self.TolFactor*self.TolFactor;

                    if (CandidateFactor_LimitedDecrease < MaxDecreaseFactor):
                        if(CandidateFactor < max_tolerance):
                            self.TolFactor = CandidateFactor
                        else:
                            self.TolFactor = max_tolerance
                    else:
                        temp = 0.0
                        if(CandidateFactor > CandidateFactor_LimitedDecrease):
                            temp = CandidateFactor
                        else:
                            temp = CandidateFactor_LimitedDecrease

                        if(temp < max_tolerance):
                            self.TolFactor  = temp;
                        else:
                            self.TolFactor  = max_tolerance;

                    if(self.TolFactor < full_tolerance):
                        self.TolFactor = full_tolerance;

                    self.pressure_linear_solver.SetTolerance(self.TolFactor)
                    print self.pressure_linear_solver
                

##            pnew = self.SavePressure(pnew)
##            [pnew,pold,resnew,resold,w] = self.Relax(pnew,pold,resnew,resold,w)
##            if(iteration == 0):
##                w = 1.0
##
##            if(w > wmax):
##                w=wmax
##            if(w < wmin):
##                w = wmin
##            print "************************* ------------------------- ************** ", w
##
##            self.WritePressure(pnew)
            

            


            

            #######################################################################
            ################## UPDATING THE PRESSURE SOLUTION #####################
            ##map pressure to the structure
            (self.mapper).FluidToStructure_ScalarMap(PRESSURE,POSITIVE_FACE_PRESSURE)
#            (self.mapper).FluidToStructure_ScalarMap(PRESSURE,NEGATIVE_FACE_PRESSURE)
            print "Pressure Map OK"
            
            #correct the structural solution
            if( self.incremental_structural_solution == True):
                (self.structural_solver).solver.perform_prediction = False
            (self.structural_solver).Solve()
            print "Structural Correction: OK"

            (self.mapper).StructureToFluid_VectorMap(DISPLACEMENT,DISPLACEMENT)

            ###checking convergence ##3
            ratio = residual / initial_residual
            print "it=",iteration, " ratio = ",ratio, "absolute residual norm",residual
            if(ratio < self.fsi_convergence_toll or residual<self.fsi_absolute_toll):
                    print "*************************************************"
                    print "**************** FSI CONVERGED ******************"
                    print "*************************************************"
                    print " "
                    fsi_is_converged = True
            else:
                fsi_is_converged = False
##
##                
##
####            ddisp_norm = 0.0
####            deltadisp_norm = 0.0
####            tmp = Vector(3)
####            tmp1 = Vector(3)
####            for node in self.interface_structure_nodes:
####                d = node.GetSolutionStepValue(RELAXED_DISP)
####                dold = node.GetSolutionStepValue(RELAXED_DISP,1)
####                dold_it = node.GetValue(RELAXED_DISP)
######                print d,dold,dold_it
####                ddisp_norm += (d[0]-dold_it[0])**2 + (d[1]-dold_it[1])**2 + (d[2]-dold_it[2])**2
####                deltadisp_norm += (d[0]-dold[0])**2 + (d[1]-dold[1])**2 + (d[2]-dold[2])**2
####
####            if(deltadisp_norm != 0.0):
####                ratio = math.sqrt(ddisp_norm/deltadisp_norm)
####            else:
####                ratio = 1.0
####            print "it=",iteration, " ratio = ",ratio, " ddisp_norm =",math.sqrt(ddisp_norm)/len(self.interface_fluid_nodes), " deltadisp_norm =",math.sqrt(deltadisp_norm)/len(self.interface_fluid_nodes)
####            if(ratio < self.fsi_convergence_toll or math.sqrt(ddisp_norm)/len(self.interface_fluid_nodes)<self.fsi_absolute_toll):
####                if(iteration <3):
####                    fsi_is_converged = False
####                else:
####                    print "*************************************************"
####                    print "**************** FSI CONVERGED ******************"
####                    print "*************************************************"
####                    print " "
####                    fsi_is_converged = True
####            else:
####                fsi_is_converged = False
##            
####            print iteration
####            if(iteration < self.max_coupled_its):
####                if(iteration <3):
####                    fsi_is_converged = False
####                else:
####                    
######                    fsi_is_converged = (self.fsi_utils).CheckPressureConvergence(self.fluid_model_part.Nodes,self.fsi_convergence_toll);
######                    fsi_is_converged = (self.fsi_utils).CheckPressureConvergence(self.interface_fluid_nodes,self.fsi_convergence_toll);
####                    pressure_toll = self.fsi_convergence_toll
####                    fsi_is_converged = (self.fsi_utils).CheckPressureConvergence(self.interface_fluid_nodes,pressure_toll);
####                    if(fsi_is_converged == True):
####                        print "*************************************************"
####                        print "**************** FSI CONVERGED ******************"
####                        print "*************************************************"
####                        print "pipoo"
####                        print " "
####            else:
####                fsi_is_converged = True
####                print "*************************************************"
####                print "** ATTENTION: FSI NOT CONVERGED *****************"
####                print "*************************************************"
####                print " "
##
##
##            
##
##            #######################################################################
##            ########## TRANSFERRING STRUCTURAL SOLUTION TO THE FLUID ##############
##            ## map displacements to the structure
##            (self.mapper).StructureToFluid_VectorMap(RELAXED_DISP,DISPLACEMENT)

            ##move the mesh -- only on the interface NOT inside the domain
            if( self.complete_mesh_move_at_iterations == False):
                (self.mesh_solver).MoveNodes();
            else:
                (self.mesh_solver).Solve()

            ## determine the corresponding velocity
            (self.utilities).CopyVectorVar(MESH_VELOCITY,VELOCITY,self.interface_fluid_nodes);
            
            #######################################################################
            ## finalizing fluid solution and getting ready for the next step ######
            ##update the projections
            (self.fluid_solver).ActOnLonelyNodes();
            if(iteration == (self.max_coupled_its-1) or fsi_is_converged == True ): 
                (self.fluid_solver).SolveStep3();
##            (self.fluid_solver).SolveStep3();

            if(self.oss_stabilization == False):
                (self.utilities).SetToZero_VectorVar(PRESS_PROJ,self.fluid_model_part.Nodes)

            #correct the velocity field by applying the pressure
            (self.fluid_solver).SolveStep4();       
            (self.fluid_solver).SavePressureIteration()
            
            if(fsi_is_converged == False and iteration < self.max_coupled_its ): 
                ## prepare the pressure for the next iteration (overwrites the pressure at the old step)
                (self.utilities).CopyVectorVar(VELOCITY,FRACT_VEL,(self.fluid_model_part).Nodes);
                

            iteration = iteration + 1

            old_residual = residual;

            

        self.step = self.step  + 1

        return iteration
                
        







        
    
