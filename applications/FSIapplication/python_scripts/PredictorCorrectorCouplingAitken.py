#importing the Kratos Library
from Kratos import *
from KratosFSIApplication import *
from KratosIncompressibleFluidApplication import *
from KratosALEApplication import *
import math

def AddVariables(fluid_model_part,structure_model_part):
    import incompressible_fluid_solver
    incompressible_fluid_solver.AddVariables(fluid_model_part);
    structure_model_part.AddNodalSolutionStepVariable(RELAXED_DISP);
    print "variables for FractionalStepCoupling added correctly"

def AddDofs(fluid_model_part,structure_model_part):
    import incompressible_fluid_solver
    incompressible_fluid_solver.AddDofs(fluid_model_part);

class PredictorCorrectorCoupling:
    
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

        #setting fsi_convergence_toll
        self.fsi_convergence_toll = 0.001;
        self.fsi_absolute_toll = 1e-6;

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        pILUPrecond = ILU0Preconditioner()
        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
        self.pressure_linear_solver =  BICGSTABSolver(1e-5, 5000,pILUPrecond)

        #definition of the neighbour search strategy
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(fluid_model_part,number_of_avg_elems,number_of_avg_nodes)

        self.step = 0

        self.max_mu = 0.0
        
    def Initialize(self):
        (self.neighbour_search).Execute()
        
        self.fluid_solver = ResidualBasedFluidStrategy(self.fluid_model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)   

        (self.fluid_solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"

    def Solve_only_fluid(self):
        (self.fluid_solver).Solve()

    def Solve_only_structure(self):
        (self.structural_solver).Solve()
        
    def SolutionStep1(self):
        normDx = Array3(); normDx[0] = 0.00; normDx[1] = 0.00; normDx[2] = 0.00;
        is_converged = False
        iteration = 0

        while(	is_converged == False and iteration < self.max_vel_its  ): 
	    (self.fluid_solver).FractionalVelocityIteration(normDx);
            is_converged = (self.fluid_solver).ConvergenceCheck(normDx,self.vel_toll);
	    #print iteration,normDx
            iteration = iteration + 1
            
    def Solve(self):

        ###################################################################################
        ############################# viscous computation  ################################
        ###################################################################################

        #predicting the pressure on the structural side
        if( self.force_prediction_order != 0):
            structure_buffer_size = self.structure_model_part.GetBufferSize()
            (self.fsi_utils).StructuralPressurePrediction(POSITIVE_FACE_PRESSURE,self.interface_structure_nodes, structure_buffer_size ,self.force_prediction_order)
        else:
            (self.mapper).FluidToStructure_ScalarMap(PRESSURE,POSITIVE_FACE_PRESSURE)
        
        #solve the structure (prediction)
        print "performing Structural Prediction"
        (self.structural_solver).solver.perform_prediction = True
        (self.structural_solver).Solve()
        print "Structural Prediction: OK"

        #saving the relaxed and unrelaxed vars at this iteration -> makes a copy 
        #from the VariableListDataValueContainer to the DataValueContainer
        #note that the first iteration is unrelaxed
        (self.utilities).SaveVectorVar(DISPLACEMENT,DISPLACEMENT,self.interface_structure_nodes )

        (self.utilities).CopyVectorVar(DISPLACEMENT,RELAXED_DISP,self.interface_structure_nodes);
        (self.utilities).SaveVectorVar(RELAXED_DISP,RELAXED_DISP,self.interface_structure_nodes )

        
        ##map displacements to the structure
        (self.mapper).StructureToFluid_VectorMap(RELAXED_DISP,DISPLACEMENT)

        ##move the mesh
        (self.mesh_solver).Solve()
        
        ##set the fluid velocity at the interface to be equal to the corresponding mesh velocity
        (self.utilities).CopyVectorVar(MESH_VELOCITY,VELOCITY,self.interface_fluid_nodes);
        print "Displacement Map: OK"

        print "Mesh Movement: OK"

        ##initialization of the fluid solution step
        (self.fluid_solver).InitializeFractionalStep(self.step, self.time_order);
	(self.fluid_solver).InitializeProjections(self.step);
	(self.fluid_solver).AssignInitialStepValues();
	
        ###################################################################################
        ################### coupled pressure-structure computation ########################
        ###################################################################################
        fsi_is_converged = False
        iteration = 0

        self.mu = 0.0; #first step unaccelerated
        #self.mu = self.max_mu

        while(	fsi_is_converged == False and iteration < self.max_coupled_its  ): 

            ##solving the fractional step - taking in account the non linearity in the convection
            self.SolutionStep1()

            #pressure solution
            dp = (self.fluid_solver).SolveStep2();

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

            #check convergence in pressure
#            fsi_is_converged = (self.fsi_utils).CheckPressureConvergence(self.interface_fluid_nodes,self.fsi_convergence_toll);
##            p_norm = 0.00
##            dp_norm = 0.00
##            tmp = 0.0
##            for node in self.interface_fluid_nodes:
##                p = node.GetSolutionStepValue(PRESSURE)
##                tmp = p - node.GetSolutionStepValue(PRESSURE_OLD_IT)
##                p_norm += p*p
##                dp_norm += tmp*tmp;
##            p_norm = math.sqrt(p_norm);
##            dp_norm = math.sqrt(dp_norm);
##            ratio = math.sqrt(dp_norm/p_norm)
##            print "dp_norm = ",dp_norm
##            print "p_norm =", p_norm
##            print "fsi convergence ratio ==================== " , ratio
##            if(ratio < 0.01):
##                fsi_is_converged = True
##            else:
##                fsi_is_converged = False
##	    print "**********************************************","coupled it = ",iteration,"dp =",dp

            
            iteration = iteration + 1
            
                
            #######################################################################
            ##relax the structural solution before transferring
##            print "old_mu = ",self.mu
            if(iteration < 2):
                new_mu = self.mu
            else:
                new_mu = self.aitken_utils.ComputeAitkenFactor(self.interface_structure_nodes,self.mu);
            omega = 1.0 - new_mu;
            print "omega = ",omega
                
##            print "old_mu = ",self.mu
##            if(iteration <= 2):
##                omega = 1
##            else:
##                omega = self.aitken_utils.ComputeAitkenFactor(self.interface_structure_nodes,self.mu);
##                omega = -1.0*omega
                
   #         print "new_mu = ",new_mu



            #saving data for the next time step
            if(self.max_mu < new_mu):
                self.max_mu = new_mu
            self.mu = new_mu
            (self.utilities).SaveVectorVar(DISPLACEMENT,DISPLACEMENT,self.interface_structure_nodes )
            (self.utilities).SaveVectorVar(RELAXED_DISP,RELAXED_DISP,self.interface_structure_nodes )


            self.aitken_utils.ComputeRelaxedDisplacement(self.interface_structure_nodes,omega)


            ddisp_norm = 0.0
            deltadisp_norm = 0.0
            tmp = Vector(3)
            tmp1 = Vector(3)
            for node in self.interface_structure_nodes:
                d = node.GetSolutionStepValue(RELAXED_DISP)
                dold = node.GetSolutionStepValue(RELAXED_DISP,1)
                dold_it = node.GetValue(RELAXED_DISP)
##                print d,dold,dold_it
                ddisp_norm += (d[0]-dold_it[0])**2 + (d[1]-dold_it[1])**2 + (d[2]-dold_it[2])**2
                deltadisp_norm += (d[0]-dold[0])**2 + (d[1]-dold[1])**2 + (d[2]-dold[2])**2

            if(deltadisp_norm != 0.0):
                ratio = math.sqrt(ddisp_norm/deltadisp_norm)
            else:
                ratio = 1.0
            print "it=",iteration, " ratio = ",ratio, " ddisp_norm =",math.sqrt(ddisp_norm)/len(self.interface_fluid_nodes), " deltadisp_norm =",math.sqrt(deltadisp_norm)/len(self.interface_fluid_nodes)
            if(ratio < self.fsi_convergence_toll or math.sqrt(ddisp_norm)/len(self.interface_fluid_nodes)<self.fsi_absolute_toll):
                if(iteration <3):
                    fsi_is_converged = False
                else:
                    print "*************************************************"
                    print "**************** FSI CONVERGED ******************"
                    print "*************************************************"
                    print " "
                    fsi_is_converged = True
            else:
                fsi_is_converged = False    



            

            #######################################################################
            ########## TRANSFERRING STRUCTURAL SOLUTION TO THE FLUID ##############
            ## map displacements to the structure
            (self.mapper).StructureToFluid_VectorMap(RELAXED_DISP,DISPLACEMENT)

            ##move the mesh -- only on the interface NOT inside the domain
##            (self.mesh_solver).solver.MoveNodes();
            (self.mesh_solver).Solve()

            ## determine the corresponding velocity
            (self.utilities).CopyVectorVar(MESH_VELOCITY,VELOCITY,self.interface_fluid_nodes);
            
            #######################################################################
            ## finalizing fluid solution and getting ready for the next step ######
            ##update the projections
            (self.fluid_solver).ActOnLonelyNodes();
            (self.fluid_solver).SolveStep3();

            if(self.oss_stabilization == False)
                (self.utilities).SetToZero_VectorVar(PRESS_PROJ,self.fluid_model_part.Nodes)

            #correct the velocity field by applying the pressure
            (self.fluid_solver).SolveStep4();       
                
            if(fsi_is_converged == False ): 
                ## prepare the pressure for the next iteration (overwrites the pressure at the old step)
                (self.utilities).CopyVectorVar(VELOCITY,FRACT_VEL,(self.fluid_model_part).Nodes);
                (self.fluid_solver).SavePressureIteration()

        self.step = self.step  + 1
                
        







        
    
