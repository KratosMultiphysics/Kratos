from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.FSIApplication import *
from KratosMultiphysics.ALEApplication import *
CheckForPreviousImport()

import math

def AddVariables(fluid_model_part,structure_model_part):
    import incompressible_fluid_solver
    incompressible_fluid_solver.AddVariables(fluid_model_part);
    print "variables for FractionalStepCoupling added correctly"

def AddDofs(fluid_model_part,structure_model_part):
    import incompressible_fluid_solver
    incompressible_fluid_solver.AddDofs(fluid_model_part);

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



        #initialize the list of interface nodes
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
        self.laplacian_form = 2; #1 = laplacian, 2 = Discrete Laplacian
        self.predictor_corrector = False;
        self.echo_level = 0

        #iterative coupling
        self.max_coupled_its = 20

        #setting fsi_convergence_toll
        self.fsi_convergence_toll = 5 * self.press_toll;

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        pILUPrecond = ILU0Preconditioner()
        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pILUPrecond)

        #definition of the neighbour search strategy
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(fluid_model_part,number_of_avg_elems,number_of_avg_nodes)

        self.step = 0
        
    def Initialize(self):
        (self.neighbour_search).Execute()
        
        self.fluid_solver = ResidualBasedFluidStrategy(self.fluid_model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)   

        (self.fluid_solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"

    def Solve_only_fluid(self):
        (self.fluid_solver).Solve()

    def Solve_only_structure(self):
        (self.fluid_solver).Solve()
        
    def SolutionStep1(self):
        normDx = Array3(); normDx[0] = 0.00; normDx[1] = 0.00; normDx[2] = 0.00;
        is_converged = False
        iteration = 0

        while(	is_converged == False and iteration < self.max_vel_its  ): 
	    (self.fluid_solver).FractionalVelocityIteration(normDx);
            is_converged = (self.fluid_solver).ConvergenceCheck(normDx,self.vel_toll);
	    print iteration,normDx
            iteration = iteration + 1
            
    def Solve(self):

        ###################################################################################
        ############################# viscous computation  ################################
        ###################################################################################
#here i should set to zero the acceleration in the DataValueContainer       
        #solve the structure (prediction)
        (self.structural_solver).Solve()
        print "Structural Prediction: OK"

        #saving the acceleration at this iteration -> makes a copy of the acceleration in the
        #from the VariableListDataValueContainer to the DataValueContainer
        (self.utilities).SaveVectorVar(ACCELERATION,ACCELERATION,(self.structure_model_part).Nodes )

        
        ##map displacements to the structure
        (self.mapper).StructureToFluid_VectorMap(DISPLACEMENT,DISPLACEMENT)
        ##set the fluid velocity at the interface to be equal to the corresponding mesh velocity
        (self.utilities).CopyVectorVar(MESH_VELOCITY,VELOCITY,self.interface_fluid_nodes);
        print "Displacement Map: OK"

        ##move the mesh
        (self.mesh_solver).Solve()
        print "Mesh Movement: OK"

        ##initialization of the fluid solution step
        (self.fluid_solver).InitializeFractionalStep(self.step, self.time_order);
	(self.fluid_solver).InitializeProjections(self.step);
	(self.fluid_solver).AssignInitialStepValues();       

        ##solving the fractional step - taking in account the non linearity in the convection
        self.SolutionStep1()

        ###################################################################################
        ################### coupled pressure-structure computation ########################
        ###################################################################################
        fsi_is_converged = False
        iteration = 0



        while(	fsi_is_converged == False and iteration < self.max_coupled_its  ): 
            #pressure solution
            dp = (self.fluid_solver).SolveStep2();

            #######################################################################
            ################## UPDATING THE PRESSURE SOLUTION #####################
            ##map pressure to the structure
            (self.mapper).FluidToStructure_ScalarMap(PRESSURE,POSITIVE_FACE_PRESSURE)
#            (self.mapper).FluidToStructure_ScalarMap(PRESSURE,NEGATIVE_FACE_PRESSURE)
            print "Pressure Map OK"
            
            #correct the structural solution
            (self.structural_solver).Solve()
            print "Structural Correction: OK"            

            #check convergence in pressure
            fsi_is_converged = (self.fsi_utils).CheckPressureConvergence(self.interface_fluid_nodes,self.fsi_convergence_toll);
##            p_norm = 0.00
##            dp_norm = 0.00
##            tmp
##            for node in self.interface_fluid_nodes:
##                p = node.GetSolutionStepValue(PRESSURE)
##                tmp = p - node.GetSolutionStepValue(PRESSURE,1)
##                p_norm += p*p
##                dp_norm += tmp*tmp;
##            p_norm = math.sqrt(p_norm);
##            dp_norm = math.sqrt(dp_norm);
##            ratio = dp/dp_norm
##            print "dp_norm = ",dp_norm
##            print "p_norm =", p_norm
##            print "fsi convergence ratio ==================== " , ratio
##            if(ratio < 0.005):
##                fsi_is_converged = True
##            else:
##                fsi_is_converged = False
            


	    print "**********************************************","coupled it = ",iteration,"dp =",dp
            iteration = iteration + 1


            #######################################################################
            ########## TRANSFERRING STRUCTURAL SOLUTION TO THE FLUID ##############
            ## map displacements to the structure
            (self.mapper).StructureToFluid_VectorMap(DISPLACEMENT,DISPLACEMENT)

            ##move the mesh -- only on the interface NOT inside the domain
            (self.mesh_solver).solver.MoveNodes();
            #(self.mesh_solver).Solve()

            ## determine the corresponding velocity
            (self.utilities).CopyVectorVar(MESH_VELOCITY,VELOCITY,self.interface_fluid_nodes);
            
            #######################################################################
            ## finalizing fluid solution and getting ready for the next step ######
            ##update the projections
            (self.fluid_solver).ActOnLonelyNodes();
            (self.fluid_solver).SolveStep3();

            #correct the velocity field by applying the pressure
            (self.fluid_solver).SolveStep4();       
                
            if(fsi_is_converged == False ): 
                ## prepare the pressure for the next iteration (overwrites the pressure at the old step)
                (self.utilities).CopyVectorVar(VELOCITY,FRACT_VEL,(self.fluid_model_part).Nodes);
                print "setting the fractional step to the value of velocity for next iteration"
                (self.fluid_solver).SavePressureIteration()

                #saving the acceleration at this iteration -> makes a copy of the acceleration in the
                #from the VariableListDataValueContainer to the DataValueContainer
                (self.utilities).SaveVectorVar(ACCELERATION,ACCELERATION,(self.structure_model_part).Nodes )

        self.step = self.step  + 1
                
        







        
    
