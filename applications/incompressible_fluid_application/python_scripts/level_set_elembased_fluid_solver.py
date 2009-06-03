#importing Kratos main library
from Kratos import *
from KratosConvectionDiffusionApplication import *
from KratosIncompressibleFluidApplication import *

import monolithic_solver_eulerian

def AddVariables(model_part):
    ##displacements 
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    ##velocities
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(FRACT_VEL);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(CONV_PROJ);
##    model_part.AddNodalSolutionStepVariable(CONVECTION_VELOCITY);
    model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
    ##pressures 
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT);
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    ##gravity
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    ##material properties
    model_part.AddNodalSolutionStepVariable(NODAL_MASS);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);

    model_part.AddNodalSolutionStepVariable(NORMAL);
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
##    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
##    model_part.AddNodalSolutionStepVariable(IS_FLUID)
##    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    model_part.AddNodalSolutionStepVariable(ARRHENIUS)
    model_part.AddNodalSolutionStepVariable(IS_DIVIDED)
    ##...aqui lista variables para utilizar
    monolithic_solver_eulerian.AddVariables(model_part)
#adding of Variables to Model Part should be here when the "very fix container will be ready"



    
    print "variables for the level set solver added correctly"

def AddDofs(model_part):
  
    for node in model_part.Nodes:

##        #adding dofs for fluid solution
        node.AddDof(PRESSURE);
        node.AddDof(FRACT_VEL_X);
        node.AddDof(FRACT_VEL_Y);
        node.AddDof(FRACT_VEL_Z);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);

        #adding dofs for convecting the distance function
        node.AddDof(DISTANCE);
##        node.AddDof(TEMPERATURE);
    monolithic_solver_eulerian.AddDofs(model_part)
    print "dofs for the levelset solver added correctly"

class ElemBasedLevelSetSolver:
    
    def __init__(self,model_part,domain_size,body_force):

            #neighbour search
            number_of_avg_elems = 10
            number_of_avg_nodes = 10
            self.mesh_neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
 
            self.model_part = model_part
            self.domain_size = domain_size
            self.body_force = body_force

            ##extrapolation tools
            self.extrapolation_tools = ElemBasedExtrapolationUtilities(model_part)

            ##distance function calculator
            self.distance_calculator = BodyDistanceCalculationUtils()
                
            ##distance tools
            self.distance_tools = ElemBasedDistanceUtilities(model_part)

            ##BC tools
            self.bc_tools = ElemBasedBCUtilities(model_part)
              
##            #assignation of parameters to be used in the strategy (go to line 161)
            self.CalculateReactions = False;
            self.CalculateNormDxFlag = True;
            self.vel_toll = 0.001;
            self.press_toll = 0.001; 
            self.max_vel_its = 3;
            self.max_press_its = 10;
            self.time_order = 1;
            self.laplacian_form = 3; #1 = laplacian, #2 = Discrete Laplacian, #3 discrete laplacian tau=Dt
            self.predictor_corrector = True;
            self.echo_level = 0

            #convection solver order
            self.convection_order = 2 #order of the time scheme of the convection solver
            self.reform_convection_matrix = True
            self.ReformDofAtEachIteration = True
##            self.predict_levelset = True
            self.correct_levelset = True

            #definition of the SOLVERS and related tools
            pConvPrecond = DiagonalPreconditioner()
            self.convection_linear_solver =  BICGSTABSolver(1e-9, 5000,pConvPrecond)

            pDiagPrecond = DiagonalPreconditioner()
##    ##        pILUPrecond = ILU0Preconditioner()
            self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
            self.pressure_linear_solver =  BICGSTABSolver(1e-4, 5000,pDiagPrecond)
##    ##        self.pressure_linear_solver =  BICGSTABSolver(1e-4, 5000,pILUPrecond)

            ##pure convection tool
            if(self.domain_size == 2):
##                self.convection_solver = PureConvectionUtilities2D();
                self.convection_solver = PureConvectionCrankNUtilities2D();
            else: 
##                self.convection_solver = PureConvectionUtilities3D();
                self.convection_solver = PureConvectionCrankNUtilities3D();

            ##redistancing settings
            self.solve_step = 0
            self.dist_recalculation_step = 2
            self.redistance_frequency  = 1
            self.reorder = True

            ##velocity extrapolation distance -- needed to accurately convect the distance function
            self.extrapolation_distance = 1
            self.number_of_extrapolation_layers = 3

    ################################################################
    ################################################################
    def Initialize(self):
        print "entered in initialization"
##        #calculate the normals to the overall domain
##        self.normal_tools.CalculateBodyNormals(self.model_part.Elements,self.domain_size);
        
        #look for neighbours on the base mesh
        (self.mesh_neighbour_search).Execute()

        #constructing the fluid solver
        self.solver = monolithic_solver_eulerian.MonolithicSolver(self.model_part, self.domain_size)
        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU,1)
        self.solver.Initialize()

##        self.solver = ResidualBasedFluidStrategy(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector) 
##        (self.solver).SetEchoLevel(self.echo_level)



        #costruct matrices for convection solver -
        #note that it should be constructed here ONLY
        #if it is fixed during the iterations!!!!!!
        if( self.reform_convection_matrix == False):
            self.convection_solver.ConstructSystem(self.model_part,DISTANCE,VELOCITY,MESH_VELOCITY);

        # Initialize distance function
        self.RecalculateDistanceFunction()

        print "finished initialization"


    
        
    ################################################################
    ################################################################
    def CalculateDistances(self):
        if(self.domain_size == 2):
            self.distance_calculator.CalculateDistances2D(self.model_part.Elements,DISTANCE, self.reorder);
        else:
            self.distance_calculator.CalculateDistances3D(self.model_part.Elements,DISTANCE, self.reorder);
        

    ################################################################
    ################################################################
    #take care! needs neighbours on the overall domain
    def RecalculateDistanceFunction(self):
        print "entered in RecalculateDistanceFunction"
        #mark all nodes outside of the fluid domain and the first layer of nodes inside
        self.distance_tools.MarkExternalAndMixedNodes()

        #change sign
        self.distance_tools.ChangeSignToDistance()

        #calculate distances towards the interior of the domain
        self.CalculateDistances();
        
        #change sign 
        self.distance_tools.ChangeSignToDistance()

        #mark as visited all of the nodes inside the fluid domain
        self.distance_tools.MarkInternalAndMixedNodes()

        #calculate distances towards the outside
        if(self.domain_size == 2):
            self.distance_calculator.CalculateDistances2D(self.model_part.Elements,DISTANCE, True);
        else:
            self.distance_calculator.CalculateDistances3D(self.model_part.Elements,DISTANCE, True);

        #save as distance of the old time step
        self.distance_tools.SaveScalarVariableToOldStep(DISTANCE)
        print "finished RecalculateDistanceFunction"

    ################################################################
    ################################################################
    #take care! needs neighbours on the overall domain
    def Convect(self):
                
        ########### convect distance function ###################
        if( self.ReformDofAtEachIteration == True):
            print "Convect changing the matrix"
            #find neighbours
            (self.mesh_neighbour_search).Execute()            
            
            # construct system -- could be done once if the mesh does not change
            self.convection_solver.ConstructSystem(self.model_part,DISTANCE,VELOCITY,MESH_VELOCITY);

            #calculate projections
            self.convection_solver.CalculateProjection(self.model_part,DISTANCE,NODAL_AREA,VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);
            
            #perform convection step
            self.convection_solver.ConvectScalarVar(self.model_part,self.convection_linear_solver,DISTANCE,VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ,self.convection_order);

            #free memory
            self.convection_solver.ClearSystem()        
            
        else:
            error
            print "Convect without changing the matrix"
            #find neighbours
            (self.mesh_neighbour_search).Execute()
            
            #calculate projections
            self.convection_solver.CalculateProjection(self.model_part,DISTANCE,NODAL_AREA,VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);
            
            #perform convection step
            self.convection_solver.ConvectScalarVar(self.model_part,self.convection_linear_solver,DISTANCE,VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ,self.convection_order);


    
    ################################################################
    ################################################################
    #take care! needs neighbours on the overall domain
    def Extrapolate(self):

        self.extrapolation_tools.ExtrapolateVelocities(self.number_of_extrapolation_layers)
        
    ################################################################
    ################################################################
    def SetModelPart(self):

        self.bc_tools.SetDividedElem_2D()

    def FreeModelPart(self):
        
        self.bc_tools.SetToZeroPressureAndVelocity(self.extrapolation_distance)


    ################################################################
    ################################################################

    ################################################################
    ################################################################
    def Solve(self):
        ############## identifying neighbours #################
        (self.mesh_neighbour_search).Execute()
        #   print "neighbour search finished"


        ## PREDICTION ##
        
        ############## extrapolating by layer v ###############
        self.SetModelPart() #to identify the extrapolation domain
        self.Extrapolate()


        #   print "extrapolation finished"
        ############## convect distance function ##############
        self.Convect()
        #   print "convection finished"
        ############## calculate distances   ##################
        for node in self.model_part.Nodes:
            node.Free(DISTANCE);
            
        if( self.solve_step > self.dist_recalculation_step):
            self.RecalculateDistanceFunction();
            self.dist_recalculation_step += self.redistance_frequency
        #   print "distance calculation finished"
####        self.ComputeSmoothedDensities(delta)

        self.SetModelPart() #to identify the fluid domain to be solved
        print "*******************  setting model part finished         *****************"

        #solve fluid domain
        (self.solver).Solve()
        print "solving procedure finished"
##        for node in self.model_part.Nodes:
##            print node.Id, "    ", node.GetSolutionStepValue(VELOCITY)
        self.FreeModelPart()
        print "freeing pressure and velocity finished"

        #updating the step
        self.solve_step = self.solve_step + 1;

        ## CORRECTION ##
        
####        ############## identifying neighbours ################
####        (self.mesh_neighbour_search).Execute()
##        ############## extrapolating by layer v ##############
        self.Extrapolate()
##        ############## convect distance function #############
##        if(self.correct_levelset == True):
##            self.Convect()
##            print "corrected level set function"
##        ############## calculate distances   ##################
##        if( self.solve_step > self.dist_recalculation_step):
##            self.RecalculateDistanceFunction();
##            self.dist_recalculation_step += self.redistance_frequency
