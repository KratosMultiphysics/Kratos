#importing Kratos main library
from Kratos import *
from KratosConvectionDiffusionApplication import *
from KratosIncompressibleFluidApplication import *
import math

import monolithic_solver_eulerian

def AddVariables(model_part):
    ##displacements 
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    ##velocities
    model_part.AddNodalSolutionStepVariable(VELOCITY);
##    model_part.AddNodalSolutionStepVariable(FRACT_VEL);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(CONV_PROJ);
    model_part.AddNodalSolutionStepVariable(CONVECTION_VELOCITY); 
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
    model_part.AddNodalSolutionStepVariable(POROSITY);
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
    model_part.AddNodalSolutionStepVariable(AUX_INDEX)
    model_part.AddNodalSolutionStepVariable(DIAMETER)
##...aqui lista variables para utilizar
    

#adding the variables of the Monolithic Solver
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
            ##redistance calculator
            if(self.domain_size == 2):
                self.distance_utils = SignedDistanceCalculationUtils2D();
            else: 
                self.distance_utils = SignedDistanceCalculationUtils3D();
                
            ##distance tools
            self.distance_tools = ElemBasedDistanceUtilities(model_part)

            ##BC tools
            self.bc_tools = ElemBasedBCUtilities(model_part)

            #convection solver setting
            self.convection_order = 2 #order of the time scheme of the convection solver
            self.reform_convection_matrix = True
            self.ReformDofAtEachIteration = True

            #definition of the CONVECTION-SOLVERS and related tools
            pConvPrecond = DiagonalPreconditioner()
            self.convection_linear_solver =  BICGSTABSolver(1e-9, 5000,pConvPrecond)

            ##pure convection tool
            if(self.domain_size == 2):
##                self.convection_solver = PureConvectionUtilities2D();
                self.convection_solver = PureConvectionCrankNUtilities2D();
            else: 
##                self.convection_solver = PureConvectionUtilities3D();
                self.convection_solver = PureConvectionCrankNUtilities3D();

            ##redistancing settings
            self.solve_step = 0
            self.dist_recalculation_step = 1
            self.redistance_frequency  = 1
            self.reorder = True

            self.number_of_extrapolation_layers = 3
            

    ################################################################
    ################################################################
    def Initialize(self):
        print "entered in initialization"
##        #calculate the normals to the overall domain
##        self.normal_tools.CalculateBodyNormals(self.model_part.Elements,self.domain_size);
##
##        print "to be removed SOOOOOOOOOOOOOOOOOOOOOOOOOONNNNNNNNNN!!!!!"
##        for node in self.model_part.Nodes:
##            node.SetSolutionStepValue(DIAMETER,0,0.01)

        #look for neighbours on the base mesh
        (self.mesh_neighbour_search).Execute()

        #constructing the fluid solver
        self.solver = monolithic_solver_eulerian.MonolithicSolver(self.model_part, self.domain_size)
        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU,0)
        self.max_iter = 10
        self.solver.Initialize()


        #costruct matrices for convection solver -
        #note that it should be constructed here ONLY
        #if it is fixed during the iterations!!!!!!
        if( self.reform_convection_matrix == False):
            self.convection_solver.ConstructSystem(self.model_part,DISTANCE,VELOCITY,MESH_VELOCITY);

        # Initialize distance function
        self.RecalculateDistanceFunction()

        print "finished initialization"


    
        
    ################################################################
##    ################################################################
##    def CalculateDistances(self):
##        if(self.domain_size == 2):
##            self.distance_calculator.CalculateDistances2D(self.model_part.Elements,DISTANCE, self.reorder);
##        else:
##            self.distance_calculator.CalculateDistances3D(self.model_part.Elements,DISTANCE, self.reorder);
        

    ################################################################
    ################################################################
    #take care! needs neighbours on the overall domain
    def RecalculateDistanceFunction(self):
        self.distance_utils.CalculateDistances(self.model_part,DISTANCE, 10000000.0)


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
            self.convection_solver.ConstructSystem(self.model_part,DISTANCE,CONVECTION_VELOCITY,MESH_VELOCITY);

##            #calculate projections version#1 
##            self.convection_solver.CalculateProjection(self.model_part,DISTANCE,NODAL_AREA,VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);
            
            #calculate projections version#2 convecting with darcy_velocity / epsilon
            self.convection_solver.CalculateProjection(self.model_part,DISTANCE,NODAL_AREA, CONVECTION_VELOCITY ,MESH_VELOCITY,TEMP_CONV_PROJ);
            
            #perform convection step
            self.convection_solver.ConvectScalarVar(self.model_part,self.convection_linear_solver,DISTANCE,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ,self.convection_order);

            #free memory
            self.convection_solver.ClearSystem()        
            
        else:
            print "Convect without changing the matrix"
            #find neighbours
            (self.mesh_neighbour_search).Execute()

##            #calculate projections version#1 
##            self.convection_solver.CalculateProjection(self.model_part,DISTANCE,NODAL_AREA,VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);
            
            #calculate projections version#2 convecting with darcy_velocity / epsilon
            self.convection_solver.CalculateProjection(self.model_part,DISTANCE,NODAL_AREA, CONVECTION_VELOCITY ,MESH_VELOCITY,TEMP_CONV_PROJ);
            
            #perform convection step
            self.convection_solver.ConvectScalarVar(self.model_part,self.convection_linear_solver,DISTANCE,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ,self.convection_order);


    
    ################################################################
    ################################################################
    #take care! needs neighbours on the overall domain
    def Extrapolate(self):

        self.extrapolation_tools.ExtrapolateVelocities(self.number_of_extrapolation_layers)
        
    ################################################################
    ################################################################
    def SetModelPart(self):

        self.bc_tools.SetDividedElem_2D()

##    def FreeModelPart(self):
##        
##        self.bc_tools.SetToZeroPressureAndVelocity(self.extrapolation_distance)
##        print "MODEL PART FREE *****************************************************************"
    ################################################################
    ################################################################

    def CalculateDelta_t(self, delta_t_in):


        for node in self.model_part.Nodes:

            velx = node.GetSolutionStepValue(VELOCITY_X,0)
            vely = node.GetSolutionStepValue(VELOCITY_Y,0)
            velz = node.GetSolutionStepValue(VELOCITY_Z,0)
            vel = math.sqrt(velx*velx + vely*vely + velz*velz)
            eps = node.GetSolutionStepValue(POROSITY,0)
            if (eps == 0.0):
                eps = 1.0
                
            density = 1000;
            mu = 1e-3;
            dp = 0.01;
            
            kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*dp*dp);
    ##        print kinv
    ##        print max_vel

            A = kinv * mu 
    ##        print A
            B = 1.75 * density /eps * vel * math.sqrt( kinv / (eps * 150.0))
    ##        print B

            temp_delta_t = 1e6
            delta_t = delta_t_in
            if((A+B) != 0.0):
                if(node.GetSolutionStepValue(DISTANCE) <= 0.0):
                    temp_delta_t = 1/(A + B)


            if (temp_delta_t < delta_t_in):
                delta_t = temp_delta_t
            
##        print delta_t

        return delta_t

        
    ################################################################
    ################################################################

    ################################################################
    ################################################################
    def Solve(self):
        ############## identifying neighbours #################
        (self.mesh_neighbour_search).Execute()
##        print "neighbour search finished"


        ## PREDICTION ##
        
        ############## extrapolating by layer v ###############
        self.SetModelPart() #to identify the extrapolation domain
##        print "1st setting model part finished"
        self.Extrapolate()


##        print "1st extrapolation finished"
        ############## convect distance function ##############
        for node in self.model_part.Nodes:
            eps = node.GetSolutionStepValue(POROSITY,0)
            if (eps == 0.0):
                eps = 1.0
            diam = node.GetSolutionStepValue(DIAMETER,0)
            if (diam == 0.0):
                node.SetSolutionStepValue(DIAMETER,0,1.0)
            
            vx = node.GetSolutionStepValue(VELOCITY_X,0)/eps
            vy = node.GetSolutionStepValue(VELOCITY_Y,0)/eps
            vz = node.GetSolutionStepValue(VELOCITY_Z,0)/eps
            node.SetSolutionStepValue(CONVECTION_VELOCITY_X,0,vx)
            node.SetSolutionStepValue(CONVECTION_VELOCITY_Y,0,vy)
            node.SetSolutionStepValue(CONVECTION_VELOCITY_Z,0,vz)
##            if (node.GetSolutionStepValue(VELOCITY_X) != 0.0):
##                if (node.GetSolutionStepValue(POROSITY) == 0.5):
##                    print node.Id
##                    print node.GetSolutionStepValue(VELOCITY,0)
##                    print node.GetSolutionStepValue(CONVECTION_VELOCITY,0)

                
        self.Convect()
##        print "convection finished"

##        for node in self.model_part.Nodes:
##            if (node.GetSolutionStepValue(DISTANCE) >= 0.0):
##                if(node.X >= 0.1):
##                    if(node.X<= 9.9):
##                        node.SetSolutionStepValue(VELOCITY_X,0,0.0)
##                        node.SetSolutionStepValue(VELOCITY_Y,0,0.0)
##                        node.SetSolutionStepValue(VELOCITY_Z,0,0.0)
        
        
        ############## calculate distances   ##################
##        for node in self.model_part.Nodes:
##            node.Free(DISTANCE);
            
        if( self.solve_step > self.dist_recalculation_step):
            self.RecalculateDistanceFunction();
            self.dist_recalculation_step += self.redistance_frequency
##        print "distance calculation finished"
####        self.ComputeSmoothedDensities(delta)

        self.SetModelPart() #to identify the fluid domain to be solved
##        print "2nd setting model part finished"

        

        #solve fluid domain
        (self.solver).Solve()
##        print "****************      solving procedure finished     ***********************"


        #updating the step
        self.solve_step = self.solve_step + 1;

        ## CORRECTION ##
        
##        ############## identifying neighbours ################
##        (self.mesh_neighbour_search).Execute()
        ############## extrapolating by layer v ##############
        self.Extrapolate()
##        print "2nd extrapolation finished"

        ############## convect distance function #############
        for node in self.model_part.Nodes:
            eps = node.GetSolutionStepValue(POROSITY,0)
            if (eps == 0.0):
                eps = 1.0
            vx = node.GetSolutionStepValue(VELOCITY_X,0)/eps
            vy = node.GetSolutionStepValue(VELOCITY_Y,0)/eps
            vz = node.GetSolutionStepValue(VELOCITY_Z,0)/eps
            node.SetSolutionStepValue(CONVECTION_VELOCITY_X,0,vx)
            node.SetSolutionStepValue(CONVECTION_VELOCITY_Y,0,vy)
            node.SetSolutionStepValue(CONVECTION_VELOCITY_Z,0,vz)
##            if (node.GetSolutionStepValue(VELOCITY_X) != 0.0):
##                if (node.GetSolutionStepValue(POROSITY) == 0.5):
##                    print node.Id
##                    print node.GetSolutionStepValue(VELOCITY,0)
##                    print node.GetSolutionStepValue(CONVECTION_VELOCITY,0)

        
        self.Convect()
        
##            print "corrected level set function"
        ############## calculate distances   ##################
        if( self.solve_step > self.dist_recalculation_step):
            self.RecalculateDistanceFunction();
            self.dist_recalculation_step += self.redistance_frequency
        
