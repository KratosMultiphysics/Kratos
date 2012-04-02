from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
CheckForPreviousImport()

def AddVariables(base_model_part):
    base_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    base_model_part.AddNodalSolutionStepVariable(VELOCITY);
    base_model_part.AddNodalSolutionStepVariable(FRACT_VEL);
    base_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    base_model_part.AddNodalSolutionStepVariable(PRESSURE);
    base_model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT);
    base_model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    base_model_part.AddNodalSolutionStepVariable(CONV_PROJ);
    base_model_part.AddNodalSolutionStepVariable(NODAL_MASS);
    base_model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    base_model_part.AddNodalSolutionStepVariable(DENSITY);
    base_model_part.AddNodalSolutionStepVariable(VISCOSITY);
    base_model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    base_model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);

    base_model_part.AddNodalSolutionStepVariable(CONVECTION_VELOCITY);
    base_model_part.AddNodalSolutionStepVariable(NORMAL);

    base_model_part.AddNodalSolutionStepVariable(DISTANCE)
    base_model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
    base_model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    base_model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    base_model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    base_model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    base_model_part.AddNodalSolutionStepVariable(IS_FLUID)
    base_model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)

    base_model_part.AddNodalSolutionStepVariable(TEMPERATURE); #to be removed!!!

    
    



    
    print "variables for the level set solver added correctly"

def AddDofs(base_model_part):
  
    for node in base_model_part.Nodes:

        #adding dofs for fluid solution
        node.AddDof(PRESSURE);
        node.AddDof(FRACT_VEL_X);
        node.AddDof(FRACT_VEL_Y);
        node.AddDof(FRACT_VEL_Z);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);

        #adding dofs for convecting the distance function
        node.AddDof(DISTANCE);

    print "dofs for the levelset solver added correctly"


class LevelSetSolver:
    
    def __init__(self,base_model_part,fluid_model_part,domain_size, extrapolation_distance, fluid1_density, fluid1_viscosity, fluid2_density, fluid2_viscosity, fluid1_bodyforce):


        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.base_mesh_neighbour_search = FindNodalNeighboursProcess(base_model_part,number_of_avg_elems,number_of_avg_nodes)
        self.fluid_mesh_neighbour_search = FindNodalNeighboursProcess(fluid_model_part,number_of_avg_elems,number_of_avg_nodes)

        self.base_model_part = base_model_part
        self.fluid_model_part = fluid_model_part
        self.domain_size = domain_size

        #fluid characteristics
        self.fluid1_density = fluid1_density;
        self.fluid1_viscosity = fluid1_viscosity
        self.fluid2_density = fluid2_density;
        self.fluid2_viscosity = fluid2_viscosity
        self.fluid1_bodyforce = fluid1_bodyforce
        

        #nodal area calculator
        self.nodal_area_calculator = CalculateNodalAreaProcess(fluid_model_part,domain_size);

        #assignation of parameters to be used
        self.vel_toll = 0.001;
        self.press_toll = 0.001;
        self.max_vel_its = 10;
        self.max_press_its = 2;
        self.time_order = 2;
        self.CalculateReactions = False;
        self.ReformDofAtEachIteration = True; 
        self.CalculateNormDxFlag = True;
        self.laplacian_form = 2; #1 = laplacian, 2 = Discrete Laplacian
        self.predictor_corrector = True;

        self.echo_level = 0

        

        #convection solver order
        self.convection_order = 2
        pConvPrecond = DiagonalPreconditioner()
        self.convection_linear_solver = linear_solver =  BICGSTABSolver(1e-6, 5000,pConvPrecond)
        self.reform_convection_matrix = False
        self.predict_levelset = True
        self.correct_levelset = True

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        pILUPrecond = ILU0Preconditioner()
        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
##        self.pressure_linear_solver =  BICGSTABSolver(1e-5, 5000,pILUPrecond)
        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)

        ##level set tools
        self.level_set_tools = LevelSetUtilitiesImplicitExtrapolation()

        ##distance function calculator
        self.distance_tools = BodyDistanceCalculationUtils()

        ##calculate normals
        self.normal_tools = BodyNormalCalculationUtils()


        ##pure convection tool
        if(self.domain_size == 2):
            self.convection_solver = PureConvectionUtilities2D();
        else:
            self.convection_solver = PureConvectionUtilities3D();

        ##redistancing settings
        self.solve_step = 0
        self.dist_recalculation_step = 1
        self.redistance_frequency  = 1
        self.reorder = True

        ##velocity extrapolation distance -- needed to accurately convect the distance function
        self.extrapolation_distance = extrapolation_distance
        self.extrapolation_type = 1  #2 = implicit extrapolation -- #1 = explicit  - ATTENTION 2 DOES NOT WORK SATISFACTORELY

    ################################################################
    ################################################################
    def Initialize(self):

        #calculate the normals to the overall domain
        self.normal_tools.CalculateBodyNormals(self.base_model_part,self.domain_size);
        
        #look for neighbours on the base mesh
        (self.base_mesh_neighbour_search).Execute()

        ## identify the fluid domain
        self.level_set_tools.RegenerateFluidModelPart(self.base_model_part,self.fluid_model_part, DISTANCE, self.extrapolation_distance, self.domain_size)

        #constructing the fluid solver
        self.solver = ResidualBasedFluidStrategy(self.fluid_model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)   
        (self.solver).SetEchoLevel(self.echo_level)
        
        print "finished initialization"

    ################################################################
    ################################################################
    ## assumes distance is negative inside and positive outside
    def InitialDistanceCalculation(self):
        (self.base_mesh_neighbour_search).Execute()
        self.level_set_tools.RegenerateFluidModelPart(self.base_model_part,self.fluid_model_part, DISTANCE,self.extrapolation_distance, self.domain_size)

        for node in self.fluid_model_part.Nodes:
            if(node.GetSolutionStepValue(IS_BOUNDARY) == 1.0):
                node.SetValue(IS_VISITED,1.0)
                node.SetSolutionStepValue(DISTANCE,0,0.0)

        print "entered in RecalculateDistanceFunction"
        #prepare to calculate distance inside the fluid
        self.level_set_tools.FluidDistanceComputation_FromBoundary(self.fluid_model_part.Nodes,IS_VISITED,DISTANCE);

        #calculate the distance in the fluid only
        self.CalculateDistances()
        
        #change sign inside the fluid
        self.level_set_tools.SetDistanceToNegative(self.fluid_model_part.Nodes, DISTANCE);

        #mark all the fluid domain as visited and free the rest
        self.level_set_tools.MarkNodesAsVisited(self.base_model_part.Nodes, IS_VISITED);

        #calculate the distance inside the rest of the domain
        self.CalculateDistances()
       
    ################################################################
    ################################################################
    def CalculateDistances(self):
        if(self.domain_size == 2):
            self.distance_tools.CalculateDistances2D(self.base_model_part.Elements,DISTANCE, self.reorder);
        else:
            self.distance_tools.CalculateDistances3D(self.base_model_part.Elements,DISTANCE, self.reorder);
        

    ################################################################
    ################################################################
    #take care! needs neighbours on the overall domain
    def RecalculateDistanceFunction(self):
        print "entered in RecalculateDistanceFunction"
        #prepare to calculate distance inside the fluid
        self.level_set_tools.FluidDistanceComputation_FromBoundary(self.fluid_model_part.Nodes,IS_VISITED,DISTANCE);

        #calculate the distance in the fluid only
        self.CalculateDistances()

        #change sign inside the fluid
        self.level_set_tools.SetDistanceToNegative(self.fluid_model_part.Nodes, DISTANCE);

        #mark all the fluid domain as visited and free the rest
        self.level_set_tools.MarkNodesAsVisited(self.base_model_part.Nodes, IS_VISITED);

        #calculate the distance inside the rest of the domain
        self.CalculateDistances()
        print "finished  RecalculateDistanceFunction"

    ################################################################
    ################################################################
    #take care! needs neighbours on the overall domain
    def Convect(self):
                
        ########### convect distance function ###################
        print "QUIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"        
        # construct system -- could be done once if the mesh does not change
        self.convection_solver.ConstructSystem(self.fluid_model_part,DISTANCE,CONVECTION_VELOCITY,MESH_VELOCITY);

        #calculate projections
        self.convection_solver.CalculateProjection(self.fluid_model_part,DISTANCE,NODAL_AREA,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);
        
        #perform convection step
        self.convection_solver.ConvectScalarVar(self.fluid_model_part,self.convection_linear_solver,DISTANCE,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ,self.convection_order);

        #free memory
        self.convection_solver.ClearSystem();

 
                
        
    ################################################################
    ################################################################
    def Solve(self):

        (self.fluid_mesh_neighbour_search).Execute()
        
        ############## convect distance function #################
        if(self.predict_levelset == True):
            self.Convect()
            print "predicted level set function"
      
        ########### identify & solve fluid ######################
        #identify fluid domain
        (self.base_mesh_neighbour_search).Execute()
        self.level_set_tools.RegenerateFluidModelPart(self.base_model_part,self.fluid_model_part, DISTANCE, self.extrapolation_distance, self.domain_size)

        #apply the fluid properties before solution
        if(self.laplacian_form == 1):
            print "FIXING PRESSURES"
            fix_pressure = True;
        else:
            print "NOT FIXING PRESSURES"
            fix_pressure = False;
        self.level_set_tools.ApplyFluidProperties(self.base_model_part.Nodes, fix_pressure, self.fluid1_density,self.fluid1_viscosity,self.fluid2_density,self.fluid2_viscosity,self.fluid1_bodyforce)


        #detect fluid neighbours
        (self.fluid_mesh_neighbour_search).Execute()

        print self.fluid_model_part

        #solve fluid domain
        (self.solver).ApplyFractionalVelocityFixity()
        (self.solver).Solve()

        print "just before extrapolation:"
        print self.fluid_model_part
        #extrapolate the velocity to the boundaries
        self.level_set_tools.ExtrapolateVelocities(self.base_model_part,self.fluid_model_part,self.extrapolation_distance, DISTANCE, VELOCITY, CONVECTION_VELOCITY)
        
        #updating the step
        self.solve_step = self.solve_step + 1;

        ############## convect distance function #################
        if(self.correct_levelset == True):
            self.Convect()
            print "corrected level set function"
   
        #if required recalculate distance function - 
        if( self.solve_step > self.dist_recalculation_step):
            (self.base_mesh_neighbour_search).Execute()
            self.RecalculateDistanceFunction();
            self.dist_recalculation_step += self.redistance_frequency

        

