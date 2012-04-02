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
    base_model_part.AddNodalSolutionStepVariable(IS_STRUCTURE) ##perche e presente 2 volte?!?!?
    base_model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    base_model_part.AddNodalSolutionStepVariable(IS_FLUID)
    base_model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)


    
    



    
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


class NDLevelSetSolver:
    
    def __init__(self,base_model_part,fluid_model_part,extrapolation_model_part,domain_size, extrapolation_distance,body_force):

        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.base_mesh_neighbour_search = FindNodalNeighboursProcess(base_model_part,number_of_avg_elems,number_of_avg_nodes)
        self.fluid_mesh_neighbour_search = FindNodalNeighboursProcess(fluid_model_part,number_of_avg_elems,number_of_avg_nodes)
        self.extrapolation_mesh_neighbour_search = FindNodalNeighboursProcess(extrapolation_model_part,number_of_avg_elems,number_of_avg_nodes)

        self.base_model_part = base_model_part
        self.fluid_model_part = fluid_model_part
        self.extrapolation_model_part = extrapolation_model_part
        self.domain_size = domain_size
        self.body_force = body_force

        #nodal area calculator
        self.nodal_area_calculator = CalculateNodalAreaProcess(fluid_model_part,domain_size);

        #assignation of parameters to be used
        self.vel_toll = 0.001;
        self.press_toll = 0.001; 
        self.max_vel_its = 3;
        self.max_press_its = 2;
        self.time_order = 1;
        self.CalculateReactions = False;
        self.ReformDofAtEachIteration = True; 
        self.CalculateNormDxFlag = True;
        self.laplacian_form = 2; #1 = laplacian, 2 = Discrete Laplacian
        self.predictor_corrector = True;

        self.echo_level = 0

        

        #convection solver order
        self.convection_order = 1
        pConvPrecond = DiagonalPreconditioner()
        self.convection_linear_solver = linear_solver =  BICGSTABSolver(1e-9, 5000,pConvPrecond)
        self.reform_convection_matrix = False
        self.predict_levelset = True
        self.correct_levelset = True

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
##        pILUPrecond = ILU0Preconditioner()
        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
        self.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pDiagPrecond)
##        self.pressure_linear_solver =  BICGSTABSolver(1e-4, 5000,pILUPrecond)

        ##level set tools
        self.level_set_tools = LevelSetUtilities()

        ##distance function calculator
        self.distance_tools = BodyDistanceCalculationUtils()

        ##calculate normals
        self.normal_tools = BodyNormalCalculationUtils()


        ##pure convection tool
        if(self.domain_size == 2):
            self.convection_solver = PureConvectionCrankNUtilities2D();
        else:
            self.convection_solver = PureConvectionCrankNUtilities3D();

        ##redistancing settings
        self.solve_step = 0
        self.dist_recalculation_step = 1
        self.redistance_frequency  = 1
        self.reorder = True

        ##velocity extrapolation distance -- needed to accurately convect the distance function
        self.extrapolation_distance = extrapolation_distance
        self.extrapolation_type = 2  #2 = implicit extrapolation -- #1 = explicit
        self.extrapolation_max_press_its = 10;

        self.acceptable_distance = 0.0

    ################################################################
    ################################################################
    def Initialize(self):

        #calculate the normals to the overall domain
        self.normal_tools.CalculateBodyNormals(self.base_model_part,self.domain_size);
        
        #look for neighbours on the base mesh
        (self.base_mesh_neighbour_search).Execute()

        ## identify the fluid domain
        self.level_set_tools.RegenerateFluidModelPart(self.base_model_part,self.fluid_model_part, DISTANCE, self.acceptable_distance, self.domain_size)

        #constructing the fluid solver
        self.solver = ResidualBasedNDFluidStrategy(self.fluid_model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)   
        (self.solver).SetEchoLevel(self.echo_level)


        #costruct matrices for convection solver -
        #note that it should be constructed here ONLY
        #if it is fixed during the iterations!!!!!!
        if( self.reform_convection_matrix == False):
            (self.base_mesh_neighbour_search).Execute()
        
            # construct system -- could be done once if the mesh does not change
            self.convection_solver.ConstructSystem(self.base_model_part,DISTANCE,CONVECTION_VELOCITY,MESH_VELOCITY);

        #define the extrapolation solver if needed
        if( self.extrapolation_type == 2):
            vel_toll = 0.001;
            press_toll = 0.001;
            max_vel_its = 5;
            time_order = 1;

            CalculateReactions = False;
            ReformDofAtEachIteration = True; 
            CalculateNormDxFlag = True;
            laplacian_form = 2; #1 = laplacian, 2 = Discrete Laplacian
            predictor_corrector = True;
            self.extrapolation_solver = ResidualBasedNDFluidStrategy(self.extrapolation_model_part,self.velocity_linear_solver,self.pressure_linear_solver, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag, vel_toll, press_toll, max_vel_its, self.extrapolation_max_press_its, time_order, self.domain_size, laplacian_form, predictor_corrector)   
            (self.extrapolation_solver).SetEchoLevel(0)
        else:
            self.number_of_extrapolation_layers = 10;

        
##        print "finished initialization"

    ################################################################
    ################################################################
    ## assumes distance is negative inside and positive outside
    def InitialDistanceCalculation(self):
        (self.base_mesh_neighbour_search).Execute()
        self.level_set_tools.RegenerateFluidModelPart(self.base_model_part,self.fluid_model_part, DISTANCE, self.acceptable_distance, self.domain_size)

        for node in self.fluid_model_part.Nodes:
            if(node.GetSolutionStepValue(IS_BOUNDARY) == 1.0):
                node.SetValue(IS_VISITED,1.0)
                node.SetSolutionStepValue(DISTANCE,0,0.0)

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
        if( self.reform_convection_matrix == True):
            print "QUIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
            #identify domain for convection
            max_dist = 1.0*self.extrapolation_distance
            min_dist = -1.0*self.extrapolation_distance
            self.level_set_tools.GenerateModelPart(self.base_model_part,self.extrapolation_model_part,DISTANCE,IS_VISITED,min_dist,max_dist)
##REMEMBER: GenerateModelPart is defining the extrapolation_model_part in function of extrapolation distance...
##leaving the following Execute, ConstructSystem and CalculateProjection acting on the whole model_part, the above
##is useless.*******************************************************************************************************
            #find neighbours  ##base 
            (self.base_mesh_neighbour_search).Execute()            
            print "find neighbours"

            
            # construct system -- could be done once if the mesh does not change ##base
            self.convection_solver.ConstructSystem(self.base_model_part,DISTANCE,CONVECTION_VELOCITY,MESH_VELOCITY);
            print "construct system"

##            for node in self.extrapolation_model_part.Nodes:
##                print "IDquiiiiiiiiiiiii ", node.Id
##                print "distance ",    node.GetSolutionStepValue(DISTANCE,0)
##                print "CONVEVTION_VELOCITY ",               node.GetSolutionStepValue(CONVECTION_VELOCITY,0)
##                print "MESH_VELOCITY ",                node.GetSolutionStepValue(MESH_VELOCITY,0)
##                print "TEMP_CONV_PROJ ",                node.GetSolutionStepValue(TEMP_CONV_PROJ,0)
##                if (node.GetSolutionStepValue(NODAL_AREA) < 1e-7):
##                    print "NODAL_AREA ",                node.GetSolutionStepValue(NODAL_AREA,0)
##                    print "ID ", node.Id
##                    node.SetSolutionStepValue(NODAL_AREA,0, 1.0);
##                    print "NODAL_AREA ",                node.GetSolutionStepValue(NODAL_AREA,0)


            
##            #calculate projections  ##extrapolate
##            self.convection_solver.CalculateProjection(self.extrapolation_model_part,DISTANCE,NODAL_AREA,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);
####            print "projections"

                #calculate projections   ##base
            self.convection_solver.CalculateProjection(self.base_model_part,DISTANCE,NODAL_AREA,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);
##            print "projections"


##            for element in self.extrapolation_model_part.Elements:
##                  print "ID_elem ", element.Id  
##
##            for node in self.extrapolation_model_part.Nodes:
##               if (node.GetSolutionStepValue(NODAL_AREA) < 1e-7):
##                   print "ID ", node.Id
##                   print "NODAL_AREA ",    node.GetSolutionStepValue(NODAL_AREA,0)

                    
            #perform convection step
            self.convection_solver.ConvectScalarVar(self.extrapolation_model_part,self.convection_linear_solver,DISTANCE,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ,self.convection_order);
##            print "conv step"
##            #Repeat the last 2 steps:
##            #calculate projections
##            self.convection_solver.CalculateProjection(self.extrapolation_model_part,DISTANCE,NODAL_AREA,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);
##            
##            #perform convection step
##            self.convection_solver.ConvectScalarVar(self.extrapolation_model_part,self.convection_linear_solver,DISTANCE,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ,self.convection_order);

                
            #free memory
            self.convection_solver.ClearSystem();

            #recalculate global neighbours
##            (self.base_mesh_neighbour_search).Execute()

            
            
        else:
            #find neighbours   
            (self.base_mesh_neighbour_search).Execute()

##            for node in self.extrapolation_model_part.Nodes:
##                if (node.Id == 1732 ):
##                    print "ID ", node.Id
##                    print "NODAL_AREA ",       node.GetSolutionStepValue(NODAL_AREA,0)
##                    print "distance ",    node.GetSolutionStepValue(DISTANCE,0)
##                if (node.Id == 1741 ):
##                    print "ID ", node.Id
##                    print "NODAL_AREA ",       node.GetSolutionStepValue(NODAL_AREA,0)
##                    print "distance ",    node.GetSolutionStepValue(DISTANCE,0)
##                if (node.Id == 1779 ):
##                    print "ID ", node.Id
##                    print "NODAL_AREA ",       node.GetSolutionStepValue(NODAL_AREA,0)
##                    print "distance ",    node.GetSolutionStepValue(DISTANCE,0)
##                if (node.Id == 1830 ):
##                    print "ID ", node.Id
##                    print "NODAL_AREA ",       node.GetSolutionStepValue(NODAL_AREA,0)
##                    print "distance ",    node.GetSolutionStepValue(DISTANCE,0)
##                if (node.Id == 1848 ):
##                    print "ID ", node.Id
##                    print "NODAL_AREA ",       node.GetSolutionStepValue(NODAL_AREA,0)
##                    print "distance ",    node.GetSolutionStepValue(DISTANCE,0)
##                if (node.Id == 1873 ):
##                    print "ID ", node.Id
##                    print "NODAL_AREA ",       node.GetSolutionStepValue(NODAL_AREA,0)
##                    print "distance ",    node.GetSolutionStepValue(DISTANCE,0)
##                if (node.Id == 1835):
##                    print "ID ", node.Id
##                    print "NODAL_AREA ",       node.GetSolutionStepValue(NODAL_AREA,0)
##                    print "distance ",    node.GetSolutionStepValue(DISTANCE,0)
##                if (node.Id == 1855 ):
##                    print "ID ", node.Id
##                    print "NODAL_AREA ",       node.GetSolutionStepValue(NODAL_AREA,0)
##                    print "distance ",    node.GetSolutionStepValue(DISTANCE,0)
##                if (node.Id == 1893 ):
##                    print "ID ", node.Id
##                    print "NODAL_AREA ",       node.GetSolutionStepValue(NODAL_AREA,0)
##                    print "distance ",    node.GetSolutionStepValue(DISTANCE,0)
##                if (node.Id == 1918 ):
##                    print "ID ", node.Id
##                    print "NODAL_AREA ",       node.GetSolutionStepValue(NODAL_AREA,0)
##                    print "distance ",    node.GetSolutionStepValue(DISTANCE,0)

##                print "CONVEVTION_VELOCITY ",               node.GetSolutionStepValue(CONVECTION_VELOCITY,0)
##                print "MESH_VELOCITY ",                node.GetSolutionStepValue(MESH_VELOCITY,0)
##                print "TEMP_CONV_PROJ ",                node.GetSolutionStepValue(TEMP_CONV_PROJ,0)
##                if (node.GetSolutionStepValue(NODAL_AREA) < 1e-7):
##                    print "NODAL_AREA ",                node.GetSolutionStepValue(NODAL_AREA,0)
##                    print "ID ", node.Id
##                    node.SetSolutionStepValue(NODAL_AREA,0, 1.0);
##                    print "NODAL_AREA ",                node.GetSolutionStepValue(NODAL_AREA,0)


            
            #calculate projections
            self.convection_solver.CalculateProjection(self.base_model_part,DISTANCE,NODAL_AREA,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);
            
            #perform convection step
            self.convection_solver.ConvectScalarVar(self.base_model_part,self.convection_linear_solver,DISTANCE,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ,self.convection_order);

##            #Repeat the last 2 steps:
##            #calculate projections
##            self.convection_solver.CalculateProjection(self.extrapolation_model_part,DISTANCE,NODAL_AREA,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);
##            
##            #perform convection step
##            self.convection_solver.ConvectScalarVar(self.extrapolation_model_part,self.convection_linear_solver,DISTANCE,CONVECTION_VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ,self.convection_order);



    
    ################################################################
    ################################################################
    #take care! needs neighbours on the overall domain
    def Extrapolate(self):

        if(self.extrapolation_type == 1): #explicit extrapolation
            #calculating neighbours on the overall model part
            (self.base_mesh_neighbour_search).Execute()
        
            self.level_set_tools.ExtrapolateVelocities(self.base_model_part,DISTANCE,CONVECTION_VELOCITY,VELOCITY,self.extrapolation_distance);
            
##            self.level_set_tools.ExtrapolateVelocitiesByLayer(self.base_model_part,DISTANCE,CONVECTION_VELOCITY,VELOCITY,self.number_of_extrapolation_layers);
            
        elif(self.extrapolation_type == 2):
            print "BEGINNING EXTRAPOLATION"
            
            #generate a new model part including the extrapolation domain
            self.level_set_tools.ImplicitExtrapolation_PreProcess(self.base_model_part,self.extrapolation_model_part,DISTANCE,VELOCITY,CONVECTION_VELOCITY,self.extrapolation_distance)


##            print "extrapolation model part", self.extrapolation_model_part

            #look for the neighbours in this model part
            (self.extrapolation_mesh_neighbour_search).Execute()

            print "11111111111111111"
            self.level_set_tools.ApplyMinimumExtrapolationPressureFix(self.extrapolation_model_part)
          
            #solve the fluid domain
            (self.extrapolation_solver).ApplyFractionalVelocityFixity()
            (self.extrapolation_solver).Solve()

            #finalize calculations
            self.level_set_tools.ImplicitExtrapolation_PostProcess(self.base_model_part,DISTANCE,CONVECTION_VELOCITY,VELOCITY,self.body_force)

            #calculating neighbours on the overall model part
            (self.base_mesh_neighbour_search).Execute()

            print "FINISHING EXTRAPOLATION"



                
        
    ################################################################
    ################################################################
    def Solve(self):

        (self.base_mesh_neighbour_search).Execute()
        
        ############## PREDICTION: convect distance function #################
        
        if(self.predict_levelset == True):
            self.Convect()
##            print "predicted level set function"
      
        ############# identify & solve fluid #################################
        
        #identify fluid domain
        self.level_set_tools.RegenerateFluidModelPart(self.base_model_part,self.fluid_model_part, DISTANCE, self.acceptable_distance, self.domain_size)
##        print "identify fluid domain   "
        
        #detect fluid neighbours
        (self.fluid_mesh_neighbour_search).Execute()
##        print "detect fluid neighbours   "
         
        print self.fluid_model_part

##        #calculate pressure projection
##        if(self.solve_step > 3):
##            self.nodal_area_calculator.Execute();
##            (self.solver).SolveStep3();

        #solve FLUID DOMAIN
        (self.solver).ApplyFractionalVelocityFixity()
        (self.solver).Solve()
##        print "detect fluid neighbours   "
        
        #updating the step
        self.solve_step = self.solve_step + 1;
##        print "updating the step   "

        #extrapolate velocities 
        self.Extrapolate()
##        print "extrapolate velocities   "
        
        #calculating neighbours on the overall model part
        (self.base_mesh_neighbour_search).Execute()
##        print "calculating neighbours on the overall model part   "
        
        ############## convect distance function #################
        
        if(self.correct_levelset == True):
            self.Convect()
##            print "corrected level set function"
    
        

        
        #if required recalculate distance function - 
        if( self.solve_step > self.dist_recalculation_step):
            self.RecalculateDistanceFunction();
            self.dist_recalculation_step += self.redistance_frequency
##        print " Recalculation   "
