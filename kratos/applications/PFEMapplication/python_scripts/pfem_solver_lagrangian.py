#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.PFEMApplication import *

def AddVariables(model_part):
    for node in model_part.Nodes:
        #adding variables

#        node.GetSolutionStepValue(IS_LAGRANGIAN_INLET)
        node.GetSolutionStepValue(IS_STRUCTURE)
        node.GetSolutionStepValue(IS_BOUNDARY)
        node.GetSolutionStepValue(IS_FREE_SURFACE)
        node.GetSolutionStepValue(NODAL_H)

        node.GetSolutionStepValue(IS_FLUID)
        
        node.GetSolutionStepValue(NORMAL)
        node.GetSolutionStepValue(NORMAL_TO_WALL)
        
        node.GetSolutionStepValue(ACCELERATION)
        node.GetSolutionStepValue(VELOCITY)
        node.GetSolutionStepValue(DISPLACEMENT)
        node.GetSolutionStepValue(FRACT_VEL)
        node.GetSolutionStepValue(MESH_VELOCITY)
        node.GetSolutionStepValue(PRESSURE)
        node.GetSolutionStepValue(PRESSURE,1)
        node.GetSolutionStepValue(PRESS_PROJ)
        node.GetSolutionStepValue(CONV_PROJ)
        node.GetSolutionStepValue(NODAL_AREA)  
                                          
        #adding dofs
        node.AddDof(PRESSURE)
        node.AddDof(FRACT_VEL_X)
        node.AddDof(FRACT_VEL_Y)
        node.AddDof(FRACT_VEL_Z)
        node.AddDof(VELOCITY_X)
        node.AddDof(VELOCITY_Y)
        node.AddDof(VELOCITY_Z)
        

    print "variables for the incompressible fluid solver added correctly"


class PFEMSolver:
    
    def __init__(self,model_part,OuputName,box_corner1,box_corner2,domain_size):

        model_part.ProcessInfo.SetValue(DYNAMIC_TAU, 0.001);

        self.domain_size = domain_size
        self.model_part = model_part
        self.output_name = OuputName
        self.pure_lagrangian = True
        self.close_result_file = True
       
        if(domain_size == 2):
            self.Mesher = TriGenModeler()
            self.neigh_finder = FindNodalNeighboursProcess(model_part,9,18)
        elif (domain_size == 3):
            self.Mesher = TetGenModeler()
            self.neigh_finder = FindNodalNeighboursProcess(model_part,20,30)

        self.Hfinder  = FindNodalHProcess(model_part);  
        self.ActOnWalls  = ActOnWallsNodalProcess(model_part)

        number_of_smoothing_loops = 20
        self.CoordinateSmoother = CoordinateLaplacianSmootherProcess(model_part,number_of_smoothing_loops);
            
        self.MeshMover= MoveMeshProcess(model_part);
        self.EraseNodes = NodeEraseProcess(model_part);

        self.NormalTools = NormalCalculationUtils()
        self.NormalToWall = NormalToWallCalculationUtils()
        self.PfemUtils = PfemUtils()

        self.VolumeCorrectionTools =  VolumeCorrectionUtils()
        self.oldvol = 0.0

        self.alpha_shape = 1.5

        ##saving the limits of the box (all the nodes external to this will be erased)
        self.box_corner1 = box_corner1
        self.box_corner2 = box_corner2
        
        #assignation of parameters to be used
        self.vel_toll = 0.0001;
        self.press_toll = 0.001;
        self.max_vel_its = 7;
        self.max_press_its = 3;  
        self.time_order = 1;
        self.CalculateReactions = False;
        self.ReformDofAtEachIteration = True; 
        self.CalculateNormDxFlag = True;
        self.laplacian_form = 2; #1 = laplacian, 2 = Discrete Laplacian
        self.predictor_corrector = False; #do not set it to true!!!!
        self.smooth = True

        self.step = 0

        

 
        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)

        pILUPrecond = ILU0Preconditioner() 
        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pILUPrecond)
        #self.pressure_linear_solver =  TFQMRSolver(1e-9, 5000,pDiagPrecond)
        #pressure_linear_solver = SkylineLUFactorizationSolver()

        


    def EstimateDeltaTime(self,min_dt,max_dt):
        return (self.PfemUtils).EstimateDeltaTime(min_dt,max_dt,self.model_part)



    def Initialize(self,DeltaTime,output_time_increment):

        print "Initializing PFEM solver"
        
        #finding nodal connectivity
        (self.neigh_finder).Execute();

        domain_size = self.domain_size
        
        #time increment for output 
        self.output_time_increment = output_time_increment
        self.next_output_time = self.output_time_increment
        
        #prediction of the structural motion
        conditions = (self.model_part).Conditions
        (self.NormalTools).CalculateOnSimplex(conditions,domain_size);
        (self.NormalToWall).CalculateNormalToWall(conditions,domain_size);

        #calculating the Hmap after reading the gid mesh
        (self.Hfinder).Execute()

        ##erase the nodes out of the bounding box
        print "bounding box limits"
        print self.box_corner1
        print self.box_corner2
        (self.PfemUtils).IdentifyFluidNodes(self.model_part);
        (self.PfemUtils).MarkOuterNodes(self.box_corner1,self.box_corner2,(self.model_part).Nodes );
        (self.EraseNodes).Execute();

        ##remesh
        (self.Mesher).ReGenerateMesh(self.model_part,self.alpha_shape)

        print "*************************************************************"
        
        ##recalculating neighbours 
        (self.neigh_finder).Execute();

        ##apply the boundary conditions as needed
        (self.PfemUtils).ApplyBoundaryConditions(self.model_part,self.laplacian_form)      

        #initializing the tool to allow the lagrangian inlet
 #       self.LagrangianInlet = LagrangianInletProcess(self.model_part,DeltaTime);


        #initializing the fluid solver
#        self.fluid_solver = ResidualBasedFluidStrategyCoupled(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll, self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)   
        self.fluid_solver = ResidualBasedFluidStrategy(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll, self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)   
        (self.fluid_solver).SetEchoLevel(0)



        (self.lagrangian_tools) = LagrangianUtils()
        (self.VariableUtils) = VariableUtils()

        print "Intialize PFEM solver terminated"
       
        
           
    ######################################################################
    def Solve(self,time,gid_io):

        self.PredictionStep(time)
        self.FluidSolutionStep()

        self.PostSolutionStep()

        self.OutputStep(time,gid_io)

    ######################################################################
    def SolutionStep1(self):
        #step initialization for the fluid solution
        (self.fluid_solver).InitializeFractionalStep(self.step, self.time_order);
	(self.fluid_solver).InitializeProjections(self.step);
	(self.fluid_solver).AssignInitialStepValues();

	normDx = Array3(); normDx[0] = 0.00; normDx[1] = 0.00; normDx[2] = 0.00;
        is_converged = False
        iteration = 0

        while(	is_converged == False and iteration < self.max_vel_its  ):
            (self.VariableUtils).SaveVectorVar(FRACT_VEL,VAUX,(self.model_part).Nodes )
	    (self.fluid_solver).FractionalVelocityIteration(normDx);
            (self.lagrangian_tools).CalculateStep1DisplacementCorrection((self.model_part))
            (self.MeshMover).Execute();
            is_converged = (self.fluid_solver).ConvergenceCheck(normDx,self.vel_toll);
            iteration = iteration + 1
            
    ######################################################################
    #in this step the node position is predicted, the mesh is regenerated
    #and the boundary conditions are imposed
    def PredictionStep(self,time):
        domain_size = self.domain_size
               
        #calculate normals
        conditions = (self.model_part).Conditions
        (self.NormalTools).CalculateOnSimplex(conditions,domain_size)
        (self.NormalToWall).CalculateNormalToWall(conditions,domain_size)



        #lagrangian_prediction
        (self.lagrangian_tools).ExplicitLagrangianPrediction((self.model_part))

        (self.PfemUtils).MoveLonelyNodes(self.model_part)

        (self.MeshMover).Execute();

        ##ensure that no node gets too close to the walls
        (self.ActOnWalls).Execute();

        ##move the mesh
        (self.MeshMover).Execute();

        if( self.pure_lagrangian == False):

            ##ensure that no node gets too close to the walls
            (self.ActOnWalls).Execute();

            ##move the mesh
            (self.MeshMover).Execute();

            ##smooth the final position of the nodes to homogeneize the mesh distribution
            if(self.smooth == True):
                (self.CoordinateSmoother).Execute();

            ##move the mesh
            (self.MeshMover).Execute();

        #regenerate the mesh
        self.Remesh()

    ######################################################################
    def Remesh(self):
        ##erase all conditions and elements prior to remeshing
        ((self.model_part).Elements).clear();
        ((self.model_part).Conditions).clear();

        ##erase the nodes out of the bounding box
        (self.PfemUtils).MarkOuterNodes(self.box_corner1,self.box_corner2,(self.model_part).Nodes );
        (self.EraseNodes).Execute();

        ##remesh
        (self.Mesher).ReGenerateMesh(self.model_part,self.alpha_shape)

        print self.model_part

        ##recalculating neighbours 
        (self.neigh_finder).Execute();

        ##apply the boundary conditions as needed
        (self.PfemUtils).ApplyBoundaryConditions(self.model_part,self.laplacian_form);
        
    ######################################################################
    #ALE solution of the fluid on the moving domain
    def FluidSolutionStep(self):

        
        ##solve the fluid problem

        (self.fluid_solver).ApplyFractionalVelocityFixity()

        if(self.step < self.time_order + 1):
            (self.fluid_solver).InitializeFractionalStep(self.step, self.time_order);
            (self.fluid_solver).InitializeProjections(self.step);
            (self.fluid_solver).AssignInitialStepValues();
            (self.fluid_solver.IterativeSolve())
        else:
            if(self.predictor_corrector == False):
                self.PureFractionalStepFluidSolutionStep()
            else:
                self.PredictorCorrectorFluidSolutionStep()
       
        self.step = self.step + 1

        (self.fluid_solver).Clear()

    ######################################################################
    def PureFractionalStepFluidSolutionStep(self):
        self.SolutionStep1()
        
        (self.fluid_solver).SolveStep2();
        (self.fluid_solver).ActOnLonelyNodes();
        (self.fluid_solver).SolveStep3();

        (self.VariableUtils).SaveVectorVar(FRACT_VEL,VAUX,(self.model_part).Nodes )
        (self.fluid_solver).SolveStep4();
        (self.lagrangian_tools).CalculateFinalDisplacementCorrection((self.model_part))
        (self.MeshMover).Execute();
    
      
    ######################################################################
    #ALE solution of the fluid on the moving domain
    def PredictorCorrectorFluidSolutionStep(self):

        it = 1
        pressure_is_converged = False
        while(pressure_is_converged == False and it <= self.max_press_its):

##            inverted_elements = False
##            vol = (self.PfemUtils).CalculateVolume(self.model_part,self.domain_size,inverted_elements)
##            print vol
##            if(inverted_elements == True):
##                print "**********************************************************"
##                print "emergency remesh"
##                print "**********************************************************"
##                self.Remesh()
            
            self.SolutionStep1()
            
            dp = (self.fluid_solver).SolveStep2();
            (self.fluid_solver).ActOnLonelyNodes();
            (self.fluid_solver).SolveStep3();

            (self.VariableUtils).SaveVectorVar(FRACT_VEL,VAUX,(self.model_part).Nodes )
            (self.fluid_solver).SolveStep4();
            (self.lagrangian_tools).CalculateFinalDisplacementCorrection((self.model_part))
            (self.MeshMover).Execute();

            p_norm = (self.fluid_solver).SavePressureIteration()

            it = it + 1
            print "********************* ",it, " ******", dp/p_norm , " **************"
            if(dp/p_norm < self.press_toll):
                pressure_is_converged = True
       
        self.step = self.step + 1

        (self.fluid_solver).Clear()



    ######################################################################
    #nodes in the inlet introduced and ALE acceleration calculated
    def PostSolutionStep(self):
        print "doing nothing in the post solution step"
 
        ##estimate error
        ##error_estimator.Execute();

        ##lagrangian fluid step
        ##(self.LagrangianInlet).Execute();


    ######################################################################
    def OutputStep(self,time,gid_io):
        if(time >= self.next_output_time):
            self.next_output_time = self.next_output_time + self.output_time_increment;
            
            file_name = self.output_name
            file_name = file_name + str(time)

            gid_io.ChangeOutputName(file_name,GiDPostMode.GiD_PostBinary);
            gid_io.WriteMesh((self.model_part).GetMesh(),self.domain_size,GiDPostMode.GiD_PostBinary);

            gid_io.WriteNodalResults(NODAL_AREA, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(NORMAL, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(NORMAL_TO_WALL, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(PRESSURE, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(IS_FREE_SURFACE, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(IS_BOUNDARY, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(IS_STRUCTURE, (self.model_part).Nodes, time, 0);
  #  	    gid_io.WriteNodalResults(IS_LAGRANGIAN_INLET, (self.model_part).Nodes, time, 0);
    	    gid_io.WriteNodalResults(IS_FLUID, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(VELOCITY, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(MESH_VELOCITY, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(ACCELERATION, (self.model_part).Nodes, time, 0);

            if( self.close_result_file == True ):
                gid_io.Flush()
                gid_io.CloseResultFile();

    


