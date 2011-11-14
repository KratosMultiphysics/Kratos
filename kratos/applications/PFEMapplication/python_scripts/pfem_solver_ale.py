# -*- coding: utf-8 -*-
#importing the Kratos Library
from Kratos import *
from KratosIncompressibleFluidApplication import *
from KratosPFEMApplication import *
from KratosMeshingApplication import *
import time

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(NORMAL);
    model_part.AddNodalSolutionStepVariable(NORMAL_TO_WALL);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(FRACT_VEL);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT);
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    model_part.AddNodalSolutionStepVariable(CONV_PROJ);
##    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(NODAL_MASS);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(IS_FLUID);
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(ARRHENIUS);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    print "variables for the incompressible fluid solver added correctly"

def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(PRESSURE);
        node.AddDof(FRACT_VEL_X);
        node.AddDof(FRACT_VEL_Y);
        node.AddDof(FRACT_VEL_Z);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);


        
        node.AddDof(IS_STRUCTURE);
    print "dofs for the incompressible fluid solver added correctly"


class PFEMSolver:
    
    def __init__(self,model_part,OuputName,box_corner1,box_corner2,domain_size):
        model_part.ProcessInfo.SetValue(DYNAMIC_TAU, 0.001);

        self.domain_size = domain_size
        self.model_part = model_part
        self.output_name = OuputName
        self.pure_lagrangian = False
        self.correct_volume = True
        self.prediction_order = 2
        self.echo_level = 0

        self.mesh_every_nsteps = 10;
        self.substep = 0;
       
        if(domain_size == 2):
            self.Mesher = TriGenPFEMModeler()
#            self.Mesher = TriGenModeler()
            self.neigh_finder = FindNodalNeighboursProcess(model_part,9,18)
        elif (domain_size == 3):
            self.Mesher = TetGenPfemModeler()
            self.neigh_finder = FindNodalNeighboursProcess(model_part,20,30)

        self.Hfinder  = FindNodalHProcess(model_part);  
        self.ActOnWalls  = ActOnWallsNodalProcess(model_part)

        self.add_nodes = True
        self.remove_nodes = True
        number_of_smoothing_loops = 20
        reduction_factor = 1.0
        self.CoordinateSmoother = CoordinateLaplacianSmootherProcess(model_part,number_of_smoothing_loops,reduction_factor);
            
        self.MeshMover= MoveMeshProcess(model_part);
        self.EraseNodes = NodeEraseProcess(model_part);

        self.NormalTools = NormalCalculationUtils()
        self.NormalToWall = NormalToWallCalculationUtils()
        self.PfemUtils = PfemUtils()

        if(self.correct_volume == True):
            self.VolumeCorrectionTools =  VolumeCorrectionUtils()

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
        self.laplacian_form = 3; #1 = laplacian, 2 = Discrete Laplacian; 3 = discrete laplacian + tau=dt
        self.predictor_corrector = False; 
        self.smooth = True
        self.h_factor = 0.2 #nodes closer than h_factor * H will be erased

        self.step = 0

        self.displacement_correction_afterstep4 = False

        self.eulerian_lagrangian = True

        self.uzawa_term_is_active = False

        #aux variable
        self.normDx = 0.0
	self.total_remeshing_time = 0.00
        

 
        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)

##        pILUPrecond = ILU0Preconditioner() 
##        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pILUPrecond)

        pDiagPrecond = DiagonalPreconditioner()
        self.pressure_linear_solver =  BICGSTABSolver(1e-4, 5000,pDiagPrecond)
        #self.pressure_linear_solver =  TFQMRSolver(1e-9, 5000,pDiagPrecond)
        #pressure_linear_solver = SkylineLUFactorizationSolver()

        self.projections_are_initialized = True;

        


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

        #identify the fluid nodes
        (self.PfemUtils).ApplyBoundaryConditions(self.model_part,self.laplacian_form)      
        (self.PfemUtils).IdentifyFluidNodes(self.model_part);
        (self.PfemUtils).ApplyMinimalPressureConditions(self.model_part);

        ##erase the nodes out of the bounding box
        print "bounding box limits"
        print self.box_corner1
        print self.box_corner2
        
        (self.PfemUtils).MarkOuterNodes(self.box_corner1,self.box_corner2,(self.model_part).Nodes );
        (self.EraseNodes).Execute();

        ##remesh
        if(self.domain_size==2):
            (self.Mesher).ReGenerateMesh("Fluid2D","Condition2D",self.model_part,self.EraseNodes, self.remove_nodes, self.add_nodes, self.alpha_shape, self.h_factor)
        else:
            (self.Mesher).ReGenerateMesh("Fluid3D","Condition3D",self.model_part,self.EraseNodes, self.remove_nodes, self.add_nodes, self.alpha_shape, self.h_factor)
            
#        (self.Mesher).ReGenerateMesh(self.model_part,self.alpha_shape)
        print "remeshing in initalize performed succesfully"

        print "*************************************************************"
        
        ##recalculating neighbours 
        (self.neigh_finder).Execute();

        ##apply the boundary conditions as needed
        (self.PfemUtils).ApplyBoundaryConditions(self.model_part,self.laplacian_form)      
        (self.PfemUtils).IdentifyFluidNodes(self.model_part);
        (self.PfemUtils).ApplyMinimalPressureConditions(self.model_part);

        #initializing the tool to allow the lagrangian inlet
 #       self.LagrangianInlet = LagrangianInletProcess(self.model_part,DeltaTime);

        print "time order = ", self.time_order
        #initializing the fluid solver
        print self.vel_toll
        self.domain_size = int(self.domain_size)
        self.laplacian_form = int(self.laplacian_form)
        solver_configuration = FractionalStepConfiguration(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.domain_size,self.laplacian_form )

        self.ReformDofAtEachIteration = bool(self.ReformDofAtEachIteration)
        self.vel_toll = float(self.vel_toll)
        self.press_toll = float(self.press_toll)
        self.max_vel_its = int(self.max_vel_its)
        self.max_press_its = int(self.max_press_its)
        self.time_order = int(self.time_order)
        self.domain_size = int(self.domain_size)
        self.predictor_corrector = bool(self.predictor_corrector)
        self.fluid_solver = FractionalStepStrategy( self.model_part, solver_configuration, self.ReformDofAtEachIteration, self.vel_toll, self.press_toll, self.max_vel_its, self.max_press_its, self.time_order, self.domain_size,self.predictor_corrector)

#        self.fluid_solver = ResidualBasedFluidStrategyCoupled(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll, self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)   
##        self.fluid_solver = ResidualBasedFluidStrategy(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll, self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)
##        self.fluid_solver = FractionalStepStrategy(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll, self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)
##        (self.fluid_solver).SetEchoLevel(self.echo_level)
        #determine the original volume
        if(self.correct_volume == True):
            [inverted_elements,vol] = self.CheckForInvertedElements()
            self.originalvolume = vol
            
        (self.lagrangian_tools) = LagrangianUtils()
        (self.VariableUtils) = VariableUtils()

        self.utilities = VariableUtils()


        print "Intialize PFEM solver terminated"
       
        
           
    ######################################################################
    def Solve(self,time,gid_io):
        self.FindNeighbours()

        if(self.eulerian_lagrangian == True):
            print "SONO QUI"
            self.EulerianLagrangianPredictionStep(time)
        else:
            print "SONO LI"
            self.PredictionStep(time)
            
        self.FluidSolutionStep()

        self.PostSolutionStep()

        self.OutputStep(time,gid_io)

            
    ######################################################################
    def CheckForInvertedElements(self):
        volume = (self.PfemUtils).CalculateVolume(self.model_part,self.domain_size)
        inverted_elements = False
        if(volume < 0.0):
            volume = - volume
            inverted_elements = True
        return [inverted_elements,volume]
        
    ######################################################################
    def ReduceTimeStep(self,model_part):
        (self.PfemUtils).ReduceTimeStep(model_part);
        
    ######################################################################
    def FindNeighbours(self):
        (self.neigh_finder).Execute();
        
    ######################################################################
    #in this step the node position is predicted, the mesh is regenerated
    #and the boundary conditions are imposed
    def PredictionStep(self,time):

        if(self.step > 2):
            print "*************** second order prediction ***************"  
            domain_size = self.domain_size
                  
            #calculate normals
            conditions = (self.model_part).Conditions
            print "pfem_solver_ale line 257"
            (self.NormalTools).CalculateOnSimplex(conditions,domain_size)
            print "pfem_solver_ale line 259"
            (self.NormalToWall).CalculateNormalToWall(conditions,domain_size)
            print "pfem_solver_ale line 261"
            #performing a first order prediction of the fluid displacement
            (self.PfemUtils).Predict(self.model_part)
            (self.MeshMover).Execute();

            #move free nodes and treat nodes close to the boundary
            (self.PfemUtils).MoveLonelyNodes(self.model_part)
            (self.ActOnWalls).Execute();        
            (self.MeshMover).Execute();

            #remesh
            self.Remesh()
            print self.model_part

            #improve the prediction to reach second order
            if(self.prediction_order == 2 and self.step > 2):
                #lagrangian improvement of the correction
                self.LagrangianCorrection()

            #smooth the mesh position
            if(self.smooth == True):
                (self.CoordinateSmoother).Execute();
                (self.MeshMover).Execute();
        
            #verify if the mesh is suitable to go on calculations and eventually remesh
            [inverted_elements,vol] = self.CheckForInvertedElements()
            print "vol = " , vol
            if(inverted_elements == True):
                self.Remesh()
                print "********************************************* "
                print "**** PERFORMED EMERGENCY REMESH ************* "
                print "********************************************* "

            self.CalculateProjections()

            
   
            

        
    ######################################################################
    def Remesh(self):
        remeshing_start_time=time.time()
        
        #clearing the sytem matrices as the connectivity will change
        (self.fluid_solver).Clear()
        print "fluid solver"

        ##erase the nodes out of the bounding box
        print "h_factor=", self.h_factor        
        (self.PfemUtils).MarkOuterNodes(self.box_corner1,self.box_corner2,(self.model_part).Nodes );
        (self.PfemUtils).MarkExcessivelyCloseNodes((self.model_part).Nodes,self.h_factor);
        (self.EraseNodes).Execute();
        print "erase is executed"

        (self.neigh_finder).ClearNeighbours();
        print "neighbours are cleared"
        
        ##erase all conditions and elements prior to remeshing
        ((self.model_part).Elements).clear();
        ((self.model_part).Conditions).clear();
        print "elements are cleared"
        
        ##remesh
        #(self.Mesher).ReGenerateMesh(self.model_part,self.alpha_shape)
        if(self.domain_size==2):
            (self.Mesher).ReGenerateMesh("Fluid2D","Condition2D",self.model_part,self.EraseNodes, self.remove_nodes, self.add_nodes, self.alpha_shape,self.h_factor)
        else:
            (self.Mesher).ReGenerateMesh("Fluid3D","Condition3D",self.model_part,self.EraseNodes, self.remove_nodes, self.add_nodes, self.alpha_shape, self.h_factor)

        print self.model_part

        print "a"

        ##recalculating neighbours 
        (self.neigh_finder).Execute();

        ##apply the boundary conditions as needed
        print "b"
        (self.PfemUtils).ApplyBoundaryConditions(self.model_part,self.laplacian_form);
        print "c"
        (self.PfemUtils).IdentifyFluidNodes(self.model_part);
        print "d"
        (self.PfemUtils).ApplyMinimalPressureConditions(self.model_part);
        print "e"

        self.substep = 0;
        
        remeshing_end_time=time.time()
        
        print "Remeshing time :", remeshing_end_time - remeshing_start_time
        self.total_remeshing_time += remeshing_end_time - remeshing_start_time
        print "Total remeshing time :", self.total_remeshing_time
        
        
    ######################################################################
    def SolutionStep1(self):
       #step initialization for the fluid solution
        (self.fluid_solver).InitializeFractionalStep(self.step, self.time_order);
        (self.fluid_solver).InitializeProjections(self.step,self.projections_are_initialized);
        self.projections_are_initialized = True;
        (self.fluid_solver).AssignInitialStepValues();

        self.normDx = 0.0 #[0] = 0.00; self.normDx[1] = 0.00; self.normDx[2] = 0.00;
        is_converged = False
        iteration = 0

        #iterative solution of the velocity
        while(  is_converged == False and iteration < self.max_vel_its  ):
            self.normDx = (self.fluid_solver).FractionalVelocityIteration();
            is_converged = (self.fluid_solver).ConvergenceCheck(self.normDx,self.vel_toll);
            iteration = iteration + 1
            
    ######################################################################
    def LagrangianCorrection(self):
        print "Lagrangian prediction"

        #step initialization for the fluid solution
        (self.fluid_solver).ApplyFractionalVelocityFixity()
        (self.fluid_solver).InitializeFractionalStep(self.step, self.time_order);
        (self.fluid_solver).InitializeProjections(self.step,self.projections_are_initialized);
        self.projections_are_initialized = True;
        (self.fluid_solver).AssignInitialStepValues();

        print "before lagrangian solution"

        (self.VariableUtils).SaveVectorVar(FRACT_VEL,VAUX,(self.model_part).Nodes )
        self.normDx = (self.fluid_solver).FractionalVelocityIteration();

        is_converged = (self.fluid_solver).ConvergenceCheck(self.normDx,self.vel_toll);

        (self.lagrangian_tools).CalculateStep1DisplacementCorrection((self.model_part))
        print "after lagrangian solution"

        #deleting unnecessary memory
##      (self.fluid_solver).Clear()

        #return normDx
            
    ######################################################################
    #ALE solution of the fluid on the moving domain
    def FluidSolutionStep(self):
        (self.fluid_solver).ApplyFractionalVelocityFixity()

        if(self.step < self.time_order + 2):
            print "IterativeSolve"
#            self.PredictorCorrectorFluidSolutionStep()
            (self.fluid_solver).InitializeFractionalStep(self.step, self.time_order);
            (self.fluid_solver).InitializeProjections(self.step,self.projections_are_initialized);
            self.projections_are_initialized = True;

            (self.fluid_solver).AssignInitialStepValues();
            print "adsdfsafafds"
            (self.fluid_solver.IterativeSolve())
            print "qqqq"
            
        else:
            if(self.predictor_corrector == False):
                self.PureFractionalStepFluidSolutionStep()
            else:
                self.PredictorCorrectorFluidSolutionStep()
       
#        self.step = self.step + 1
        self.substep = self.substep + 1



##        (self.fluid_solver).Clear()
#        if(self.substep == self.mesh_every_nsteps):
#            (self.fluid_solver).Clear()

    ######################################################################
    def PureFractionalStepFluidSolutionStep(self):
        print "PureFractionalStepFluidSolutionStep beginning"
        self.SolutionStep1()
        print "PureFractionalStepFluidSolutionStep after step1"
        
        (self.fluid_solver).SolveStep2();
        (self.fluid_solver).ActOnLonelyNodes();
        if(self.uzawa_term_is_active == True):
            (self.fluid_solver).SolveStep2_Mp();
        (self.fluid_solver).SolveStep3();

        if( self.displacement_correction_afterstep4 == True):
            (self.VariableUtils).SaveVectorVar(FRACT_VEL,VAUX,(self.model_part).Nodes )
            (self.fluid_solver).SolveStep4();
            (self.lagrangian_tools).CalculateFinalDisplacementCorrection((self.model_part))
            (self.MeshMover).Execute();
        else:
            (self.fluid_solver).SolveStep4();        


    
            
    ######################################################################
    #ALE solution of the fluid on the moving domain
    def PredictorCorrectorFluidSolutionStep(self):

        print "PredictorCorrectorFluidSolutionStep beginning"
        it = 1
        pressure_is_converged = False
        while(pressure_is_converged == False and it <= self.max_press_its):
            
            self.SolutionStep1()

            dp = (self.fluid_solver).SolveStep2();
            (self.fluid_solver).ActOnLonelyNodes();
            if(self.uzawa_term_is_active == True):
                (self.fluid_solver).SolveStep2_Mp();
            (self.fluid_solver).SolveStep3();

            if( self.displacement_correction_afterstep4 == True):
                (self.VariableUtils).SaveVectorVar(FRACT_VEL,VAUX,(self.model_part).Nodes )
                (self.fluid_solver).SolveStep4();
                (self.lagrangian_tools).CalculateFinalDisplacementCorrection((self.model_part))
                (self.MeshMover).Execute();
            else:
                (self.fluid_solver).SolveStep4();
                
            p_norm = (self.fluid_solver).SavePressureIteration()

            it = it + 1

            ratio = 0.0
            if( p_norm > 1e-10):
                ratio = dp/p_norm
            print "********************* ",it, " ******", ratio , " **************"
            if(ratio < self.press_toll):
                pressure_is_converged = True
       


    ######################################################################


    ######################################################################
    #nodes in the inlet introduced and ALE acceleration calculated
    def PostSolutionStep(self):
        self.step = self.step + 1
        print "doing nothing in the post step"
        ##lagrangian fluid step
        #(self.LagrangianInlet).Execute();

        ##calculate the spatial (ALE) acceleration
##        if(self.domain_size == 2):
##            (self.PfemUtils).CalculateSpatialALEAcceleration2D(self.model_part)
##        else:
##            (self.PfemUtils).CalculateSpatialALEAcceleration3D(self.model_part)



    ######################################################################
    def OutputStep(self,time,gid_io):
##        if(time >= self.next_output_time):
        if(time >= self.next_output_time):
            self.next_output_time = self.next_output_time + self.output_time_increment;

            #writing mesh 
            gid_io.InitializeMesh( time );
            gid_io.WriteNodeMesh((self.model_part).GetMesh());
            gid_io.WriteMesh((self.model_part).GetMesh());
            gid_io.FinalizeMesh();

            gid_io.InitializeResults(time, (self.model_part).GetMesh());
            
            gid_io.WriteNodalResults(VELOCITY, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(PRESSURE, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(IS_FREE_SURFACE, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(IS_BOUNDARY, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(IS_STRUCTURE, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(VISCOSITY, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(DENSITY, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(IS_FLUID, (self.model_part).Nodes, time, 0);

            gid_io.Flush()
            gid_io.FinalizeResults()
            
            
##            file_name = self.output_name
##            file_name = file_name + str(time)
##
##            print "before changing name"
##            gid_io.ChangeOutputName(file_name,GiDPostMode.GiD_PostBinary);
##            print "name changed"
##            my_mesh = (self.model_part).GetMesh()
##            print my_mesh
##            print self.model_part
##            gid_io.WriteMesh((self.model_part).GetMesh(),self.domain_size,GiDPostMode.GiD_PostBinary);
##            print "mesh written"
####            gid_io.WriteNodalResults(NODAL_AREA, (self.model_part).Nodes, time, 0);
####            gid_io.WriteNodalResults(NORMAL, (self.model_part).Nodes, time, 0);
####            gid_io.WriteNodalResults(NORMAL_TO_WALL, (self.model_part).Nodes, time, 0);
##            gid_io.WriteNodalResults(PRESSURE, (self.model_part).Nodes, time, 0);
##            gid_io.WriteNodalResults(IS_FREE_SURFACE, (self.model_part).Nodes, time, 0);
##            gid_io.WriteNodalResults(IS_BOUNDARY, (self.model_part).Nodes, time, 0);
##            gid_io.WriteNodalResults(IS_STRUCTURE, (self.model_part).Nodes, time, 0);
####  #       gid_io.WriteNodalResults(IS_LAGRANGIAN_INLET, (self.model_part).Nodes, time, 0);
####          gid_io.WriteNodalResults(NODAL_H, (self.model_part).Nodes, time, 0);
####          gid_io.WriteNodalResults(IS_FLUID, (self.model_part).Nodes, time, 0);
##            gid_io.WriteNodalResults(VELOCITY, (self.model_part).Nodes, time, 0);
####            gid_io.WriteNodalResults(MESH_VELOCITY, (self.model_part).Nodes, time, 0);
####            gid_io.WriteNodalResults(ACCELERATION, (self.model_part).Nodes, time, 0);
####            gid_io.WriteNodalResults(BODY_FORCE, (self.model_part).Nodes, time, 0);
####            gid_io.WriteNodalResults(PRESS_PROJ, (self.model_part).Nodes, time, 0);
####            gid_io.WriteNodalResults(CONV_PROJ, (self.model_part).Nodes, time, 0);
##
##            print "results finished"
##            gid_io.Flush()
####                gid_io.CloseResultFile();

    
    ######################################################################
    #in this step the node position is predicted, the mesh is regenerated
    #and the boundary conditions are imposed
    def EulerianLagrangianPredictionStep(self,time):

        print self.step
        if(self.step > 2):
            print "*************** eulerian-lagrangian prediction ***************"  
            domain_size = self.domain_size
                  
            #calculate normals
            conditions = (self.model_part).Conditions
            (self.NormalTools).CalculateOnSimplex(conditions,domain_size)
            (self.NormalToWall).CalculateNormalToWall(conditions,domain_size)

            print "ln566"
            #performing a first order prediction of the fluid displacement
            #without moving the mesh
            (self.PfemUtils).Predict(self.model_part)

            print "ln571"
            #lagrangian improvement of the displacements
            if(self.step > 2):
                #copying the velocity to the mesh_velocity to set to zero the convective part)
                (self.utilities).CopyVectorVar(VELOCITY,MESH_VELOCITY,(self.model_part).Nodes);
                self.LagrangianCorrection()

            print "ln522"
            (self.MeshMover).Execute();

            print "ln525"

            #move free nodes and treat nodes close to the boundary
            (self.PfemUtils).MoveLonelyNodes(self.model_part)
            print "ln529"
            (self.ActOnWalls).Execute();

            print "ln 532"
            (self.MeshMover).Execute();

            print "ln529"

             #smooth the mesh position
            if( self.smooth == True):
                (self.CoordinateSmoother).Execute();
                (self.MeshMover).Execute();           

            #remesh
            self.Remesh()
            print self.model_part

            self.CalculateProjections()

    ######################################################################
    def CalculateProjections(self):
        (self.PfemUtils).CalculateNodalMass(self.model_part,self.domain_size);
        (self.fluid_solver).SolveStep3()
        print "Projections were calculated"


##    ######################################################################
##    #in this step the node position is predicted, the mesh is regenerated
##    #and the boundary conditions are imposed
##    def PredictionStep(self,time):
##        domain_size = self.domain_size
##              
##        #calculate normals
##        conditions = (self.model_part).Conditions
##        (self.NormalTools).CalculateOnSimplex(conditions,domain_size)
##        (self.NormalToWall).CalculateNormalToWall(conditions,domain_size)
##
##      #performing a first order prediction of the fluid displacement
##        print "****************************" , self.prediction_order
##        if( self.prediction_order == 12):
##            (self.PfemUtils).Predict(self.model_part)
##            (self.utilities).SetToZero_VectorVar(MESH_VELOCITY,(self.model_part).Nodes);
##            if(self.step > 2):
##                print "eulerian second order prediction"
##                self.LagrangianCorrection()
##        elif( self.prediction_order == 2):
##            print "second order prediction"
##            #(self.PfemUtils).QuasiLagrangianMove(self.model_part)
##            if(self.step > 2):
##                ##
##                (self.PfemUtils).Predict(self.model_part)
##                (self.MeshMover).Execute();
##                
##                ##calculating projections on the updated mesh position
##                (self.utilities).CopyVectorVar(VELOCITY,FRACT_VEL,(self.model_part).Nodes);
##                (self.PfemUtils).CalculateNodalArea(self.model_part,self.domain_size);
##                (self.fluid_solver).ActOnLonelyNodes();
##                (self.fluid_solver).SolveStep3()
##        
##                if(self.step > 2):
##                    self.LagrangianCorrection()
##        elif( self.prediction_order == 1):
##            print "first order prediction"
##            (self.PfemUtils).Predict(self.model_part)
##        else:
##            print "zeroth order prediction - no motion!!"
##
##        (self.MeshMover).Execute();
##            
##        (self.PfemUtils).MoveLonelyNodes(self.model_part)
##        (self.MeshMover).Execute();
##
##
##        ##ensure that no node gets too close to the walls
##        (self.ActOnWalls).Execute();
##
##        ##move the mesh
##        (self.MeshMover).Execute();
##
##        ##smooth the final position of the nodes to homogeneize the mesh distribution
##        (self.CoordinateSmoother).Execute();
##
##        ##move the mesh
##        (self.MeshMover).Execute();
##            
##        [inverted_elements,vol] = self.CheckForInvertedElements()
##
##        ##calculating projections on the updated mesh
##        if(self.step > 2):
##            (self.utilities).CopyVectorVar(VELOCITY,FRACT_VEL,(self.model_part).Nodes);
##            (self.PfemUtils).CalculateNodalArea(self.model_part,self.domain_size);
##            (self.fluid_solver).ActOnLonelyNodes();
##            (self.fluid_solver).SolveStep3()
##        
##        print "pfem_solver_ale line 262"
##        print self.model_part
##        
##        #regenerate the mesh
##        print self.substep
##        self.Remesh()


        
        
##
##        if(self.correct_volume == True):
##            #calculate the new volume (after moving and alpha shape)
##            [inverted_elements,vol] = self.CheckForInvertedElements()
##
##            #correct the volume to cure the volume loss
##            (self.VolumeCorrectionTools).CorrectVolume(self.originalvolume, vol, (self.model_part).Nodes)
##            (self.MeshMover).Execute();        


