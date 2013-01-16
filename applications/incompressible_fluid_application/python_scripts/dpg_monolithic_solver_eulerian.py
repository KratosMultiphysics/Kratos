#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

import levelset_solver

#settings for the convection solver
distance_settings = ConvectionDiffusionSettings()
distance_settings.SetUnknownVariable(DISTANCE)
distance_settings.SetConvectionVariable(VELOCITY)
distance_settings.SetMeshVelocityVariable(MESH_VELOCITY)
#distance_settings.SetVolumeSourceVariable(HEAT_FLUX)
#distance_settings.SetDiffusionVariable(ARRHENIUSAUX)
#distance_settings.SetDensityVariable(AUX_INDEX) # For level set solver rho and C are assigned to 1

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE);
    model_part.AddNodalSolutionStepVariable(IS_FLUID);
    model_part.AddNodalSolutionStepVariable(IS_POROUS);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(POROSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY_AIR);
    model_part.AddNodalSolutionStepVariable(AIR_SOUND_VELOCITY);
    model_part.AddNodalSolutionStepVariable(SOUND_VELOCITY);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(ADVPROJ);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);
    model_part.AddNodalSolutionStepVariable(THAWONE);
    model_part.AddNodalSolutionStepVariable(THAWTWO); 
    model_part.AddNodalSolutionStepVariable(REACTION); 
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_DT);
    model_part.AddNodalSolutionStepVariable(ARRHENIUS); 
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    model_part.AddNodalSolutionStepVariable(NORMAL);
    model_part.AddNodalSolutionStepVariable(DISTANCE);
    model_part.AddNodalSolutionStepVariable(IS_SLIP);
    
    #variables needed for the distance solver
    levelset_solver.AddVariables(model_part,distance_settings)

    print "variables for the MONOLITHIC_SOLVER_EULERIAN added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(VELOCITY_X,REACTION_X);
        node.AddDof(VELOCITY_Y,REACTION_Y);
        node.AddDof(VELOCITY_Z,REACTION_Z);
        node.AddDof(PRESSURE,REACTION_WATER_PRESSURE);
	node.AddDof(AIR_PRESSURE,REACTION_AIR_PRESSURE);
	
    levelset_solver.AddDofs(model_part,distance_settings)        
    print "dofs for the monolithic solver added correctly"

class MonolithicSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.domain_size = domain_size

        self.alpha = -0.3
        self.move_mesh_strategy = 0
        self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeDPGEnriched( self.alpha,self.move_mesh_strategy,self.domain_size )
        #self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakScheme( self.alpha,self.move_mesh_strategy )
        
        #definition of the solvers
        #self.linear_solver =  SkylineLUFactorizationSolver()
##        self.linear_solver =SuperLUSolver()
##        self.linear_solver = MKLPardisoSolver()

        #pPrecond = DiagonalPreconditioner()
##        pPrecond = ILU0Preconditioner()
        #self.linear_solver =  BICGSTABSolver(1e-6, 5000,pPrecond)
        
	#gmres_size = 30
	#ilu_level_of_fill = 2
	#tol = 1e-5
	#verbosity = 0
        #self.linear_solver = PastixSolver(tol,gmres_size,ilu_level_of_fill,verbosity,False)         
        #self.linear_solver = PastixSolver(verbosity,False)
        
        #new solvers
	gmres_size = 50
	tol = 1e-7
	verbosity = 0
	self.linear_solver = AMGCLSolver(AMGCLSmoother.DAMPED_JACOBI,AMGCLIterativeSolverType.BICGSTAB,tol,200,verbosity,gmres_size)         

        #definition of the convergence criteria
        self.rel_vel_tol = 1e-5
        self.abs_vel_tol = 1e-7
        self.rel_pres_tol = 1e-5
        self.abs_pres_tol = 1e-7

       # self.conv_criteria = UPCriteria(1e-12,1e-14,1e-15,1e-17)
        #self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, 0.001);

        self.dynamic_tau_levelset = 0.01
        self.dynamic_tau_fluid = 1.0        
        self.oss_switch  = 0

        #non newtonian setting
        self.regularization_coef = 1000
        
        self.max_iter = 30
                            
        #default settings
        self.echo_level = 0
        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False
    
##        print "Construction monolithic solver finished"
        
        # Creat Lavel_set solver
        #construct the model part 
        if(domain_size == 2):
            raise "error, still not implemented in 2D"
            conv_elem = "SUPGConv2D"
            conv_cond = "Condition2D"
        else:
            conv_elem = "SUPGConv3D"
            conv_cond = "Condition3D"
        self.level_set_model_part = ModelPart("level_set_model_part")
        self.conv_generator = ConnectivityPreserveModeler()
        (self.conv_generator).GenerateModelPart(self.model_part,self.level_set_model_part,conv_elem,conv_cond)        
        #(ParallelFillCommunicator(self.level_set_model_part)).Execute();   

        #constructing the convection solver for the distance
        self.level_set_solver = levelset_solver.Solver(self.level_set_model_part,domain_size,distance_settings)
        self.level_set_solver.max_iterations = 8
        
        ################################################
        #properties of the two fluids
        self.rho1 = 2400.0 #applied on the negative part of the domain 1000.0
        self.conductivity1 = 1.0
        
        self.rho2 = 1.0 #applied to the positive part of the domain#1.0
        self.conductivity2 = 1.0 
        
	self.mu   = 3.0e-3
	self.divergence_clearance_performed = False
        ################################################ 
        
        #Distance utilities
        
         ################################################        
        if(self.domain_size == 2):
            self.redistance_utils = ParallelDistanceCalculator2D()
        else:
            self.redistance_utils = ParallelDistanceCalculator3D()

        self.max_levels = 5
	self.redistance_frequency = 1
        self.max_edge_size = self.redistance_utils.FindMaximumEdgeSize(self.level_set_model_part)
        self.max_distance = self.max_edge_size * 5.0;


        self.max_ns_iterations = 8
	self.internal_step_counter = 1  
	
	###Slip condition
	self.use_slip_conditions = False
    #######################################################################	
    def ApplyFluidProperties(self):
        #apply density
        mu1 = 1.0*self.mu/self.rho1
        #mu1 = self.mu
        #mu2 = 0.01*self.mu/self.rho2
        mu2 = mu1

        for node in self.model_part.Nodes:
            dist = node.GetSolutionStepValue(DISTANCE)
            if(dist < 0):
                node.SetSolutionStepValue(DENSITY,0,self.rho1)
                node.SetSolutionStepValue(VISCOSITY,0,mu1)
            else:
                node.SetSolutionStepValue(DENSITY,0,self.rho2)
                node.SetSolutionStepValue(VISCOSITY,0,mu2)	
    #######################################################################
    def Initialize(self):
        #creating the solution strategy
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol,self.abs_vel_tol,\
                                           self.rel_pres_tol,self.abs_pres_tol)
##        self.conv_criteria = UPCriteria(self.rel_vel_tol,self.abs_vel_tol,
##                                        self.rel_pres_tol,self.abs_pres_tol)
	builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        #self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,builder_and_solver,self.max_iter,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)   

        self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,builder_and_solver,self.max_iter,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)   
        (self.solver).SetEchoLevel(self.echo_level)
        print ">>>>>>>>>>>>>>>", self.oss_switch
        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau_fluid);
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch );
        self.model_part.ProcessInfo.SetValue(M, self.regularization_coef );

##        print "Initialization monolithic solver finished"
        # LEvel_set solver initialization
        self.level_set_solver.dynamic_tau =self.dynamic_tau_levelset
        self.redistance_utils.CalculateDistances(self.model_part,DISTANCE,NODAL_AREA,self.max_levels,self.max_distance)
        self.level_set_solver.Initialize()

        self.ApplyFluidProperties()
        
        # assigne to 1 density(AUX_INDEX) and specific_heat and to 0 conductivitu for the level_set solver
        #for node in self.model_part.Nodes:
            #node.SetSolutionStepValue(AUX_INDEX,0,1.0) 
            #node.SetSolutionStepValue(ARRHENIUSAUX,0,0.0)            
            #node.SetSolutionStepValue(SPECIFIC_HEAT,0,1.0)
            
        ############### saving inlet nodes
        #self.inlet_nodes = []
	#for node in self.model_part.Nodes:
	  #if(node.IsFixed(DISTANCE)):
	    #self.inlet_nodes.append(node);

        self.next_redistance = self.redistance_frequency
        ###         ###         ###        
        ### FOR SLIP
        ###         ###         ###
        # Manullay assign!
        for cond in self.model_part.Conditions:
            cond.SetValue(IS_STRUCTURE,1.0)
        # if we use slip conditions, calculate normals on the boundary
        if (self.use_slip_conditions == True):
	    (FindConditionsNeighboursProcess(self.model_part, 3, 20)).ClearNeighbours()          
	    (FindConditionsNeighboursProcess(self.model_part, 3, 20)).Execute()	  
            self.normal_util = NormalCalculationUtils()
            self.normal_util.CalculateOnSimplex(self.model_part,self.domain_size,IS_STRUCTURE,0.0,35.0)#,0.0,35.0
            
            #for node in self.model_part.Nodes:
		#if (node.GetSolutionStepValue(IS_SLIP) > 0.0):
		   #node.SetValue(IS_STRUCTURE,1.0)
		#else:
		   #node.SetValue(IS_STRUCTURE,0.0)		  
		   #node.SetSolutionStepValue(IS_STRUCTURE,0,1.0)
		   #node.Fix(IS_STRUCTURE)

	############## saving inlet nodes
	self.inlet_nodes = []
        for cond in self.model_part.Conditions:
	  if(cond.GetValue(IS_INLET) > 0):
	    for node in cond.GetNodes():
	      self.inlet_nodes.append(node);	      
		   
    ################################################################
    ################################################################
    def WettenNodes(self):
      for node in self.inlet_nodes:
	if(node.GetSolutionStepValue(DISTANCE) > 0):
	  node.SetSolutionStepValue(DISTANCE,0,-1.0e-3); 
	  
    #######################################################################   
    def DoRedistance(self):
	#self.WettenNodes()
	
	#redistance if required
        print "beginning recalculation of distances"
        self.redistance_utils.CalculateDistances(self.model_part,DISTANCE,NODAL_AREA,self.max_levels,self.max_distance)

	print "finished recalculation of distances"    
    
     #######################################################################   
    def ConvectDistance(self):
        print "beginning convection step for the distance function"
        self.level_set_model_part.ProcessInfo = self.model_part.ProcessInfo
        (self.level_set_model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS,distance_settings)
        (self.level_set_model_part.ProcessInfo).SetValue(DYNAMIC_TAU,self.dynamic_tau_levelset)#self.dynamic_tau
	#self.WettenNodes()
	(self.level_set_solver).Solve()
        print "finished convection step for the distance function"     
     #######################################################################                 
      #######################################################################      
    def Solve(self):
 	#at the beginning of the calculations do a div clearance step  
  	if(self.divergence_clearance_performed == False):
	  for node in self.model_part.Nodes:	    
	    node.SetSolutionStepValue(DISTANCE,1,node.GetSolutionStepValue(DISTANCE))
	  self.divergence_clearance_performed = True    

        ##recompute distance function as needed
        #if(self.internal_step_counter >= self.next_redistance):
	  #self.DoRedistance()
	  #for node in self.model_part.Nodes:	    
	    #node.SetSolutionStepValue(DISTANCE,1,node.GetSolutionStepValue(DISTANCE))
	  #self.next_redistance = self.internal_step_counter + self.redistance_frequency

	#convect distance function
        self.ConvectDistance()
        #recompute distance function as needed
        if(self.internal_step_counter >= self.next_redistance):
	  self.DoRedistance()
	  self.next_redistance = self.internal_step_counter + self.redistance_frequency	  

	self.ApplyFluidProperties()
        #Recompute normals if necessary
	if(self.ReformDofSetAtEachStep == True):
           if self.use_slip_conditions == True:
	      self.normal_util.CalculateOnSimplex(self.model_part,self.domain_size,IS_STRUCTURE,0.0,35.0)#,0.0,35.0
        
        (self.solver).Solve()
##        print "solving step monolithic solver finished"
        self.internal_step_counter += 1       

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
    
    ########################################################################

        




