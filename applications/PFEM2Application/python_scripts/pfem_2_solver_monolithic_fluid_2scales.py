from KratosMultiphysics import *
from KratosMultiphysics.PFEM2Application import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.OpenCLApplication import *
#from KratosMultiphysics.GPUSolversApplication import *

from math import sqrt
import time as timer
#dir(pure_diffusion_application)

def AddVariables(model_part,linea_model_part):
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(DISTANCE);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(FRACT_VEL);
    model_part.AddNodalSolutionStepVariable(G_VALUE);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(YP);

    model_part.AddNodalSolutionStepVariable(RHS);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(NORMAL);	
    #model_part.AddNodalSolutionStepVariable(DENSITY);	

    model_part.AddNodalSolutionStepVariable(DELTA_VELOCITY) 
    model_part.AddNodalSolutionStepVariable(MEAN_SIZE) 
    model_part.AddNodalSolutionStepVariable(NODAL_AREA) 
    model_part.AddNodalSolutionStepVariable(NODAL_MASS) 
    model_part.AddNodalSolutionStepVariable(MASS)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(TAU)
    model_part.AddNodalSolutionStepVariable(TEMPERATURE)



def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);
        node.AddDof(DISTANCE);
        node.AddDof(TEMPERATURE);


    print "DoFs for the Poisson solver added correctly"

class PFEM2Solver:
    #######################################################################
    #def __init__(self,model_part,linea_model_part,domain_size):
    def __init__(self,model_part,model_part_topo,domain_size):
        self.model_part = model_part
        self.model_part_topo = model_part_topo
        
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        #definition of the solvers
	gmres_size = 50
	tol = 1e-5
	verbosity = 1
	pDiagPrecond = DiagonalPreconditioner()
	#self.monolitic_linear_solver= GPUBICGSTABSolverWithDiagonalPreconditioner(tol, 5000)
	#self.monolitic_linear_solver= ViennaCLSolver(tol,3000,OpenCLPrecision.Double,OpenCLSolverType.CG,OpenCLPreconditionerType.AMG_DAMPED_JACOBI) 
	#self.monolitic_linear_solver = BICGSTABSolver(1e-5, 1000,pDiagPrecond) # SkylineLUFactorizationSolver() 
	#self.monolitic_linear_solver = BICGSTABSolver(1e-5, 1000,pDiagPrecond) # SkylineLUFactorizationSolver() 
	self.monolitic_linear_solver =  AMGCLSolver(AMGCLSmoother.ILU0,AMGCLIterativeSolverType.CG,tol,500,verbosity,gmres_size)      #BICGSTABSolver(1e-7, 5000) # SkylineLUFactorizationSolver() 
	
	#self.viscosity_linear_solver = AMGCLSolver(AMGCLSmoother.DAMPED_JACOBI,AMGCLIterativeSolverType.CG,tol,200,verbosity,gmres_size)      #BICGSTABSolver(1e-7, 5000) # SkylineLUFactorizationSolver() #
	#self.pressure_linear_solver = DeflatedCGSolver(1e-5, 3000, True,1000)
	#self.viscosity_linear_solver = DeflatedCGSolver(1e-6, 3000, True,1000)
	self.temperature_linear_solver =  BICGSTABSolver(1e-4, 5000,pDiagPrecond) #
	#self.temperature_linear_solver =  AMGCLSolver(AMGCLSmoother.DAMPED_JACOBI,AMGCLIterativeSolverType.CG,tol,200,verbosity,gmres_size)      #BICGSTABSolver(1e-7, 5000) # SkylineLUFactorizationSolver() #
        #self.pressure_linear_solver =  SkylineLUFactorizationSolver()  #we set the type of solver that we want 
        self.conv_criteria = DisplacementCriteria(1e-9,1e-15)  #tolerance for the solver 
        
        self.domain_size = domain_size
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part)
        (self.neighbour_search).Execute()
	self.neighbour_elements_search= FindElementalNeighboursProcess(model_part,domain_size,number_of_avg_elems)
	(self.neighbour_elements_search).Execute()
	self.neighbour_elements_search_topo= FindElementalNeighboursProcess(model_part_topo,domain_size,number_of_avg_elems)
	(self.neighbour_elements_search_topo).Execute()
        ##calculate normals
        self.normal_tools = BodyNormalCalculationUtils()
        
        
    #######################################################################
    #######################################################################
    def Initialize(self):
        #creating the solution strategy
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False	
        MoveMeshFlag = False
        pDiagPrecond = DiagonalPreconditioner()
        #build the edge data structure
	#if self.domain_size==2:
	#	self.matrix_container = MatrixContainer2D()
	#else: 
	#	self.matrix_container = MatrixContainer3D()
	#maximum_number_of_particles= 10*self.domain_size
	maximum_number_of_particles= 8*self.domain_size

	self.ExplicitStrategy=PFEM2_Explicit_Strategy(self.model_part,self.domain_size, MoveMeshFlag)

	self.VariableUtils = VariableUtils()
	
	if self.domain_size==2:
		self.moveparticles = MoveParticleUtilityDiffFluidOnly2D(self.model_part,maximum_number_of_particles)
	else:
		self.moveparticles = MoveParticleUtilityDiffFluidOnly3D(self.model_part,maximum_number_of_particles)

	print "self.domain_size = ", self.domain_size
	if self.domain_size==2:
		self.calculatewatervolume = CalculateWaterFraction2D(self.model_part)
	else:	
		self.calculatewatervolume = CalculateWaterFraction3D(self.model_part)
	print "hola0"
	self.moveparticles.MountBinDiff()
        #self.matrix_container.ConstructCSRVector(self.model_part)
        #self.matrix_container.BuildCSRData(self.model_part)
	self.water_volume=0.0  #we initialize it at zero
	self.water_initial_volume=0.0 #we initialize it at zero
	self.water_initial_volume_flag=True #we initialize it at zero
	self.mass_correction_factor=0.0

        

	self.normal_tools.CalculateBodyNormals(self.model_part,self.domain_size);  
	condition_number=1
	'''
	if self.domain_size==2:
		self.addBC = AddFixedVelocityCondition2D(self.model_part)
	else:
		self.addBC = AddFixedVelocityCondition3D(self.model_part)

	(self.addBC).AddThem()
	'''	
	#RICC TOOLS
	self.convection_solver = PureConvectionUtilities2D();
	self.convection_solver.ConstructSystem(self.model_part,DISTANCE,VELOCITY,FORCE);
	if self.domain_size==2:
		self.distance_utils = SignedDistanceCalculationUtils2D()
	else:
		self.distance_utils = SignedDistanceCalculationUtils3D()
	self.redistance_step = 0

	print "hola1"
	import strategy_python #implicit solver

	self.monolitic_solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.monolitic_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
	self.temperature_solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.temperature_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
	#pDiagPrecond = DiagonalPreconditioner()
	self.linear_solver =  BICGSTABSolver(1e-7, 5000,pDiagPrecond)
	#(self.moveparticles).ReplaceParticlesVelocityAndDistance();
	self.calculateinitialdrag=0.0
	self.streamlineintegration=0.0
	self.accelerateparticles=0.0
	self.navierstokes=0.0
	self.lagrangiantoeulerian=0.0
	self.reseed=0.0
	self.prereseed=0.0
	self.erasing=0.0
	self.total=0.0
	self.nodaltasks=0.0
	self.thermal=0.0
	self.implicitviscosity=0.0
	print "hola2"
        
	overwrite_particle_data=False
	intial_offset=Vector(3);
	intial_offset[0]=0.0
	intial_offset[1]=0.0
	intial_offset[2]=0.0
	(self.moveparticles).IntializeTransferTool(self.model_part_topo, intial_offset,overwrite_particle_data);
      
                 
    #######################################################################   
    def Solve(self):
		combustion_problem=False;
		


		add_gravity=False #in the monolitic solver we do not add the gravity, it is added directly in the implicit monolitic system
		#add_gravity=True #in the monolitic solver we do not add the gravity, it is added directly in the implicit monolitic system

		transfer_pressure=False
		temperature_implicit_calculation=False;
		viscosity_streamline_integrate=False; #True means explicit integration! False means we will have to solve the fractional velocity system implicetely later, once information has been passed to mesh
		#pressure_gradient_integrate=False; #we use the pressure of the previous step. not reccomended for 2 fluid flows. REMOVED FROM THE CODE	
		calculate_viscous_contribution_after_everything=False;		
		#implicit_velocity_correction=False;	# as expected this can not be active the implicit correction is done:
		#if (calculate_viscous_contribution_after_everything):
		#	implicit_velocity_correction=False;
		if combustion_problem:
			temperature_implicit_calculation=True;
		t1 = timer.time()

		#calculating RHS by viscous forces using information of the previous time step:	
		#self.CalculateExplicitViscosityContribution();
		t2 = timer.time()
		self.calculateinitialdrag = self.calculateinitialdrag + t2-t1

		#in order to define a reasonable number of substeps we calculate a mean courant in each element
		(self.moveparticles).CalculateVelOverElemSize();
		t2a = timer.time()

		#streamline integration:
		(self.moveparticles).MoveParticlesDiff(viscosity_streamline_integrate,add_gravity);	
		t3 = timer.time()
		self.streamlineintegration = self.streamlineintegration + t3-t2a
		t3a = timer.time()

		#Reseeding using streamlines in inverse way in case there are too few particles. Not very accurate since particles follow streamlines, no matter if it's water or air.
		#(for 1 fluid should be accurate though
		pre_minimum_number_of_particles=1#self.domain_size;
		added_displacement=Vector(3);
		added_displacement[0]=-0.0
		added_displacement[1]=-0.02
		added_displacement[2]=-0.0
		(self.moveparticles).PreReseedUsingTopographicDomain(pre_minimum_number_of_particles,added_displacement)
		t4 = timer.time()
		self.reseed = self.reseed + t4-t3a
		self.prereseed = self.prereseed + t4-t3a

		#self.distance_utils.CalculateDistances(self.model_part,DISTANCE,0.1)

		#for node in self.model_part.Nodes:
		#	node.SetSolutionStepValue(TEMP_CONV_PROJ,node.GetSolutionStepValue(DISTANCE))

		#transfering data from the particles to the mesh:
		(self.moveparticles).TransferLagrangianToEulerian(transfer_pressure);

		#for node in self.model_part.Nodes:
		#	node.SetSolutionStepValue(MESH_VELOCITY,node.GetSolutionStepValue(VELOCITY))

		#for node in self.model_part.Nodes:
		#	node.SetSolutionStepValue(DISTANCE,node.GetSolutionStepValue(TEMP_CONV_PROJ))
		#self.convection_solver.CalculateProjection(self.model_part,DISTANCE,NODAL_AREA,VELOCITY,FORCE,TEMP_CONV_PROJ);
		#self.convection_solver.ConvectScalarVar(self.model_part,self.linear_solver,DISTANCE,VELOCITY,FORCE,TEMP_CONV_PROJ,1);
		
		print "hola1"

		#value  = elem_size*0.02


		t5 = timer.time()
		self.lagrangiantoeulerian = self.lagrangiantoeulerian + t5-t4
		print "hola2"
		#Flagging splitted elements (the ones who need special integration)
		#if combustion_problem:
		#	(self.moveparticles).FlagSplittedElementsAndTheirNodesForCombustion();
		#else:
		(self.moveparticles).FlagSplittedElementsAndTheirNodes();
		print "hola4"
		#if (implicit_velocity_correction==False): #already solved in the streamlines, no need to do it now.
		#	(self.moveparticles).MoveDataFromMESH_VELOCITYtoFRACT_VEL(); #insted of first step,. when moving through streamlines the first step is already performed
		#else: #we did not do it. doing it explicitely
		#	(self.moveparticles).CalculateFirstStep();
			#self.model_part.ProcessInfo.SetValue(NL_ITERATION_NUMBER, 0) #1= semiimplicit, =0 fully implciit
			#(self.solver).Solve()
		
		t6 = timer.time()

		if (temperature_implicit_calculation):
			self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 5)
			print "calculating temperature implicitely"	
			(self.temperature_solver).Solve()
		t7=timer.time()
		self.thermal = self.thermal + t7-t6
		print "hola5"
		
		(self.VariableUtils).CopyVectorVar(MESH_VELOCITY,VELOCITY,self.model_part.Nodes)
		#self.CalculateExplicitPressureViscousCorrection()
		#pressure iterations
		
		
		full_reset=False;
		#self.CalculateMassMatrix();
		#(self.moveparticles).ResetBoundaryConditions(full_reset) 
		
		#implicit everything
		print "hola6"
		full_reset=True;
		(self.moveparticles).ResetBoundaryConditions(full_reset) 

		(self.monolitic_solver).Solve() #implicit resolution of the pressure system. All the other tasks are explicit
		print "hola7"
		full_reset=True;
		(self.moveparticles).ResetBoundaryConditions(full_reset) 
		#delta_velocity= Velocity(final) - MeshVelocity(from the particles), so we add to the particles the correction done in the mesh.
		(self.moveparticles).CalculateDeltaVelocity();		
		t11 = timer.time()
		self.navierstokes = self.navierstokes + t11-t6
		#transfering the information to the mesh:
		modify_particle_pressure=False
		(self.moveparticles).AccelerateParticlesWithoutMovingUsingDeltaVelocity(add_gravity,modify_particle_pressure);
		t12 = timer.time()
		self.accelerateparticles = self.accelerateparticles + t12-t11

		#reseeding in elements that have few particles to avoid having problems in next iterations:
		post_minimum_number_of_particles=self.domain_size*2;
		self.water_volume = (self.calculatewatervolume).Calculate()
		if (self.water_initial_volume_flag): 
			print "calculated water volume for the first time"
			self.water_initial_volume=self.water_volume
			self.water_initial_volume_flag=False
		print self.water_volume
		water_fraction= self.water_volume/(self.water_initial_volume+1e-9)
		self.mass_correction_factor = (1.0 - water_fraction) * 100.0 * 0.0
		if self.mass_correction_factor>0.5:
			 self.mass_correction_factor
		print "mass correction factor: ", self.mass_correction_factor
		print "current mass loss is : " , (1.0 - water_fraction) * 100.0 , " % "
		bounding_box_lower_corner=Vector(3);
		bounding_box_lower_corner[0]=0.0
		bounding_box_lower_corner[1]=0.1
		bounding_box_lower_corner[2]=0.0
		bounding_box_upper_corner=Vector(3);
		bounding_box_upper_corner[0]=1.0
		bounding_box_upper_corner[1]=2.0
		bounding_box_upper_corner[2]=0.0
		(self.moveparticles).PostReseedOnlyInBoundingBox(post_minimum_number_of_particles,self.mass_correction_factor,bounding_box_lower_corner,bounding_box_upper_corner);		
		t13 = timer.time()
		self.reseed=self.reseed+ t13-t12
		#self.nodaltasks = self.nodaltasks + t11-t9

		#print ". erasing  = ", t14-t13
		self.total = self.total + t13-t1 
		
		if combustion_problem==False:
			print "----------TIMES----------"
			print "self.calculateinitialdrag  " , self.calculateinitialdrag , "in % = ",  100.0*(self.calculateinitialdrag)/(self.total)
			print "self.streamlineintegration " ,  self.streamlineintegration , "in % = ", 100.0*(self.streamlineintegration)/(self.total)
			print "self.accelerateparticles " ,  self.accelerateparticles , "in % = ", 100.0*(self.accelerateparticles)/(self.total)
			print "self.navierstokes " ,  self.navierstokes , "in % = ", 100.0*(self.navierstokes)/(self.total)
			print "self.nodaltasks ", self.nodaltasks , (self.nodaltasks)/(self.total)
			print "self.lagrangiantoeulerian " ,  self.lagrangiantoeulerian , "in % = ", 100.0*(self.lagrangiantoeulerian)/(self.total)
			print "self.reseed " ,  self.reseed , "in % = ", 100.0*(self.reseed)/(self.total)
			print "self.prereseed " ,  self.prereseed , "in % = ", 100.0*(self.prereseed)/(self.total)
			print "TOTAL ----- " ,  self.total
			print "THIS TIME STEP = " , t13-t1 
		else:
			print "----------TIMES----------"
			print "self.calculateinitialdrag  " , self.calculateinitialdrag , "in % = ",  100.0*(self.calculateinitialdrag)/(self.total)
			print "self.streamlineintegration " ,  self.streamlineintegration , "in % = ", 100.0*(self.streamlineintegration)/(self.total)
			print "self.accelerateparticles " ,  self.accelerateparticles , "in % = ", 100.0*(self.accelerateparticles)/(self.total)
			pressure_iterations=self.navierstokes-self.implicitviscosity-self.thermal
			print "pressureiterations " , pressure_iterations , "in % = ", (100.0*pressure_iterations)/(self.total)
			print "self.implicitviscosity " ,  self.implicitviscosity , "in % = ", 100.0*(self.implicitviscosity)/(self.total)
			print "self.allimplicitsystems " ,  self.navierstokes , "in % = ", 100.0*(self.navierstokes)/(self.total)
			print "self.thermal " ,  self.thermal , "in % = ", 100.0*(self.thermal)/(self.total)
			print "self.nodaltasks ", self.nodaltasks , (self.nodaltasks)/(self.total)
			print "self.lagrangiantoeulerian " ,  self.lagrangiantoeulerian , "in % = ", 100.0*(self.lagrangiantoeulerian)/(self.total)
			print "self.reseed " ,  self.reseed , "in % = ", 100.0*(self.reseed)/(self.total)
			print "self.prereseed " ,  self.prereseed , "in % = ", 100.0*(self.prereseed)/(self.total)
			print "TOTAL ----- " ,  self.total
			print "THIS TIME STEP = " , t13-t1 
		

    #######################################################################   
    def CalculatePressureIteration(self, iteration_number):
		self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 2)
		self.model_part.ProcessInfo.SetValue(NL_ITERATION_NUMBER, iteration_number)
		
		(self.pressure_solver).Solve() #implicit resolution of the pressure system. All the other tasks are explicit
		#if iteration_number==1:
		#	for node in self.model_part.Nodes:
		#		press= node.GetSolutionStepValue(PRESSURE)+node.GetSolutionStepValue(PRESSUREAUX)
		#		node.SetSolutionStepValue(PRESSURE,press)
		#setting the Fract step number to 3 to calculate the pressure gradient effect on the velocity
		self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 3)
		(self.ExplicitStrategy).InitializeSolutionStep(self.model_part.ProcessInfo);
		(self.ExplicitStrategy).AssembleLoop(self.model_part.ProcessInfo);
		(self.ExplicitStrategy).FinalizeSolutionStep(self.model_part.ProcessInfo);


    #######################################################################   
    def CalculateExplicitViscosityContribution(self):
		self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 0) #explicit contribution by viscosity is defined as fract step = 0
		(self.ExplicitStrategy).InitializeSolutionStep(self.model_part.ProcessInfo);
		(self.ExplicitStrategy).AssembleLoop(self.model_part.ProcessInfo);
		(self.ExplicitStrategy).FinalizeSolutionStep(self.model_part.ProcessInfo);

    #######################################################################   

    def CalculateExplicitPressureViscousCorrection(self):
		self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 6) #explicit contribution by viscosity is defined as fract step = 0
		(self.ExplicitStrategy).InitializeSolutionStep(self.model_part.ProcessInfo);
		(self.ExplicitStrategy).AssembleLoop(self.model_part.ProcessInfo);
		(self.ExplicitStrategy).FinalizeSolutionStep(self.model_part.ProcessInfo);

    ####################################################################### 

    #######################################################################   

    def CalculateMassMatrix(self):
		self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 7) #explicit contribution by viscosity is defined as fract step = 0
		(self.ExplicitStrategy).InitializeSolutionStep(self.model_part.ProcessInfo);
		(self.ExplicitStrategy).AssembleLoop(self.model_part.ProcessInfo);
		#(self.ExplicitStrategy).FinalizeSolutionStep(self.model_part.ProcessInfo);

    ####################################################################### 

    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
