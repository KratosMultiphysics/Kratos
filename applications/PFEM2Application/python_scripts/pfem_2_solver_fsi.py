from __future__ import print_function
from KratosMultiphysics import *
from KratosMultiphysics.PFEM2Application import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.OpenCLApplication import *

from math import sqrt
import time as timer
#dir(pure_diffusion_application)

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(SOLID_PRESSURE);
    model_part.AddNodalSolutionStepVariable(SOLID_YP);
    model_part.AddNodalSolutionStepVariable(DISTANCE);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(FRACT_VEL);
    model_part.AddNodalSolutionStepVariable(G_VALUE);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(YP);
    #model_part.AddNodalSolutionStepVariable(ERASE_FLAG);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    #model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    model_part.AddNodalSolutionStepVariable(RHS);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    #model_part.AddNodalSolutionStepVariable(MESH_PRESSURE);
    model_part.AddNodalSolutionStepVariable(NORMAL);	
    #model_part.AddNodalSolutionStepVariable(DENSITY);	
    model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ);
    model_part.AddNodalSolutionStepVariable(PREVIOUS_ITERATION_PRESSURE);
    #model_part.AddNodalSolutionStepVariable(FIRST_ITERATION_PRESSURE);
    model_part.AddNodalSolutionStepVariable(DELTA_VELOCITY) 
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ_NO_RO) 
    model_part.AddNodalSolutionStepVariable(MEAN_SIZE) 
    model_part.AddNodalSolutionStepVariable(NODAL_AREA) 
    model_part.AddNodalSolutionStepVariable(NODAL_MASS) 
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(TAU)
    model_part.AddNodalSolutionStepVariable(TEMPERATURE)
    model_part.AddNodalSolutionStepVariable(TEMPERATURE_OLD_IT)
    model_part.AddNodalSolutionStepVariable(SPLIT_ELEMENT)
    model_part.AddNodalSolutionStepVariable(CORRECTED_DISTANCE)
    model_part.AddNodalSolutionStepVariable(PRESSUREAUX)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);



def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);
        node.AddDof(SOLID_PRESSURE);
        node.AddDof(DISTANCE);
        #node.AddDof(TEMPERATURE);


    print("DoFs for the Poisson solver added correctly")

class PFEM2Solver:
    #######################################################################
    #def __init__(self,model_part,linea_model_part,domain_size):
    def __init__(self,model_part,domain_size):
        self.model_part = model_part
        
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
        
        #definition of the solvers
        gmres_size = 50
        tol = 1e-6
        verbosity = 1
        smoother = AMGCLSmoother.ILU0
        self.monolitic_linear_solver =  AMGCLSolver(smoother,AMGCLIterativeSolverType.CG,tol,5000,verbosity,gmres_size)      #BICGSTABSolver(1e-7, 5000) # SkylineLUFactorizationSolver()   
        self.conv_criteria = DisplacementCriteria(1e-6,1e-12)  #tolerance for the solver 
        self.domain_size = domain_size
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part)
        (self.neighbour_search).Execute()
        self.neighbour_elements_search= FindElementalNeighboursProcess(model_part,domain_size,number_of_avg_elems)
        (self.neighbour_elements_search).Execute()
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
	maximum_number_of_particles= 12*self.domain_size

	self.ExplicitStrategy=PFEM2_Explicit_Strategy(self.model_part,self.domain_size, MoveMeshFlag)

	self.VariableUtils = VariableUtils()
	
	if self.domain_size==2:
           self.moveparticles = MoveParticleUtilityDiff2D(self.model_part,maximum_number_of_particles)
	else:
           self.moveparticles = MoveParticleUtilityDiff3D(self.model_part,maximum_number_of_particles)

	print("self.domain_size = ", self.domain_size)
	if self.domain_size==2:
            self.calculatewatervolume = CalculateWaterFraction2D(self.model_part)
	else:	
           self.calculatewatervolume = CalculateWaterFraction3D(self.model_part)

	self.moveparticles.MountBinDiff()
	self.water_volume=0.0  #we initialize it at zero
	self.water_initial_volume=0.0 #we initialize it at zero
	self.water_initial_volume_flag=True #we initialize it at zero
	self.mass_correction_factor=0.0

	self.normal_tools.CalculateBodyNormals(self.model_part,self.domain_size);  
	condition_number=1
        ''' 
	if self.domain_size==2:
           self.addBC = AddMonolithicFixedVelocityCondition2D(self.model_part)
	else:
           self.addBC = AddMonolithicFixedVelocityCondition3D(self.model_part)
        (self.addBC).AddThem()
        '''

	import strategy_python #implicit solver

	self.monolitic_solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.monolitic_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)

        #to measure time spent in solving 
	self.total=0.0
	self.implicit_tasks=0.0

        
      
                 
    #######################################################################   
    def Solve(self):
        add_gravity=False #in the monolitic solver we do not add the gravity, it is added directly in the implicit monolitic system
        viscosity_streamline_integrate=False; #True means explicit integration! False means we will have to solve the fractional velocity system implicetely later, once information has been passed to mesh

        t1 = timer.time()

        #calculating RHS by viscous forces using information of the previous time step:	
        #self.CalculateExplicitViscosityContribution();
        t2 = timer.time()

        #in order to define a reasonable number of substeps we calculate a mean courant in each element
        (self.moveparticles).CalculateVelOverElemSize();
        t2a = timer.time()

        #streamline integration:
        #(self.moveparticles).MoveParticlesDiff(viscosity_streamline_integrate,add_gravity);	
        t3 = timer.time()
        
        #Reseeding using streamlines in inverse way in case there are too few particles. Not very accurate since particles follow streamlines, no matter if it's water or air.
        #(for 1 fluid should be accurate though
        pre_minimum_number_of_particles=self.domain_size*1;
        #(self.moveparticles).PreReseed(viscosity_streamline_integrate,add_gravity,pre_minimum_number_of_particles);
        t4 = timer.time()

        transfer_pressure=True
        (self.moveparticles).TransferLagrangianToEulerian(transfer_pressure);
        (self.VariableUtils).CopyVectorVar(VELOCITY,MESH_VELOCITY,self.model_part.Nodes)        
        (self.VariableUtils).CopyScalarVar(EXTERNAL_PRESSURE,PRESSURE,self.model_part.Nodes)	    

        (self.VariableUtils).CopyVectorVar(MESH_VELOCITY,VELOCITY,self.model_part.Nodes)	
        (self.moveparticles).ResetBoundaryConditions(True) 
        (self.VariableUtils).CopyScalarVar(PRESSURE,SOLID_PRESSURE,self.model_part.Nodes)	
        (self.VariableUtils).CopyScalarVar(PRESSURE,PRESSUREAUX,self.model_part.Nodes)	
        (self.VariableUtils).CopyScalarVar(DISTANCE,CORRECTED_DISTANCE,self.model_part.Nodes)	
        
        (self.moveparticles).CopyVectorVarToPreviousTimeStep(VELOCITY,self.model_part.Nodes)
        (self.moveparticles).CopyScalarVarToPreviousTimeStep(PRESSURE,self.model_part.Nodes)

        t5 = timer.time()

        (self.moveparticles).CorrectFreeSurface();

        t6 = timer.time()

        #iterations:
        non_linear_iteration_number=1
        for i in range(0,1):
            self.model_part.ProcessInfo.SetValue(NL_ITERATION_NUMBER, non_linear_iteration_number)
            #implicit everything
            #self.CalculatePressureProjection()
            self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 20)
            (self.moveparticles).ComputeDeltaVelocityForNonLinearIteration();	 #saving velocity from previous iteration
            (self.monolitic_solver).Solve() #implicit resolution of the system. All the other tasks are explicit
            self.CalculatePressureProjection()
            full_reset=True;
            (self.moveparticles).ResetBoundaryConditions(full_reset) 
            
            (self.moveparticles).UpdateParticleStresses()	
            #finally we save the last pressure for future iterations. No need to do it for the velocity since it is saved in delta_velocity
            (self.VariableUtils).CopyScalarVar(PRESSURE,PRESSUREAUX,self.model_part.Nodes) #we save the pressure of this iteration in order to have the delta_pressure in the next one
            non_linear_iteration_number =  non_linear_iteration_number +1  
        
        (self.VariableUtils).CopyScalarVar(PRESSURE,EXTERNAL_PRESSURE,self.model_part.Nodes)	    
        t11 = timer.time()
        self.implicit_tasks = self.implicit_tasks + t11-t6
    
        #delta_velocity= Velocity(final) - MeshVelocity(from the particles), so we add to the particles the correction done in the mesh.
        (self.moveparticles).CalculateDeltaVelocity();	
        
        #transfering the information to the particles:
        (self.moveparticles).AccelerateParticlesWithoutMovingUsingDeltaVelocity(add_gravity);
        t12 = timer.time()

        #reseeding in elements that have few particles to avoid having problems in next iterations:
        post_minimum_number_of_particles=self.domain_size*2;
        #(self.moveparticles).PostReseed(post_minimum_number_of_particles,self.mass_correction_factor);        

        self.water_volume = (self.calculatewatervolume).Calculate()
        if (self.water_initial_volume_flag): 
        	print("calculated water volume for the first time")
        	self.water_initial_volume=self.water_volume
        	self.water_initial_volume_flag=False
        print(self.water_volume)
        water_fraction= self.water_volume/(self.water_initial_volume+1e-9)
        self.mass_correction_factor = (1.0 - water_fraction) * 100.0 * 0.0
        print("mass correction factor: ", self.mass_correction_factor)
        print("current mass loss is : " , (1.0 - water_fraction) * 100.0 , " % ")

        t13 = timer.time()
        self.total = self.total + t13-t1 
        
        print("----------TIMES----------")
        print("self.implicit_tasks " ,  self.implicit_tasks , "in % = ", 100.0*(self.implicit_tasks)/(self.total))
        print("TOTAL ----- " ,  self.total)

        

    #######################################################################   
    def CalculateExplicitViscosityContribution(self):
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 0) #explicit contribution by viscosity is defined as fract step = 0
        (self.ExplicitStrategy).InitializeSolutionStep(self.model_part.ProcessInfo);
        (self.ExplicitStrategy).AssembleLoop(self.model_part.ProcessInfo);
        (self.ExplicitStrategy).FinalizeSolutionStep(self.model_part.ProcessInfo);

    #######################################################################   
    def CalculatePressureProjection(self):
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 3) #explicit contribution by viscosity is defined as fract step = 0
        (self.ExplicitStrategy).InitializeSolutionStep(self.model_part.ProcessInfo);
        (self.ExplicitStrategy).AssembleLoop(self.model_part.ProcessInfo);
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 10)
        (self.ExplicitStrategy).FinalizeSolutionStep(self.model_part.ProcessInfo);

    ####################################################################### 
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
