from __future__ import print_function
from KratosMultiphysics import *
from KratosMultiphysics.PFEM2Application import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.OpenCLApplication import *

from math import sqrt
import time as timer
#dir(pure_diffusion_application)

def AddVariables(model_part,topographic_model_part):
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(SOLID_PRESSURE);
    model_part.AddNodalSolutionStepVariable(SOLID_YP);
    model_part.AddNodalSolutionStepVariable(DISTANCE);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(FRACT_VEL);
    model_part.AddNodalSolutionStepVariable(G_VALUE);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(YP);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(NORMAL);	
    model_part.AddNodalSolutionStepVariable(DELTA_VELOCITY) 
    model_part.AddNodalSolutionStepVariable(NODAL_AREA) 
    model_part.AddNodalSolutionStepVariable(NODAL_MASS) 
    model_part.AddNodalSolutionStepVariable(TEMPERATURE)
    model_part.AddNodalSolutionStepVariable(TEMPERATURE_OLD_IT)
    model_part.AddNodalSolutionStepVariable(PRESSUREAUX)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(RHS);
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    
    topographic_model_part.AddNodalSolutionStepVariable(PRESSURE);
    topographic_model_part.AddNodalSolutionStepVariable(VELOCITY);
    topographic_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

def AddDofs(model_part,topographic_model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);
        node.AddDof(DISTANCE);
        #node.AddDof(TEMPERATURE);
    for node in topographic_model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);
        #node.AddDof(TEMPERATURE);


    print("DoFs for the solver added correctly")

class PFEM2Solver:
    #######################################################################
    #def __init__(self,model_part,linea_model_part,domain_size):
    def __init__(self,model_part,topographic_model_part,domain_size):
        self.model_part = model_part
        self.topographic_model_part = topographic_model_part
        
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        #definition of the solvers
        gmres_size = 50
        tol = 1e-5
        verbosity = 1
        #self.monolitic_linear_solver = SkylineLUFactorizationSolver() 
        #self.monolitic_linear_solver = BICGSTABSolver(1e-5, 1000,pDiagPrecond) # SkylineLUFactorizationSolver() 
        #self.monolitic_linear_solver= ViennaCLSolver(tol,10000,OpenCLPrecision.Double,OpenCLSolverType.CG,OpenCLPreconditionerType.AMG_DAMPED_JACOBI) 
        self.topographic_monolitic_linear_solver= AMGCLSolver(AMGCLSmoother.ILU0,AMGCLIterativeSolverType.CG,tol,1000,verbosity,gmres_size) 
        self.monolitic_linear_solver =  AMGCLSolver(AMGCLSmoother.ILU0,AMGCLIterativeSolverType.CG,tol,1000,verbosity,gmres_size)      #BICGSTABSolver(1e-7, 5000) # SkylineLUFactorizationSolver() 
        self.conv_criteria = DisplacementCriteria(1e-2,1e-6)  #tolerance for the solver 
        
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
    def Initialize(self,initial_offset):
        #creating the solution strategy
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False	
        MoveMeshFlag = False
        maximum_number_of_particles= 10*self.domain_size
        
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
        
        #to search particles. 
        self.moveparticles.MountBinDiff()
        self.water_volume=0.0  #we initialize it at zero
        self.water_initial_volume=0.0 #we initialize it at zero
        self.water_initial_volume_flag=True #we initialize it at zero
        self.mass_correction_factor=0.0


        
        #body normals, useful for boundary (inlet) conditions
        self.normal_tools.CalculateBodyNormals(self.model_part,self.domain_size);  
        condition_number=1
        	
        
        #solvers:
	import strategy_python_nonlinear #implicit solver
	import strategy_python #implicit solver
        self.monolitic_solver = strategy_python_nonlinear.SolvingStrategyPython(self.model_part,self.time_scheme,self.monolitic_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
        self.topographic_monolitic_solver = strategy_python.SolvingStrategyPython(self.topographic_model_part,self.time_scheme,self.topographic_monolitic_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
        #(self.moveparticles).ReplaceParticlesVelocityAndDistance();


        #to print the total time spent on the simulation
        self.implicit_tasks=0.0
        self.total=0.0


        #finally we solve the topographic system only oonce:
	nsteps=1
	for step in range(0,nsteps):
		(self.topographic_monolitic_solver).Solve()
		(self.moveparticles).CalculateElementalMeanStress(self.topographic_model_part)
        
        overwrite_particle_data=True
        #this is the position where we place the calculation domain
        (self.moveparticles).InitializeTransferTool(self.topographic_model_part, initial_offset,overwrite_particle_data)
      
                 
    #######################################################################   
    def Solve(self,fluid_free_bounding_box_lower_corner,fluid_free_bounding_box_upper_corner):
                t1 = timer.time()

                #some misc flags
                add_gravity=False #in the monolitic solver we do not add the gravity, it is added directly in the implicit monolitic system
		viscosity_streamline_integrate=False; #True means explicit integration! False means we will have to solve the fractional velocity system implicetely later, once information has been passed to mesh

                #in order to define a reasonable number of substeps we calculate a mean courant in each element
                (self.moveparticles).CalculateVelOverElemSize();
                t2 = timer.time()

                #streamline integration:
                (self.moveparticles).MoveParticlesDiff(viscosity_streamline_integrate,add_gravity);	
                t3 = timer.time()
                
                #Reseeding using streamlines in inverse way in case there are too few particles. Not very accurate since particles follow streamlines, no matter if it's water or air.
                #we want at least 2 (or 3 in 3d) per element)
                pre_minimum_number_of_particles=self.domain_size;
                (self.moveparticles).PreReseedUsingTopographicDomain(pre_minimum_number_of_particles)
                t4 = timer.time()
		
		#transfering data from the particles to the mesh:
		transfer_pressure=True
                (self.moveparticles).TransferLagrangianToEulerian(transfer_pressure);
		
                #storing variables in the right places
                #for node in self.model_part.Nodes:
                #        node.SetSolutionStepValue(VELOCITY,node.GetSolutionStepValue(MESH_VELOCITY))
                #        node.SetSolutionStepValue(MESH_VELOCITY,node.GetSolutionStepValue(VELOCITY))
                (self.VariableUtils).CopyVectorVar(MESH_VELOCITY,VELOCITY,self.model_part.Nodes)	
                #(self.VariableUtils).CopyVectorVar(VELOCITY,MESH_VELOCITY,self.model_part.Nodes)	
                (self.moveparticles).ResetBoundaryConditions(True) 
                #(self.VariableUtils).CopyScalarVar(SOLID_PRESSURE,PRESSURE,self.model_part.Nodes)	
                #(self.VariableUtils).CopyScalarVar(SOLID_PRESSURE,PRESSUREAUX,self.model_part.Nodes)	
                (self.VariableUtils).CopyScalarVar(PRESSURE,SOLID_PRESSURE,self.model_part.Nodes)	
                #(self.VariableUtils).CopyScalarVar(CORRECTED_DISTANCE,DISTANCE,self.model_part.Nodes)		
                (self.moveparticles).CopyVectorVarToPreviousTimeStep(VELOCITY,self.model_part.Nodes)
                (self.moveparticles).CopyScalarVarToPreviousTimeStep(PRESSURE,self.model_part.Nodes)
                (self.VariableUtils).CopyScalarVar(PRESSURE,PRESSUREAUX,self.model_part.Nodes)	
                
                t5 = timer.time()
                
                full_reset=False;
                #(self.moveparticles).ResetBoundaryConditions(full_reset) 
                
                #implicit solver!
            	self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 20)
                #(self.monolitic_solver).Solve() #implicit resolution of the pressure system. All the other tasks are explicit
                
                max_iter=20
                (self.monolitic_solver).Solve(self.model_part,self.moveparticles,(self.VariableUtils),max_iter)
 	

                t11 = timer.time()
                self.implicit_tasks = self.implicit_tasks + t11-t5


                full_reset=True;
                (self.moveparticles).ResetBoundaryConditions(full_reset) 


                #transfering the information to the mesh:
                modify_particle_pressure=True
                #delta_velocity= Velocity(final) - MeshVelocity(from the particles), so we add to the particles the correction done in the mesh.
                (self.moveparticles).CalculateDeltaVelocity();	
                (self.moveparticles).AccelerateParticlesWithoutMovingUsingDeltaVelocity(add_gravity);
                t12 = timer.time()
                
                #reseeding in elements that have few particles to avoid having problems in next iterations:
                post_minimum_number_of_particles=self.domain_size*2;
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
                #(self.moveparticles).PostReseed(post_minimum_number_of_particles,self.mass_correction_factor);		
                #(self.moveparticles).PostReseedOnlyInBoundingBox(post_minimum_number_of_particles,self.mass_correction_factor,bounding_box_lower_corner,bounding_box_upper_corner);			
                max_nodal_distance=0.7 #we have to move the domain to include the nodes that have nodal_distance>0.0 and <max_nodal_distance- otherwise there is no movement	
                (self.moveparticles).ComputeCalculationDomainDisplacement(fluid_free_bounding_box_lower_corner,fluid_free_bounding_box_upper_corner,max_nodal_distance)
                
                t13 = timer.time()
                self.total = self.total + t13-t1 
                
                print ("----------TIMES----------")
                print ("self.implicit_tasks  " , self.implicit_tasks , "in % = ",  100.0*(self.implicit_tasks)/(self.total))
                print ("TOTAL ----- " ,  self.total)
		

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
