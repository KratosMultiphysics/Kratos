from KratosMultiphysics import *
from KratosMultiphysics.PFEM2Application import *
from KratosMultiphysics.ExternalSolversApplication import *
#from KratosMultiphysics.OpenCLApplication import *        #and now our application. note that we can import as many as we need to
from math import sqrt
import time as timer
#dir(pure_diffusion_application)

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(DISTANCE);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(YP);
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    model_part.AddNodalSolutionStepVariable(RHS);
    model_part.AddNodalSolutionStepVariable(PROJECTED_VELOCITY);
    model_part.AddNodalSolutionStepVariable(NORMAL);        
    model_part.AddNodalSolutionStepVariable(PREVIOUS_ITERATION_PRESSURE);
    model_part.AddNodalSolutionStepVariable(DELTA_VELOCITY) 
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ_NO_RO) 
    model_part.AddNodalSolutionStepVariable(MEAN_SIZE) 
    model_part.AddNodalSolutionStepVariable(NODAL_AREA) 
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)

def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);
        node.AddDof(DISTANCE);

class PFEM2Solver:
    #######################################################################
    #def __init__(self,model_part,linea_model_part,domain_size):
    def __init__(self,model_part,domain_size):
        self.model_part = model_part
        
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        #definition of the solvers
        gmres_size = 50
        tol = 1e-5
        verbosity = 1
        pDiagPrecond = DiagonalPreconditioner()
        #self.pressure_linear_solver = BICGSTABSolver(1e-5, 5000,pDiagPrecond) # SkylineLUFactorizationSolver() 
        #self.pressure_linear_solver =  ViennaCLSolver(tol,500,OpenCLPrecision.Double,OpenCLSolverType.CG,OpenCLPreconditionerType.NoPreconditioner) #
        self.pressure_linear_solver=AMGCLSolver(AMGCLSmoother.DAMPED_JACOBI,AMGCLIterativeSolverType.CG,tol,200,verbosity,gmres_size)      #BICGSTABSolver(1e-7, 5000) # SkylineLUFactorizationSolver() 
        self.conv_criteria = DisplacementCriteria(1e-9,1e-15)  #tolerance for the solver 
        
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
        
        #build the edge data structure
        #if self.domain_size==2:
        #        self.matrix_container = MatrixContainer2D()
        #else: 
        #        self.matrix_container = MatrixContainer3D()
        #maximum_number_of_particles= 10*self.domain_size
        maximum_number_of_particles= 8*self.domain_size

        self.ExplicitStrategy=PFEM2_Explicit_Strategy(self.model_part,self.domain_size, MoveMeshFlag)

        self.VariableUtils = VariableUtils()
        
        if self.domain_size==2:
             self.moveparticles = MoveParticleUtilityPFEM22D(self.model_part,maximum_number_of_particles)
        else:
             self.moveparticles = MoveParticleUtilityPFEM23D(self.model_part,maximum_number_of_particles)

        print( "self.domain_size = ", self.domain_size)
        if self.domain_size==2:
             self.calculatewatervolume = CalculateWaterFraction2D(self.model_part)
        else:        
             self.calculatewatervolume = CalculateWaterFraction3D(self.model_part)

        self.moveparticles.MountBin()
        self.water_volume=0.0  #we initialize it at zero
        self.water_initial_volume=0.0 #we initialize it at zero
        self.mass_correction_factor=0.0

        self.normal_tools.CalculateBodyNormals(self.model_part,self.domain_size);  
        condition_number=1
        
        if self.domain_size==2:
              self.addBC = AddFixedVelocityCondition2D(self.model_part)
        else:
              self.addBC = AddFixedVelocityCondition3D(self.model_part)

        (self.addBC).AddThem()
        
        if self.domain_size==2:
             self.distance_utils = SignedDistanceCalculationUtils2D()
        else:
             self.distance_utils = SignedDistanceCalculationUtils3D()
        self.redistance_step = 0

        import strategy_python #implicit solver
        self.pressure_solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.pressure_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
        self.streamlineintegration=0.0
        self.accelerateparticles=0.0
        self.implicitpressure=0.0
        self.lagrangiantoeulerian=0.0
        self.reseed=0.0
        self.prereseed=0.0
        self.erasing=0.0
        self.explicitviscous = 0.0
        self.total=0.0


        
      
                 
    #######################################################################   
    def Solve(self):
        t1 = timer.time()
        #in order to define a reasonable number of substeps we calculate a mean courant in each element
        (self.moveparticles).CalculateVelOverElemSize();
        t2= timer.time()

        #streamline integration:
        discriminate_streamlines=True
        (self.moveparticles).MoveParticles(discriminate_streamlines);         
        t3 = timer.time()
        self.streamlineintegration = self.streamlineintegration + t3-t2

        #Reseeding using streamlines in inverse way in case there are too few particles. Not very accurate since particles follow streamlines, no matter if it's water or air.
        #(for 1 fluid should be accurate though
        pre_minimum_number_of_particles=self.domain_size;
        (self.moveparticles).PreReseed(pre_minimum_number_of_particles);
        t4 = timer.time()
        self.reseed = self.reseed + t4-t3
        self.prereseed = self.prereseed + t4-t3

        #transfering data from the particles to the mesh:
        (self.moveparticles).TransferLagrangianToEulerian();
        t5 = timer.time()
        self.lagrangiantoeulerian = self.lagrangiantoeulerian + t5-t4

        print( "finished lagrangian to eulerian")
        (self.VariableUtils).CopyVectorVar(PROJECTED_VELOCITY,VELOCITY,self.model_part.Nodes)

        t7=timer.time()
        #explicit viscosity. automatically limited inside to ensure stability. artifically lowered to keep Fo<0.5
        self.CalculateExplicitViscosityContribution()
        t8=timer.time()
        self.explicitviscous = self.explicitviscous + t8 - t7
        #print("finished explicit viscosity calculation")
        number_of_pressure_iterations=3
        for i in range(0,number_of_pressure_iterations):
                self.CalculatePressureIteration(i+1)
        t9=timer.time()
        self.implicitpressure = self.implicitpressure + t9 - t8
        #setting the appropiate BC, with the pressure eq the boundaries are only weakly impermeable
        full_reset=True;
        (self.moveparticles).ResetBoundaryConditions(full_reset) 
        
        #delta_velocity= Velocity(final) - ProjectedVelocity(from the particles), so we add to the particles the correction done in the mesh.
        (self.moveparticles).CalculateDeltaVelocity();        
        t11 = timer.time()
        (self.moveparticles).AccelerateParticlesWithoutMovingUsingDeltaVelocity();
        t12 = timer.time()
        self.accelerateparticles = self.accelerateparticles + t12-t11

        #reseeding in elements that have few particles to avoid having problems in next iterations:
        post_minimum_number_of_particles=self.domain_size*3;
        (self.moveparticles).PostReseed(post_minimum_number_of_particles,self.mass_correction_factor);        
        t13 = timer.time()
        self.reseed=self.reseed+ t13-t12

        #calculating water volume to correct mass
        self.water_volume = (self.calculatewatervolume).Calculate()
        if (self.water_initial_volume==0.0): 
                self.water_initial_volume=self.water_volume
        water_fraction= self.water_volume/(self.water_initial_volume+1e-9)
        self.mass_correction_factor = (1.0 - water_fraction) * 100.0 * 0.0
        print( "current mass loss is : " , (1.0 - water_fraction) * 100.0 , " % ")
        print( "water fraction ", water_fraction)
        print( "water volume ", self.water_volume)
        print( "mass correction factor: ", self.mass_correction_factor)

        #self.nodaltasks = self.nodaltasks + t11-t9

        #print ". erasing  = ", t14-t13
        self.total = self.total + t13-t1 
        
        if True:
                print( "----------TIMES----------")
                print( "self.streamlineintegration " ,  self.streamlineintegration , "in % = ", 100.0*(self.streamlineintegration)/(self.total))
                print( "self.accelerateparticles " ,  self.accelerateparticles , "in % = ", 100.0*(self.accelerateparticles)/(self.total))
                print( "self.explicitviscous " ,  self.explicitviscous , "in % = ", 100.0*(self.explicitviscous)/(self.total))
                print( "self.implicitpressure " ,  self.implicitpressure , "in % = ", 100.0*(self.implicitpressure)/(self.total))
                print( "self.lagrangiantoeulerian " ,  self.lagrangiantoeulerian , "in % = ", 100.0*(self.lagrangiantoeulerian)/(self.total))
                print( "self.reseed " ,  self.reseed , "in % = ", 100.0*(self.reseed)/(self.total))
                print( "self.prereseed " ,  self.prereseed , "in % = ", 100.0*(self.prereseed)/(self.total))
                print( "TOTAL ----- " ,  self.total)
                print( "this time step" ,  t13-t1 )
                print( "current time in simulation: ", self.model_part.ProcessInfo.GetValue(TIME) , "s")
        

    #######################################################################   
    def CalculatePressureIteration(self, iteration_number):
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 2)
        self.model_part.ProcessInfo.SetValue(NL_ITERATION_NUMBER, iteration_number)
        
        (self.pressure_solver).Solve() #implicit resolution of the pressure system. All the other tasks are explicit
        #if iteration_number==1:
        #        for node in self.model_part.Nodes:
        #        press= node.GetSolutionStepValue(PRESSURE)+node.GetSolutionStepValue(PRESSUREAUX)
        #        node.SetSolutionStepValue(PRESSURE,press)
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

    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
