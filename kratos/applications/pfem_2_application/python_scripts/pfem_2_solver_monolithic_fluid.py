from __future__ import print_function
from KratosMultiphysics import *
from KratosMultiphysics.PFEM2Application import *
from KratosMultiphysics.ExternalSolversApplication import *
import KratosMultiphysics as kratoscore
#from KratosMultiphysics.OpenCLApplication import *        #in case you want to use the gpu to solve the system
from math import sqrt
import time as timer

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
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)

    model_part.AddNodalSolutionStepVariable(VISCOSITY_AIR)
    model_part.AddNodalSolutionStepVariable(VISCOSITY_WATER)
    model_part.AddNodalSolutionStepVariable(DENSITY_AIR)
    model_part.AddNodalSolutionStepVariable(DENSITY_WATER)
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
    def __init__(self,model_part,domain_size,maximum_nonlin_iterations):
        self.model_part = model_part
        
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        self.maximum_nonlin_iterations = maximum_nonlin_iterations
        #definition of the solvers
        gmres_size = 50
        tol = 1e-5
        verbosity = 0
        if (maximum_nonlin_iterations==1):
              verbosity = 1 #if it is a linear problem, we can print more information of the single iteration without filling the screen with too much info.
        pDiagPrecond = DiagonalPreconditioner()
        #self.monolithic_linear_solver = BICGSTABSolver(1e-5, 5000,pDiagPrecond) # SkylineLUFactorizationSolver() 
        #self.monolithic_linear_solver =  ViennaCLSolver(tol,500,OpenCLPrecision.Double,OpenCLSolverType.CG,OpenCLPreconditionerType.AMG_DAMPED_JACOBI) #
        #self.monolithic_linear_solver=AMGCLSolver(AMGCLSmoother.ILU0,AMGCLIterativeSolverType.BICGSTAB,tol,1000,verbosity,gmres_size)      #BICGSTABSolver(1e-7, 5000) # SkylineLUFactorizationSolver(	
        '''
        settings = Parameters("""{
                                       "solver_type" : "AMGCL_NS_Solver",
                                       "velocity_block_preconditioner" :
                                            {
                                            "tolerance" : 1e-3,
                                            "preconditioner_type" : "ilu0"
                                        },
                                        "pressure_block_preconditioner" :
                                            {
                                            "tolerance" : 1e-2,
                                            "preconditioner_type" : "ilu0"
                                        },
                                       "tolerance" : 1e-5,
                                       "krylov_type": "bicgstab",
                                       "gmres_krylov_space_dimension": 50,
                                       "coarsening_type": "aggregation",
                                       "max_iteration": 50,
                                       "verbosity" : 2,
                                       "scaling": false,
                                       "coarse_enough" : 5000
                                   } """)
        linear_solver = AMGCL_NS_Solver(settings)
        '''
        
        #construct the linear solvers
        import linear_solver_factory
        #self.monolithic_linear_solver =  linear_solver
        self.monolithic_linear_solver = AMGCLSolver(AMGCLSmoother.ILU0,AMGCLIterativeSolverType.BICGSTAB,tol,1000,verbosity,gmres_size)
        self.conv_criteria = DisplacementCriteria(1e-3,1e-3)  #tolerance for the solver 
        self.conv_criteria.SetEchoLevel(0)
        
        self.domain_size = domain_size
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
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
        #raw_input()
        print("self.domain_size = ", self.domain_size)
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
              self.addBC = AddFixedPressureCondition2D(self.model_part)
        else:
              self.addBC = AddFixedPressureCondition3D(self.model_part)

        (self.addBC).AddThem()
        
        if self.domain_size==2:
             self.distance_utils = SignedDistanceCalculationUtils2D()
        else:
             self.distance_utils = SignedDistanceCalculationUtils3D()
        self.redistance_step = 0

        import strategy_python #implicit solver
        self.monolithic_solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.monolithic_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
        self.monolithic_solver.SetMaximumIterations(self.maximum_nonlin_iterations)

        self.streamlineintegration=0.0
        self.accelerateparticles=0.0
        self.implicitsystem=0.0
        self.lagrangiantoeulerian=0.0
        self.reseed=0.0
        self.prereseed=0.0
        self.erasing=0.0
        self.total=0.0

        self.model_part.ProcessInfo.SetValue(VOLUME_CORRECTION, 0.0)

        self.print_times=False

        visc_water=self.model_part.ProcessInfo.GetValue(VISCOSITY_WATER)
        visc_air=self.model_part.ProcessInfo.GetValue(VISCOSITY_AIR)
        dens_water=self.model_part.ProcessInfo.GetValue(DENSITY_WATER)
        dens_air=self.model_part.ProcessInfo.GetValue(DENSITY_AIR)

        for node in self.model_part.Nodes:
              node.SetSolutionStepValue(VISCOSITY_WATER,visc_water)
              node.SetSolutionStepValue(VISCOSITY_AIR,visc_air)
              node.SetSolutionStepValue(DENSITY_WATER,dens_water)
              node.SetSolutionStepValue(DENSITY_AIR,dens_air)
                 
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

        (self.VariableUtils).CopyVectorVar(PROJECTED_VELOCITY,VELOCITY,self.model_part.Nodes)
        full_reset=True;
        (self.moveparticles).ResetBoundaryConditions(full_reset) 
        (self.moveparticles).CopyVectorVarToPreviousTimeStep(VELOCITY,self.model_part.Nodes)

        t6 = timer.time()
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 20)
        (self.monolithic_solver).Solve() #implicit resolution of the system.

        t7 = timer.time()
        self.implicitsystem = self.implicitsystem + t7-t6
        self.CalculatePressureProjection()
        
        #delta_velocity= Velocity(final) - ProjectedVelocity(from the particles), so we add to the particles the correction done in the mesh.
        (self.moveparticles).CalculateDeltaVelocity();        
        t11 = timer.time()
        (self.moveparticles).AccelerateParticlesWithoutMovingUsingDeltaVelocity();
        t12 = timer.time()
        self.accelerateparticles = self.accelerateparticles + t12-t11

        #reseeding in elements that have few particles to avoid having problems in next iterations:
        post_minimum_number_of_particles=self.domain_size*2;
        (self.moveparticles).PostReseed(post_minimum_number_of_particles,self.mass_correction_factor);        
        t13 = timer.time()
        self.reseed=self.reseed+ t13-t12

        #calculating water volume to correct mass
        self.water_volume = (self.calculatewatervolume).Calculate()
        if (self.water_initial_volume==0.0): 
                self.water_initial_volume=self.water_volume
        water_fraction= self.water_volume/(self.water_initial_volume)
        self.mass_correction_factor = (1.0 - water_fraction) * 100.0 * 0.1
        print("current mass loss is : " , (1.0 - water_fraction) * 100.0 , " % ")
        #print("water fraction ", water_fraction)
        #print("water volume ", self.water_volume)
        print("mass correction factor: ", self.mass_correction_factor)
        self.model_part.ProcessInfo.SetValue(VOLUME_CORRECTION, self.mass_correction_factor)



        #self.nodaltasks = self.nodaltasks + t11-t9

        #print ". erasing  = ", t14-t13
        self.total = self.total + t13-t1 
        
        if self.print_times==True:
                print( "----------TIMES----------" )
                print( "self.streamlineintegration " ,  self.streamlineintegration , "in % = ", 100.0*(self.streamlineintegration)/(self.total) )
                print( "self.accelerateparticles " ,  self.accelerateparticles , "in % = ", 100.0*(self.accelerateparticles)/(self.total) )
                print( "self.implicitsystem " ,  self.implicitsystem , "in % = ", 100.0*(self.implicitsystem)/(self.total) )
                print( "self.lagrangiantoeulerian " ,  self.lagrangiantoeulerian , "in % = ", 100.0*(self.lagrangiantoeulerian)/(self.total) )
                print( "self.reseed " ,  self.reseed , "in % = ", 100.0*(self.reseed)/(self.total) )
                print( "self.prereseed " ,  self.prereseed , "in % = ", 100.0*(self.prereseed)/(self.total) )
                print( "TOTAL ----- " ,  self.total )
                print( "this time step" ,  t13-t1  )
                print( "current time in simulation: ", self.model_part.ProcessInfo.GetValue(TIME) , "s" )
        
    
    def RotateParticlesAndDomainVelocities(self,angles):
        (self.moveparticles).RotateParticlesAndDomainVelocities(angles)
    #######################################################################   
    def CalculatePressureProjection(self):
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 10)
        (self.ExplicitStrategy).InitializeSolutionStep();
        (self.ExplicitStrategy).AssembleLoop();
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 10)
        (self.ExplicitStrategy).FinalizeSolutionStep();


    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)

    def PrintInfo(self,print_times):
        self.print_times=print_times

    
