from __future__ import print_function
from KratosMultiphysics import *
from KratosMultiphysics.PFEM2Application import *
#from KratosMultiphysics.ConvectionDiffusionApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
#from KratosMultiphysics.OpenCLApplication import *
#from KratosMultiphysics.GPUSolversApplication import *


from math import sqrt
import time as timer
#dir(pure_diffusion_application)

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(DISTANCE);
    model_part.AddNodalSolutionStepVariable(CORRECTED_DISTANCE);
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
    model_part.AddNodalSolutionStepVariable(SPLIT_ELEMENT)
    model_part.AddNodalSolutionStepVariable(PRESSUREAUX)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)


def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);
        node.AddDof(DISTANCE);


    print("DoFs for the solver added correctly")

class PFEM2Solver:
    def __init__(self,model_part,domain_size):
        self.model_part = model_part
        
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        #definition of the solvers
        gmres_size = 50
        tol = 1e-5
        verbosity = 1
        pDiagPrecond = DiagonalPreconditioner()
        #self.linear_solver = BICGSTABSolver(1e-5, 1000,pDiagPrecond) # SkylineLUFactorizationSolver() 
        self.linear_solver =  AMGCLSolver(AMGCLSmoother.ILU0,AMGCLIterativeSolverType.BICGSTAB,tol,500,verbosity,gmres_size) #SkylineLUFactorizationSolver() 
        self.monolitic_linear_solver =  AMGCLSolver(AMGCLSmoother.ILU0,AMGCLIterativeSolverType.BICGSTAB,tol,500,verbosity,gmres_size)      #BICGSTABSolver(1e-7, 5000) # SkylineLUFactorizationSolver() 
        #(self.monolitic_linear_solver).is_symmetric=True;
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
        pDiagPrecond = DiagonalPreconditioner()
        maximum_number_of_particles= 8*self.domain_size

        self.ExplicitStrategy=PFEM2_Explicit_Strategy(self.model_part,self.domain_size, MoveMeshFlag)

        self.VariableUtils = VariableUtils()
        
        if self.domain_size==2:
                self.moveparticles = MoveParticleUtilityDiffFluidOnly2D(self.model_part,maximum_number_of_particles)
        else:
                self.moveparticles = MoveParticleUtilityDiffFluidOnly3D(self.model_part,maximum_number_of_particles)

        print( "self.domain_size = ", self.domain_size)
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
        
        if self.domain_size==2:
                self.addBC = AddMonolithicFixedVelocityCondition2D(self.model_part)
        else:
                self.addBC = AddMonolithicFixedVelocityCondition3D(self.model_part)

        (self.addBC).AddThem()
        

        ##RICC TOOLS
        #if self.domain_size==2:
        #    self.convection_solver = PureConvectionUtilities2D();
        #else:
        #    self.convection_solver = PureConvectionUtilities3D();
        #self.convection_solver.ConstructSystem(self.model_part,DISTANCE,VELOCITY,FORCE);
        #if self.domain_size==2:
        #   self.distance_utils = SignedDistanceCalculationUtils2D()
        #else:
        #    self.distance_utils = SignedDistanceCalculationUtils3D()
        #self.redistance_step = 0 
        ########

        import strategy_python #implicit solver

        self.monolitic_solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.monolitic_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)

        #to test time
        self.calculateinitialdrag=0.0
        self.streamlineintegration=0.0
        self.accelerateparticles=0.0
        self.implicit_solving=0.0
        self.lagrangiantoeulerian=0.0
        self.reseed=0.0
        self.prereseed=0.0
        self.erasing=0.0
        self.total=0.0
        self.nodaltasks=0.0
        self.thermal=0.0
        self.implicitviscosity=0.0
        
        

                 
    #######################################################################   
    #######################################################################   
    def Solve(self):
                add_gravity=False #in the monolitic solver we do not add the gravity, it is added directly in the implicit monolitic system
                transfer_pressure=False
                viscosity_streamline_integrate=False; #True means explicit integration! False means we will have to solve the fractional velocity system implicetely later, once information has been passed to mesh

                t1 = timer.time()

                #calculating RHS by viscous forces using information of the previous time step: 
                #self.CalculateExplicitViscosityContribution();
                t2 = timer.time()
                self.calculateinitialdrag = self.calculateinitialdrag + t2-t1

                #in order to define a reasonable number of substeps we calculate a mean courant in each element
                (self.moveparticles).CalculateVelOverElemSize();
                t2a = timer.time()

                #streamline integration:
                self.moveparticles.MountBinDiff()
                (self.moveparticles).MoveParticlesDiff(viscosity_streamline_integrate,add_gravity); 
                t3 = timer.time()
                self.streamlineintegration = self.streamlineintegration + t3-t2a
                t3a = timer.time()

                #Reseeding using streamlines in inverse way in case there are too few particles. Not very accurate since particles follow streamlines, no matter if it's water or air.
                #(for 1 fluid should be accurate though
                pre_minimum_number_of_particles=self.domain_size;
                (self.moveparticles).PreReseed(viscosity_streamline_integrate,add_gravity,pre_minimum_number_of_particles);
                t4 = timer.time()
                self.reseed = self.reseed + t4-t3a
                self.prereseed = self.prereseed + t4-t3a


                #RICC THINGS for eulerian transport of distance function
                #self.distance_utils.CalculateDistances(self.model_part,DISTANCE,0.1)
                #(self.VariableUtils).CopyScalarVar(DISTANCE,TEMP_CONV_PROJ,self.model_part.Nodes)

                #transfering data from the particles to the mesh:
                (self.moveparticles).TransferLagrangianToEulerian();

                #MORE RICC THINGS
                #(self.VariableUtils).CopyScalarVar(TEMP_CONV_PROJ,DISTANCE,self.model_part.Nodes)
                #self.convection_solver.CalculateProjection(self.model_part,DISTANCE,NODAL_AREA,VELOCITY,FORCE,TEMP_CONV_PROJ);
                #self.convection_solver.ConvectScalarVar(self.model_part,self.linear_solver,DISTANCE,VELOCITY,FORCE,TEMP_CONV_PROJ,2);

                t5 = timer.time()
                self.lagrangiantoeulerian = self.lagrangiantoeulerian + t5-t4

                #(self.moveparticles).CorrectFreeSurface();
                (self.moveparticles).FlagSplittedElementsAndTheirNodes();
                t6 = timer.time()
                for node in self.model_part.Nodes:
                    distance= node.GetSolutionStepValue(DISTANCE)
                    node.SetSolutionStepValue(CORRECTED_DISTANCE,distance-2.0)
                #TransferLagrtoEul copied info to MESH_VEL, but the solver uses simply velocity, so we copy the variable there.
                (self.VariableUtils).CopyVectorVar(MESH_VELOCITY,VELOCITY,self.model_part.Nodes)
                
                #implicit everything
                full_reset=True;
                self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 20)
                (self.moveparticles).ResetBoundaryConditions(full_reset) 
                (self.monolitic_solver).Solve() #implicit resolution of the pressure system. All the other tasks are explicit
                self.CalculatePressureProjection()
                full_reset=True;
                (self.moveparticles).ResetBoundaryConditions(full_reset) 
                #delta_velocity= Velocity(final) - MeshVelocity(from the particles), so we add to the particles the correction done in the mesh.
                (self.moveparticles).CalculateDeltaVelocity();                
                t11 = timer.time()
                self.implicit_solving = self.implicit_solving + t11-t6
                #transfering the information to the mesh:
                modify_particle_pressure=False
                (self.moveparticles).AccelerateParticlesWithoutMovingUsingDeltaVelocity(modify_particle_pressure,add_gravity);
                t12 = timer.time()
                self.accelerateparticles = self.accelerateparticles + t12-t11

                #reseeding in elements that have few particles to avoid having problems in next iterations:
                post_minimum_number_of_particles=self.domain_size*2;

                self.water_volume = (self.calculatewatervolume).Calculate()
                if (self.water_initial_volume_flag): 
                        print("calculated water volume for the first time")
                        self.water_initial_volume=self.water_volume
                        self.water_initial_volume_flag=False
                print(self.water_volume)
                water_fraction= self.water_volume/(self.water_initial_volume+1e-9)
                self.mass_correction_factor = (1.0 - water_fraction) * 100.0 * 0.5
                if self.mass_correction_factor>0.5:
                        self.mass_correction_factor
                print("mass correction factor: ", self.mass_correction_factor)
                print("current mass loss is : " , (1.0 - water_fraction) * 100.0 , " % ")
                (self.moveparticles).PostReseed(post_minimum_number_of_particles,self.mass_correction_factor);                
                t13 = timer.time()
                self.reseed=self.reseed+ t13-t12
                #self.nodaltasks = self.nodaltasks + t11-t9

                #print ". erasing  = ", t14-t13
                self.total = self.total + t13-t1 
                

                print("----------TIMES----------")
                print("self.calculateinitialdrag  ", self.calculateinitialdrag,"in % = ", 100.0*(self.calculateinitialdrag)/(self.total))
                print("self.streamlineintegration ",  self.streamlineintegration , "in % = ", 100.0*(self.streamlineintegration)/(self.total))
                print("self.accelerateparticles ",  self.accelerateparticles , "in % = ", 100.0*(self.accelerateparticles)/(self.total))
                print("self.implicit_solving ",  self.implicit_solving , "in % = ", 100.0*(self.implicit_solving)/(self.total))
                print("self.nodaltasks ", self.nodaltasks , (self.nodaltasks)/(self.total))
                print("self.lagrangiantoeulerian ",  self.lagrangiantoeulerian , "in % = ", 100.0*(self.lagrangiantoeulerian)/(self.total))
                print("self.reseed ",  self.reseed , "in % = ", 100.0*(self.reseed)/(self.total))
                print("self.prereseed ",  self.prereseed , "in % = ", 100.0*(self.prereseed)/(self.total))
                print("TOTAL ----- ",  self.total)
                print("THIS TIME STEP = ", t13-t1)

                


    #######################################################################   
    def CalculateExplicitViscosityContribution(self):
                self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 0) #explicit contribution by viscosity is defined as fract step = 0
                (self.ExplicitStrategy).InitializeSolutionStep(self.model_part.ProcessInfo);
                (self.ExplicitStrategy).AssembleLoop(self.model_part.ProcessInfo);
                (self.ExplicitStrategy).FinalizeSolutionStep(self.model_part.ProcessInfo);

    def CalculatePressureProjection(self):
                self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 3)
                (self.ExplicitStrategy).InitializeSolutionStep(self.model_part.ProcessInfo);
                (self.ExplicitStrategy).AssembleLoop(self.model_part.ProcessInfo);
                self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 10)
                (self.ExplicitStrategy).FinalizeSolutionStep(self.model_part.ProcessInfo);

    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
