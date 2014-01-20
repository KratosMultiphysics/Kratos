from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.PFEM2Application import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
from math import sqrt
import time as timer
# dir(pure_diffusion_application)


def AddVariables(model_part, linea_model_part):
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(FRACT_VEL)
    model_part.AddNodalSolutionStepVariable(G_VALUE)
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(YP);
    # model_part.AddNodalSolutionStepVariable(ERASE_FLAG);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    # model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    model_part.AddNodalSolutionStepVariable(RHS);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(NORMAL);
    # model_part.AddNodalSolutionStepVariable(DENSITY);
    # model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ);
    model_part.AddNodalSolutionStepVariable(PREVIOUS_ITERATION_PRESSURE);
    # model_part.AddNodalSolutionStepVariable(FIRST_ITERATION_PRESSURE);
    model_part.AddNodalSolutionStepVariable(DELTA_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ_NO_RO)
    model_part.AddNodalSolutionStepVariable(MEAN_SIZE)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(MASS)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(TAU)


    # linea_model_part.AddNodalSolutionStepVariable(ERASE_FLAG);
    # linea_model_part.AddNodalSolutionStepVariable(DISTANCE);
    # linea_model_part.AddNodalSolutionStepVariable(VELOCITY);
    # linea_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    # linea_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    # linea_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    # linea_model_part.AddNodalSolutionStepVariable(NEIGHBOUR_ELEMENTS);
    # linea_model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    # linea_model_part.AddNodalSolutionStepVariable(PRESSURE);
def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);
        node.AddDof(FRACT_VEL_X);
        node.AddDof(FRACT_VEL_Y);
        node.AddDof(FRACT_VEL_Z);
        node.AddDof(DISTANCE);

    print("variables for the Poisson solver added correctly")


class PFEM2Solver:
    #
    # def __init__(self,model_part,linea_model_part,domain_size):

    def __init__(self, model_part, domain_size):
        self.model_part = model_part

        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        # definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        self.pressure_linear_solver = BICGSTABSolver(1e-6, 5000, pDiagPrecond)
        # self.pressure_linear_solver =  SkylineLUFactorizationSolver()  #we set the type of solver that we want
        self.conv_criteria = DisplacementCriteria(1e-9, 1e-15)  # tolerance for the solver

        self.domain_size = domain_size
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)
        (self.neighbour_search).Execute()
        self.neighbour_elements_search = FindElementalNeighboursProcess(model_part, domain_size, number_of_avg_elems)
        (self.neighbour_elements_search).Execute()
        # calculate normals
        self.normal_tools = BodyNormalCalculationUtils()

    #
    #
    def Initialize(self):
        # creating the solution strategy
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = True
        MoveMeshFlag = False

        # build the edge data structure
        # if self.domain_size==2:
        #	self.matrix_container = MatrixContainer2D()
        # else:
        #	self.matrix_container = MatrixContainer3D()
        maximum_number_of_particles = 10 * self.domain_size

        self.ExplicitStrategy = PFEM2_Explicit_Strategy(self.model_part, self.domain_size, MoveMeshFlag)

        self.VariableUtils = VariableUtils()

        if self.domain_size == 2:
            self.moveparticles = MoveParticleUtilityDiff2D(self.model_part, maximum_number_of_particles)
        else:
            self.moveparticles = MoveParticleUtilityDiff3D(self.model_part, maximum_number_of_particles)

        self.calculatewatervolume = CalculateWaterFraction(self.model_part)

        self.moveparticles.MountBinDiff()
        # self.matrix_container.ConstructCSRVector(self.model_part)
        # self.matrix_container.BuildCSRData(self.model_part)
        self.water_volume = 0.0  # we initialize it at zero
        self.water_initial_volume = 0.0  # we initialize it at zero
        self.mass_correction_factor = 0.0
        # creating the solution strategy
        # CalculateReactionFlag = False
        # ReformDofSetAtEachStep = False
        # MoveMeshFlag = True
        # import strategy_python
        # self.solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.poisson_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)

        # constructing the solver
        # if self.domain_size==2:
        #   self.solver = DifussionSolver2D(self.matrix_container,self.model_part)
        # else:
        #   self.solver = DifussionSolver3D(self.matrix_container,self.model_part)
        # print "ln90"
        # self.solver.Initialize(h,radius)
        # print "ln91"

        self.normal_tools.CalculateBodyNormals(self.model_part, self.domain_size);
        condition_number = 1

        if self.domain_size == 2:
            self.addBC = AddFixedVelocityCondition2D(self.model_part)
        else:
            self.addBC = AddFixedVelocityCondition3D(self.model_part)

        (self.addBC).AddThem()

        # RICC TOOLS
        # self.convection_solver = PureConvectionUtilities2D();
        # self.convection_solver.ConstructSystem(self.model_part,DISTANCE,VELOCITY,MESH_VELOCITY);
        # self.distance_utils = SignedDistanceCalculationUtils2D()
        # self.redistance_step = 0

        import strategy_python  # implicit solver
        self.solver = strategy_python.SolvingStrategyPython(self.model_part, self.time_scheme, self.pressure_linear_solver, self.conv_criteria, CalculateReactionFlag, ReformDofSetAtEachStep, MoveMeshFlag)
        # pDiagPrecond = DiagonalPreconditioner()
        # self.linear_solver =  BICGSTABSolver(1e-9, 5000,pDiagPrecond)
        #(self.moveparticles).ReplaceParticlesVelocityAndDistance();
        self.calculateinitialdrag = 0.0
        self.streamlineintegration = 0.0
        self.accelerateparticles = 0.0
        self.navierstokes = 0.0
        self.lagrangiantoeulerian = 0.0
        self.reseed = 0.0
        self.prereseed = 0.0
        self.erasing = 0.0
        self.total = 0.0
        self.nodaltasks = 0.0

    #
    def Solve(self):
        viscosity_streamline_integrate = False;  # True means explicit integration! False means we will have to solve the fractional velocity system implicetely later, once information has been passed to mesh
        # pressure_gradient_integrate=False; #we use the pressure of the previous step. not reccomended for 2 fluid flows. REMOVED FROM THE CODE
        calculate_viscous_contribution_after_everything = True;
        # implicit_velocity_correction=False;	# as expected this can not be active the implicit correction is done:
        # if (calculate_viscous_contribution_after_everything):
        #	implicit_velocity_correction=False;

        t1 = timer.time()

        # calculating RHS by viscous forces using information of the previous time step:
        self.CalculateExplicitViscosityContribution();
        t2 = timer.time()
        self.calculateinitialdrag = self.calculateinitialdrag + t2 - t1

        # in order to define a reasonable number of substeps we calculate a mean courant in each element
        (self.moveparticles).CalculateVelOverElemSize();
        t2a = timer.time()

        # streamline integration:
        (self.moveparticles).MoveParticlesDiff(viscosity_streamline_integrate);
        t3 = timer.time()
        self.streamlineintegration = self.streamlineintegration + t3 - t2a
        t3a = timer.time()

        # Reseeding using streamlines in inverse way in case there are too few particles. Not very accurate since particles follow streamlines, no matter if it's water or air.
        #(for 1 fluid should be accurate though
        pre_minimum_number_of_particles = 3;
        (self.moveparticles).PreReseed(viscosity_streamline_integrate, pre_minimum_number_of_particles);
        t4 = timer.time()
        self.reseed = self.reseed + t4 - t3a
        self.prereseed = self.prereseed + t4 - t3a

        # transfering data from the particles to the mesh:
        (self.moveparticles).TransferLagrangianToEulerian();
        t5 = timer.time()
        self.lagrangiantoeulerian = self.lagrangiantoeulerian + t5 - t4

        # Flagging splitted elements (the ones who need special integration
        (self.moveparticles).FlagSplittedElementsAndTheirNodes();

        # if (implicit_velocity_correction==False): #already solved in the streamlines, no need to do it now.
        # (self.moveparticles).MoveDataFromMESH_VELOCITYtoFRACT_VEL(); #insted of first step,. when moving through streamlines the first step is already performed
        # else: #we did not do it. doing it explicitely
        #	(self.moveparticles).CalculateFirstStep();
                # self.model_part.ProcessInfo.SetValue(NL_ITERATION_NUMBER, 0) #1= semiimplicit, =0 fully implciit
                #(self.solver).Solve()
        t6 = timer.time()
        # pressure iterations
        number_of_pressure_iterations = 1
        for i in range(0, number_of_pressure_iterations):
            print("pressure iteration number: ", i + 1)
            self.CalculatePressureIteration(i + 1)

        # resetting boundary conditions
        (self.moveparticles).ResetBoundaryConditions()

        # implicit viscosity
        if (calculate_viscous_contribution_after_everything):
            self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 4)
            print("doing viscous implict step after everything!!!!!!!!!!!!!!!")
            (self.solver).Solve()

        number_of_pressure_iterations = 2
        for i in range(1, number_of_pressure_iterations):
            print("pressure iteration number: ", i + 1)
            self.CalculatePressureIteration(i + 1)

        (self.moveparticles).ResetBoundaryConditions()
        # implicit viscosity
        # if (calculate_viscous_contribution_after_everything):
        #	self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 4)
        #	print "doing viscous implict step after everything!!!!!!!!!!!!!!!"
        #	(self.solver).Solve()

        # delta_velocity= Velocity(final) - MeshVelocity(from the particles), so we add to the particles the correction done in the mesh.
        (self.moveparticles).CalculateDeltaVelocity();
        t11 = timer.time()
        self.navierstokes = self.navierstokes + t11 - t6
        # transfering the information to the mesh:
        (self.moveparticles).AccelerateParticlesWithoutMovingUsingDeltaVelocity();
        t12 = timer.time()
        self.accelerateparticles = self.accelerateparticles + t12 - t11

        # reseeding in elements that have few particles to avoid having problems in next iterations:
        post_minimum_number_of_particles = self.domain_size * 3;
        self.water_volume = (self.calculatewatervolume).Calculate()
        if (self.water_initial_volume == 0.0):
            self.water_initial_volume = self.water_volume
        water_fraction = self.water_volume / (self.water_initial_volume + 1e-9)
        self.mass_correction_factor = (1.0 - water_fraction) * 100.0 * 0.001
        print("mass correction factor: ", self.mass_correction_factor)
        (self.moveparticles).PostReseed(post_minimum_number_of_particles, self.mass_correction_factor);
        t13 = timer.time()
        self.reseed = self.reseed + t13 - t12
        # self.nodaltasks = self.nodaltasks + t11-t9

        # print ". erasing  = ", t14-t13
        self.total = self.total + t13 - t1
        '''
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
		'''

    #
    def CalculatePressureIteration(self, iteration_number):
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 2)
        self.model_part.ProcessInfo.SetValue(NL_ITERATION_NUMBER, iteration_number)
        if (iteration_number == 1):
            (self.VariableUtils).CopyVectorVar(MESH_VELOCITY, VELOCITY, self.model_part.Nodes)
        (self.solver).Solve()  # implicit resolution of the pressure system. All the other tasks are explicit

        # setting the Fract step number to 3 to calculate the pressure gradient effect on the velocity
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 3)
        (self.ExplicitStrategy).InitializeSolutionStep(self.model_part.ProcessInfo);
        (self.ExplicitStrategy).AssembleLoop(self.model_part.ProcessInfo);
        (self.ExplicitStrategy).FinalizeSolutionStep(self.model_part.ProcessInfo);

    #
    def CalculateExplicitViscosityContribution(self):
        self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 0)  # explicit contribution by viscosity is defined as fract step = 0
        (self.ExplicitStrategy).InitializeSolutionStep(self.model_part.ProcessInfo);
        (self.ExplicitStrategy).AssembleLoop(self.model_part.ProcessInfo);
        (self.ExplicitStrategy).FinalizeSolutionStep(self.model_part.ProcessInfo);

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)
