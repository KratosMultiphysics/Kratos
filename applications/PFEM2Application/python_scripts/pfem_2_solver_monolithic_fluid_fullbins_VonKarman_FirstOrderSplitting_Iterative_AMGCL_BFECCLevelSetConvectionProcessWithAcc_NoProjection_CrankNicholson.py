from __future__ import print_function
from KratosMultiphysics import *
from KratosMultiphysics.PFEM2Application import *
#from KratosMultiphysics.LinearSolversApplication import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
import KratosMultiphysics
#from KratosMultiphysics.OpenCLApplication import *        #in case you want to use the gpu to solve the system
from math import sqrt
from math import fabs
from math import sin
from math import cos
import time as timer
import math


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(VELOCITY_AUX)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(ACCELERATION_AUX)
    model_part.AddNodalSolutionStepVariable(DELTA_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(YP)
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ)
    model_part.AddNodalSolutionStepVariable(RHS)
    model_part.AddNodalSolutionStepVariable(PROJECTED_VELOCITY)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(PREVIOUS_ITERATION_PRESSURE)
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
    
    ########
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(DISTANCE_GRADIENT)
    ########


    ########
    #model_part.AddNodalSolutionStepVariable(YP_DISTANCE)
    #model_part.AddNodalSolutionStepVariable(YP_NITSCHE)
    model_part.AddNodalSolutionStepVariable(PRESSURE_NITSCHE)
    model_part.AddNodalSolutionStepVariable(NITSCHE_ALPHA)
    model_part.AddNodalSolutionStepVariable(NITSCHE_DELTA)
    model_part.AddNodalSolutionStepVariable(TIMEDIS_THETA)
    model_part.AddNodalSolutionStepVariable(CASE_NUMBER)
    model_part.AddNodalSolutionStepVariable(ELEMENT_ID)
    model_part.AddNodalSolutionStepVariable(PROJECTED_DISTANCE)
    model_part.AddNodalSolutionStepVariable(PROJECTED_DELTA_DISTANCE)
    model_part.AddNodalSolutionStepVariable(LUMPED_MASS_VALUE)
    #model_part.AddNodalSolutionStepVariable(LUMPED_MASS_VALUE_NITSCHE)
    model_part.AddNodalSolutionStepVariable(R_NODE_DISTANCE)
    #model_part.AddNodalSolutionStepVariable(R_NODE_DISTANCE_NITSCHE)
    model_part.AddNodalSolutionStepVariable(SCALARPROJECTEDVEL_X)
    model_part.AddNodalSolutionStepVariable(SCALARPROJECTEDVEL_Y)
    model_part.AddNodalSolutionStepVariable(SCALARPROJECTEDVEL_Z)
    model_part.AddNodalSolutionStepVariable(PRESSURE_N)
    #model_part.AddNodalSolutionStepVariable(DELTA_DISTANCE_NITSCHE)
    model_part.AddNodalSolutionStepVariable(TEMPERATURE)
    ########

    ########
    #model_part.AddNodalSolutionStepVariable(DELTA_VELOCITY_NITSCHE)
    model_part.AddNodalSolutionStepVariable(VELOCITY_NITSCHE)
    #model_part.AddNodalSolutionStepVariable(PROJECTED_VELOCITY_NITSCHE)
    #model_part.AddNodalSolutionStepVariable(VELOCITY_BEFORESOLUTION)
    #model_part.AddNodalSolutionStepVariable(BODY_FORCE_NITSCHE)
    model_part.AddNodalSolutionStepVariable(PROJECTED_DELTA_VELOCITY)
    #model_part.AddNodalSolutionStepVariable(PROJECTED_DELTA_VELOCITY_NITSCHE)
    model_part.AddNodalSolutionStepVariable(R_NODE_VELOCITY)
    #model_part.AddNodalSolutionStepVariable(R_NODE_VELOCITY_NITSCHE)
    model_part.AddNodalSolutionStepVariable(MESH_DELTA_VELOCITY)
    model_part.AddNodalSolutionStepVariable(VELOCITY_PIMPLE)
    model_part.AddNodalSolutionStepVariable(VELOCITY_N)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE_N)
    ########

    #####DEBUGVARIABLES#####
    model_part.AddNodalSolutionStepVariable(ISFOUNDA)
    model_part.AddNodalSolutionStepVariable(ISFOUNDB)
    model_part.AddNodalSolutionStepVariable(ISFOUNDRK4A)
    model_part.AddNodalSolutionStepVariable(ISFOUNDRK4B)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUBA)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUBB)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUB2A)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUB2B)

    model_part.AddNodalSolutionStepVariable(ISFOUNDAFIRST)
    model_part.AddNodalSolutionStepVariable(ISFOUNDBFIRST)
    model_part.AddNodalSolutionStepVariable(ISFOUNDASECOND)
    model_part.AddNodalSolutionStepVariable(ISFOUNDBSECOND)
    model_part.AddNodalSolutionStepVariable(ISFOUNDATHIRD)
    model_part.AddNodalSolutionStepVariable(ISFOUNDBTHIRD)
    model_part.AddNodalSolutionStepVariable(ISFOUNDAFOURTH)
    model_part.AddNodalSolutionStepVariable(ISFOUNDBFOURTH)
    model_part.AddNodalSolutionStepVariable(ISFOUNDRK4AFIRST)
    model_part.AddNodalSolutionStepVariable(ISFOUNDRK4BFIRST)
    model_part.AddNodalSolutionStepVariable(ISFOUNDRK4ASECOND)
    model_part.AddNodalSolutionStepVariable(ISFOUNDRK4BSECOND)
    model_part.AddNodalSolutionStepVariable(ISFOUNDRK4ATHIRD)
    model_part.AddNodalSolutionStepVariable(ISFOUNDRK4BTHIRD)
    model_part.AddNodalSolutionStepVariable(ISFOUNDRK4AFOURTH)
    model_part.AddNodalSolutionStepVariable(ISFOUNDRK4BFOURTH)    
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUBAFIRST)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUBBFIRST)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUBASECOND)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUBBSECOND)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUBATHIRD)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUBBTHIRD)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUBAFOURTH)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUBBFOURTH)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUB2AFIRST)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUB2BFIRST)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUB2ASECOND)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUB2BSECOND)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUB2ATHIRD)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUB2BTHIRD)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUB2AFOURTH)
    model_part.AddNodalSolutionStepVariable(ISFOUNDSUB2BFOURTH)

    model_part.AddNodalSolutionStepVariable(PROJECTED_VELOCITY_FIRST)
    model_part.AddNodalSolutionStepVariable(PROJECTED_VELOCITY_SECOND)
    #####DEBUGVARIABLES#####

def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);
        node.AddDof(DISTANCE);
        #node.AddDof(PRESSURE_NITSCHE);
        #node.AddDof(VELOCITY_NITSCHE_X);
        #node.AddDof(VELOCITY_NITSCHE_Y);
        #node.AddDof(VELOCITY_NITSCHE_Z);


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
            "velocity_block_preconditioner" : {
                "tolerance" : 1e-3,
                "preconditioner_type" : "ilu0"
            },
            "pressure_block_preconditioner" : {
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


        solver_settings = Parameters(

          '''{
             "preconditioner_type"            : "amg",
             "solver_type"                    : "AMGCL",
             "smoother_type"                  : "ilu0",
             "krylov_type"                    : "cg",
             "coarsening_type"                : "aggregation",
             "max_iteration"                  : 1000,
             "verbosity"                      : 1,
             "tolerance"                      : 1e-8,
             "scaling"                        : false,
             "block_size"                     : 1,
             "use_block_matrices_if_possible" : true,
             "coarse_enough"                  : 1000,
             "max_levels"                     : -1,
             "pre_sweeps"                     : 1,
             "post_sweeps"                    : 1
          }''')



        solver_settings = Parameters(

          '''{
            "solver_type":         "BICGSTABSolver",
            "tolerance":           1.0e-8,
            "max_iteration":       300,
            "scaling":             false,
            "preconditioner_type": "ILU0Preconditioner"
          }''')


        solver_settings = Parameters(

          '''{
            "solver_type":         "bicgstab",
            "tolerance":           1.0e-8,
            "max_iteration":       500,
            "scaling":             false,
            "preconditioner_type": "ilu0"
          }''')


        solver_settings = Parameters(

          '''{
             "preconditioner_type"            : "amg",
             "solver_type"                    : "amgcl",
             "smoother_type"                  : "ilu0",
             "krylov_type"                    : "gmres",
             "coarsening_type"                : "aggregation",
             "max_iteration"                  : 1000,
             "verbosity"                      : 2,
             "tolerance"                      : 1e-8,
             "scaling"                        : false,
             "block_size"                     : 1,
             "use_block_matrices_if_possible" : true,
             "coarse_enough"                  : 1000,
             "max_levels"                     : -1,
             "pre_sweeps"                     : 1,
             "post_sweeps"                    : 1
          }''')  
        #construct the solver

        

        #construct the linear solvers
        #import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        #import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
        #self.monolithic_linear_solver = linear_solver_factory.ConstructSolver(solver_settings)
        #self.monolithic_linear_solver =  linear_solver
        #self.monolithic_linear_solver = AMGCLSolver(AMGCLSmoother.ILU0,AMGCLIterativeSolverType.BICGSTAB,tol,1000,verbosity,gmres_size)
        #self.monolithic_linear_solver=SkylineLUFactorizationSolver()
        self.monolithic_linear_solver = linear_solver_factory.ConstructSolver(solver_settings)

        self.conv_criteria = DisplacementCriteria(1e-3,1e-3)  #tolerance for the solver
        self.conv_criteria.SetEchoLevel(0)

        self.domain_size = domain_size
        ##calculate normals
        self.normal_tools = BodyNormalCalculationUtils()

        self.time = 0.0
        self.timestep = 0.0
        self.bdfcalled=False

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
        maximum_number_of_particles= 10*1*self.domain_size
        #maximum_number_of_particles= 8*self.domain_size

        #self.ExplicitStrategy=PFEM2_Explicit_Strategy(self.model_part,self.domain_size, MoveMeshFlag)
        self.VariableUtils = VariableUtils()

        #if self.domain_size==2:
        #     self.moveparticles = MoveParticleUtilityFullBinsPFEM22D(self.model_part,maximum_number_of_particles)
        #else:
        #     self.moveparticles = MoveParticleUtilityFullBinsPFEM23D(self.model_part,maximum_number_of_particles)
        #raw_input()
        print("self.domain_size = ", self.domain_size)
        #if self.domain_size==2:
        #     self.calculatewatervolume = CalculateWaterFraction2D(self.model_part)
        #else:
        #     self.calculatewatervolume = CalculateWaterFraction3D(self.model_part)

        #self.moveparticles.MountBin()
        self.water_volume=0.0  #we initialize it at zero
        self.water_initial_volume=0.0 #we initialize it at zero
        self.mass_correction_factor=0.0

        self.normal_tools.CalculateBodyNormals(self.model_part,self.domain_size);
        condition_number=1

        #if self.domain_size==2:
        #      self.addBC = AddFixedPressureCondition2D(self.model_part)
        #else:
        #      self.addBC = AddFixedPressureCondition3D(self.model_part)

        #(self.addBC).AddThem()

        if self.domain_size==2:
             self.distance_utils = SignedDistanceCalculationUtils2D()
        else:
             self.distance_utils = SignedDistanceCalculationUtils3D()
        self.redistance_step = 0

        import KratosMultiphysics.PFEM2Application.strategy_python as strategy_python #implicit solver
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

        self.print_times=True

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

        print("Simulation Time=",self.time)
        print("BFECC START")
        #mount the search structure
        locator = BinBasedFastPointLocator2D(self.model_part)
        locator.UpdateSearchDatabase()
        #construct the utility to move the points
        bfecc_utility = BFECCConvectionRK42D(locator)
        bfecc_utility.TransferOldVelocityToBFECC(self.model_part.Nodes)
        bfecc_utility.TransferOldVelocityToOldVelocityAux(self.model_part.Nodes)

        
        print("FIRST")
        bfecc_utility.CopyAccComponentXToWatPrAcc(self.model_part.Nodes)
        KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(
            self.model_part,
            SCALARPROJECTEDVEL_X,
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.NODAL_AREA).Execute()

        #levelset_convection_settings = KratosMultiphysics.Parameters("""{
        #    "levelset_variable_name" : "SCALARPROJECTEDVEL_X",
        #    "levelset_convection_variable_name" : "VELOCITY",
        #    "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
        #    "max_CFL" : 1.0,
        #    "max_substeps" : 0,
        #    "levelset_splitting" : false,
        #    "eulerian_error_compensation" : true,
        #    "cross_wind_stabilization_factor" : 0.7
        #}""")
        levelset_convection_settings = KratosMultiphysics.Parameters("""{
            "levelset_variable_name" : "SCALARPROJECTEDVEL_X",
            "levelset_convection_variable_name" : "VELOCITY",
            "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
            "eulerian_error_compensation" : true,
            "element_type" : "levelset_convection_supg"
        }""")
        KratosMultiphysics.LevelSetConvectionProcess2D(
            self.model_part,
            self.monolithic_linear_solver,
            levelset_convection_settings).Execute()
        print("SECOND")
        bfecc_utility.CopyAccComponentYToWatPrAcc(self.model_part.Nodes)
        KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(
            self.model_part,
            SCALARPROJECTEDVEL_Y,
            KratosMultiphysics.DISTANCE_GRADIENT,
            KratosMultiphysics.NODAL_AREA).Execute()

        #levelset_convection_settings = KratosMultiphysics.Parameters("""{
        #    "levelset_variable_name" : "SCALARPROJECTEDVEL_Y",
        #    "levelset_convection_variable_name" : "VELOCITY",
        #    "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
        #    "max_CFL" : 1.0,
        #    "max_substeps" : 0,
        #    "levelset_splitting" : false,
        #    "eulerian_error_compensation" : true,
        #    "cross_wind_stabilization_factor" : 0.7
        #}""")
        levelset_convection_settings = KratosMultiphysics.Parameters("""{
            "levelset_variable_name" : "SCALARPROJECTEDVEL_Y",
            "levelset_convection_variable_name" : "VELOCITY",
            "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
            "eulerian_error_compensation" : true,
            "element_type" : "levelset_convection_supg"
        }""")
        KratosMultiphysics.LevelSetConvectionProcess2D(
            self.model_part,
            self.monolithic_linear_solver,
            levelset_convection_settings).Execute()
        

        bfecc_utility.TransferBFECCToVelocity(self.model_part.Nodes)
        (self.VariableUtils).CopyVectorVar(VELOCITY,PROJECTED_VELOCITY_FIRST,self.model_part.Nodes)
        

        for node in (self.model_part.Nodes):
         if node.X>15.49999999:
          if (node.GetSolutionStepValue(VELOCITY_X,0)<0.0):
           node.SetSolutionStepValue(VELOCITY_X,0,0.00000001)

        full_reset=True;
        bfecc_utility.ResetBoundaryConditions(self.model_part,full_reset)
        bfecc_utility.CopyVectorVarToPreviousTimeStep(VELOCITY,self.model_part.Nodes)
        print("first BFECC END")


        
        print("Stokes Solve START")
        (self.monolithic_solver).Solve() #implicit resolution of the system.
        for node in (self.model_part.Nodes):
         if node.X>15.49999999:
          if (node.GetSolutionStepValue(VELOCITY_X,0)<0.0):
           node.SetSolutionStepValue(VELOCITY_X,0,0.00000001)
        (self.VariableUtils).CopyVectorVar(VELOCITY,PROJECTED_VELOCITY,self.model_part.Nodes)
        print("Stokes Solve END")



        
        bfecc_utility.CalculateAccelerationOnTheMeshFirstOrder(self.model_part)
        #bfecc_utility.CalculateAccelerationOnTheMeshSecondOrder(self.model_part)


         #self.nodaltasks = self.nodaltasks + t11-t9

         #print ". erasing  = ", t14-t13
         #self.total = self.total + t13-t1

         #if self.print_times==True:
         #   print( "----------TIMES----------" )
         #   print( "self.streamlineintegration " ,  self.streamlineintegration , "in % = ", 100.0*(self.streamlineintegration)/(self.total) )
         #   print( "self.accelerateparticles " ,  self.accelerateparticles , "in % = ", 100.0*(self.accelerateparticles)/(self.total) )
         #   print( "self.implicitsystem " ,  self.implicitsystem , "in % = ", 100.0*(self.implicitsystem)/(self.total) )
         #   print( "self.lagrangiantoeulerian " ,  self.lagrangiantoeulerian , "in % = ", 100.0*(self.lagrangiantoeulerian)/(self.total) )
         #   print( "self.reseed " ,  self.reseed , "in % = ", 100.0*(self.reseed)/(self.total) )
         #   print( "self.prereseed " ,  self.prereseed , "in % = ", 100.0*(self.prereseed)/(self.total) )
         #   print( "TOTAL ----- " ,  self.total )
         #   print( "this time step" ,  t13-t1  )
         #   print( "current time in simulation: ", self.model_part.ProcessInfo.GetValue(TIME) , "s" )


    #def RotateParticlesAndDomainVelocities(self,angles):
    #    (self.moveparticles).RotateParticlesAndDomainVelocities(angles)
    #######################################################################
    #def CalculatePressureProjection(self):
        #self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 10)
        #(self.ExplicitStrategy).InitializeSolutionStep();
        #(self.ExplicitStrategy).AssembleLoop();
        #self.model_part.ProcessInfo.SetValue(FRACTIONAL_STEP, 10)
        #(self.ExplicitStrategy).FinalizeSolutionStep();


    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)

    def PrintInfo(self,print_times):
        self.print_times=print_times

    def WriteRestartFile(self,FileName):
        restart_file = open(FileName + ".mdpa",'w')
        import KratosMultiphysics.PFEM2Application.new_restart_utilities as new_restart_utilities
        #import new_restart_utilities
        new_restart_utilities.PrintProperties(restart_file)
        new_restart_utilities.PrintNodes(self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintElements("MonolithicPFEM22DDenizNitsche",self.model_part.Elements,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_X,"VELOCITY_X",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_Y,"VELOCITY_Y",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_Z,"VELOCITY_Z",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(PRESSURE,"PRESSURE",self.model_part.Nodes,restart_file)
        restart_file.close()


    def CopyTimeStep(self,timestep):
        self.timestep=timestep
    def CopyTime(self,time):
        self.time=time
    def ReduceTimeStep(self):
        newdt=0.5*self.timestep
        TimeCoeff = 1.0 / (newdt)
        BDFcoeffs = Vector(3)
        BDFcoeffs[0] = TimeCoeff 
        BDFcoeffs[1] = -TimeCoeff 
        BDFcoeffs[2] = 0.0
        self.model_part.ProcessInfo.SetValue(BDF_COEFFICIENTS, BDFcoeffs)
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, newdt)
    def IncreaseTimeStep(self):
        olddt=self.timestep
        TimeCoeff = 1.0 / (olddt)
        BDFcoeffs = Vector(3)
        BDFcoeffs[0] = TimeCoeff 
        BDFcoeffs[1] = -TimeCoeff 
        BDFcoeffs[2] = 0.0
        self.model_part.ProcessInfo.SetValue(BDF_COEFFICIENTS, BDFcoeffs)
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, olddt)
