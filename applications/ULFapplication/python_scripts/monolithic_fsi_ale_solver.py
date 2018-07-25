from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
#from KratosMultiphysics.SolidMechanicsApplication import *
#from KratosMultiphysics.ConstitutiveModelsApplication import *


from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
#from KratosConstitutiveModelsApplication import *

#import KratosMultiphysics.ConstitutiveModelsApplication as CMApp
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

variables_dictionary = {"PRESSURE" : PRESSURE,
                        "VELOCITY" : VELOCITY,
                        "REACTION" : REACTION,
                        "DISTANCE" : DISTANCE,
			 "AUX_VEL" : AUX_VEL,                        
                        "DISPLACEMENT" : DISPLACEMENT,
                        "IS_INTERFACE" : IS_INTERFACE,
                        "IS_STRUCTURE" : IS_STRUCTURE,
                        "VISCOUS_STRESSX": VISCOUS_STRESSX,
                        "VISCOUS_STRESSY": VISCOUS_STRESSY,
                        "IS_WATER": IS_WATER,
                        "DENSITY": DENSITY,
                        "VISCOSITY": VISCOSITY}


def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(NODAL_H)
    #model_part.AddNodalSolutionStepVariable(NODAL_LENGTH)
    model_part.AddNodalSolutionStepVariable(ADVPROJ)
    model_part.AddNodalSolutionStepVariable(DIVPROJ)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(Y_WALL)
    model_part.AddNodalSolutionStepVariable(PATCH_INDEX)
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    model_part.AddNodalSolutionStepVariable(IS_FLUID)
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    model_part.AddNodalSolutionStepVariable(BULK_MODULUS)
    model_part.AddNodalSolutionStepVariable(CAUCHY_STRESS_TENSOR)    
    #CHECK! FOR SOME STRANGE REASON THE SOLVER DOESNT WORK IF TRIPLE_POINT is not added, even though it is not used anywhere
    #model_part.AddNodalSolutionStepVariable(TRIPLE_POINT)



def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(PRESSURE, REACTION_WATER_PRESSURE)


    print("dofs for the ALE FSI monolithic solver added correctly")


class FsiAleMonolithicSolver:
    def __init__(self, model_part, domain_size):
        self.model_part = model_part        
        self.domain_size = domain_size
        # eul_model_part can be 0 (meaning that the model part is lagrangian) or 1 (eulerian)
        

        self.alpha = -0.3
        #2 Lagr 0 Eul #1 ALE
        #self.move_mesh_strategy = 2
        self.move_mesh_strategy = 1

        # definition of the convergence criteria
        self.rel_vel_tol = 1e-3
        self.abs_vel_tol = 1e-6
        self.rel_pres_tol = 1e-3
        self.abs_pres_tol = 1e-6
        self.dynamic_tau = 0.0
        self.oss_switch = 0

        # non newtonian setting
        self.regularization_coef = 1000
        self.max_iter = 30
        

        # default settings
        self.echo_level = 0
        self.compute_reactions = False
        self.ReformDofSetAtEachStep = False
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = True
        self.use_slip_conditions = False

        #self.time_scheme = None
        #self.builder_and_solver = None

        self.turbulence_model = None
        self.use_spalart_allmaras = False
        self.use_des = False
        self.Cdes = 1.0
        self.wall_nodes = list()
        #self.spalart_allmaras_linear_solver = None
        
        pDiagPrecond = DiagonalPreconditioner()
        self.linear_solver = BICGSTABSolver(1e-6, 5000, pDiagPrecond)

        self.divergence_clearance_steps = 0

        print("Construction monolithic solver finished")

        print("after reading all the model contains:")
        print(self.model_part)

        self.UlfUtils = UlfUtils()
        #self.PfemUtils = PfemUtils()
        self.mark_outer_nodes_process = MarkOuterNodesProcess(model_part);
        self.node_erase_process = NodeEraseProcess(self.model_part);

        if (self.domain_size==2):
            self.neigh_finder = FindNodalNeighboursProcess(self.model_part,9,18)
        elif (domain_size == 3):
            self.neigh_finder = FindNodalNeighboursProcess(self.model_part,20,30)
        #this is needed if we want to also store the conditions a node belongs to
        #self.cond_neigh_finder = FindConditionsNeighboursProcess(self.model_part,2, 10)

        #(self.fluid_neigh_finder).Execute();
        Hfinder  = FindNodalHProcess(self.model_part);
        Hfinder.Execute();
        #for computing pressure in hypoelastic element
        #self.pressure_calculate_process = PressureCalculateProcess(self.model_part, self.domain_size);
        self.hypoelastic_solid_stress_tensor_calculate_process=HypoelasticStressCalculateProcess(self.model_part, self.domain_size)
        
    #
    def Initialize(self):
                
        #self.Remesh()
        # creating the solution strategy
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol, self.abs_vel_tol, self.rel_pres_tol, self.abs_pres_tol)
        #self.conv_criteria = ResidualCriteria(0.0001, 0.0000001)


        #self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(self.alpha, self.move_mesh_strategy, self.domain_size)
        self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeAleFsi(self.alpha, self.move_mesh_strategy, self.domain_size)
        
        builder_and_solver = ResidualBasedBlockBuilderAndSolver(
            self.linear_solver)

        self.solver = ResidualBasedNewtonRaphsonStrategy(
            self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria,
            builder_and_solver, self.max_iter, self.compute_reactions, self.ReformDofSetAtEachStep, self.MoveMeshFlag)
        (self.solver).SetEchoLevel(self.echo_level)
        self.solver.Check()   

        ##########################################################
        #self.Remesh()        
        #(self.append_model_part_process).AppendPart(self.model_part, solid_model_part);        
        (self.neigh_finder).Execute();       
        #we need normals to prescribe the inlet velocity
        self.normal_util = NormalCalculationUtils()
        
        self.normal_util.CalculateOnSimplex(self.model_part.Conditions, self.domain_size)
        #NormalCalculationUtils().CalculateOnSimplex(self.fluid_model_part.Conditions, self.domain_size)     
        #initializes Cachy stress to zero
        self.hypoelastic_solid_stress_tensor_calculate_process.Execute()
        print("Lalalal")
# print "Initialization monolithic solver finished"
    #
    def Solve(self):

        (self.solver).Solve() #it dumps in this line... 20151020
      
        (self.UlfUtils).CalculateNodalArea(self.model_part, self.domain_size);
        #self.pressure_calculate_process.Execute()
        self.hypoelastic_solid_stress_tensor_calculate_process.Execute()
        print("Lalalal222222222222222")



    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)


    def Clear(self):
        (self.solver).Clear()

    def FindNeighbours(self):
        (self.neigh_finder).Execute();



def CreateSolver(model_part, config): #FOR 3D!
    fluid_solver = FsiAleMonolithicSolver(model_part, config.domain_size)

    if(hasattr(config, "alpha")):
        fluid_solver.alpha = config.alpha
 
    # definition of the convergence criteria
    if(hasattr(config, "velocity_relative_tolerance")):
        fluid_solver.rel_vel_tol = config.velocity_relative_tolerance
    if(hasattr(config, "velocity_absolute_tolerance")):
        fluid_solver.abs_vel_tol = config.velocity_absolute_tolerance
    if(hasattr(config, "pressure_relative_tolerance")):
        fluid_solver.rel_pres_tol = config.pressure_relative_tolerance
    if(hasattr(config, "pressure_absolute_tolerance")):
        fluid_solver.abs_pres_tol = config.pressure_absolute_tolerance
    if(hasattr(config, "dynamic_tau")):
        fluid_solver.dynamic_tau = config.dynamic_tau
    if(hasattr(config, "max_iteration")):
        fluid_solver.max_iter = config.max_iteration
    if(hasattr(config, "echo_level")):
        fluid_solver.echo_level = config.echo_level
    if(hasattr(config, "compute_reactions")):
        fluid_solver.compute_reactions = config.compute_reactions
    if(hasattr(config, "ReformDofSetAtEachStep")):
        fluid_solver.ReformDofSetAtEachStep = config.ReformDofSetAtEachStep
    if(hasattr(config, "divergence_cleareance_step")):
        fluid_solver.divergence_clearance_steps = config.divergence_cleareance_step

    #import linear_solver_factory
    #if(hasattr(config, "linear_solver_config")):
    #    fluid_solver.linear_solver = linear_solver_factory.ConstructSolver(
    #        config.linear_solver_config)

    return fluid_solver
