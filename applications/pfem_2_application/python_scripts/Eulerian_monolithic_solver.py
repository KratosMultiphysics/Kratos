from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.PFEM2Application import *
#from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

variables_dictionary = {"PRESSURE" : PRESSURE,
                        "VELOCITY" : VELOCITY,
                        "REACTION" : REACTION,
                        "DISTANCE" : DISTANCE,
			 #"AUX_VEL" : AUX_VEL,
                        "DISPLACEMENT" : DISPLACEMENT,
                        "IS_INTERFACE" : IS_INTERFACE,
                        "IS_STRUCTURE" : IS_STRUCTURE,
                        #"VISCOUS_STRESSX": VISCOUS_STRESSX,
                        #"VISCOUS_STRESSY": VISCOUS_STRESSY,
                        #"IS_WATER": IS_WATER,
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
    model_part.AddNodalSolutionStepVariable(EMBEDDED_VELOCITY)
    model_part.AddNodalSolutionStepVariable(IS_POROUS)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(Y_WALL)
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    model_part.AddNodalSolutionStepVariable(IS_FLUID)
    model_part.AddNodalSolutionStepVariable(FRACT_VEL)
    model_part.AddNodalSolutionStepVariable(IS_LAGRANGIAN_INLET)


def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(PRESSURE, REACTION_WATER_PRESSURE)


    print("dofs for the SurfaceTension monolithic solver added correctly")


class STMonolithicSolver:
    def __init__(self, model_part, domain_size):
        self.model_part = model_part
        self.domain_size = domain_size
        # eul_model_part can be 0 (meaning that the model part is lagrangian) or 1 (eulerian)
        #self.eul_model_part = 0

        self.alpha = 0.0
        #if(self.eul_model_part == 0):
        self.move_mesh_strategy = 2
        #else:
        #    self.move_mesh_strategy = 0
        #self.move_mesh_strategy = 2
        # definition of the solvers
        try:
            from KratosMultiphysics.ExternalSolversApplication import SuperLUIterativeSolver
            self.linear_solver = SuperLUIterativeSolver()
        except:
            self.linear_solver = SkylineLUFactorizationSolver()
        


        pDiagPrecond = DiagonalPreconditioner()

        #self.linear_solver = BICGSTABSolver(1e-4, 5000, pDiagPrecond)
        #self.linear_solver = BICGSTABSolver(1e-5, 5000, pILU0)
        pILU0 = ILU0Preconditioner()
        self.linear_solver = BICGSTABSolver(1e-5, 5000, pILU0)
        self.linear_solver = BICGSTABSolver(1e-4, 5000, pILU0)

        #self.linear_solver = SkylineLUFactorizationSolver()
        #self.linear_solver = BICGSTABSolver(1e-6, 5000, pDiagPrecond)


        #self.linear_solver = SkylineLUFactorizationSolver()
        # definition of the convergence criteria
        self.rel_vel_tol = 1e-3
        self.abs_vel_tol = 1e-5 #1e-6
        self.rel_pres_tol = 1e-3
        self.abs_pres_tol = 1e-5 #1e-6
        self.abs_pres_tol = 1e-1 #1e-6
        self.dynamic_tau = 0.0
        self.oss_switch = 0

        # non newtonian setting
        self.regularization_coef = 0
        self.max_iter = 30

        # default settings
        self.echo_level = 0
        self.compute_reactions = True
        self.ReformDofSetAtEachStep = False
        self.CalculateNormDxFlag = False
        self.MoveMeshFlag = False
        self.use_slip_conditions = False

        #self.time_scheme = None
        #self.builder_and_solver = None

        self.turbulence_model = None
        self.use_spalart_allmaras = False
        self.use_des = False
        self.Cdes = 1.0
        self.wall_nodes = list()
        self.spalart_allmaras_linear_solver = None

        self.divergence_clearance_steps = 0

        print("Construction monolithic solver finished")

        print("after reading all the model contains:")
        print(self.model_part)

        #if(self.eul_model_part == 0):
	    
        self.Streamline = Streamline()
        self.Pfem2Utils = Pfem2Utils()
	     
          #self.PfemUtils = PfemUtils()
        self.mark_outer_nodes_process = MarkOuterNodesProcess(model_part);
        self.node_erase_process = NodeEraseProcess(model_part);
        self.add_nodes=True
        self.alpha_shape = 1.4;
            #self.mark_free_surface_process = MarkFreeSurfaceProcess(model_part);

	  #saving the limits of the box (all the nodes external to this will be erased)
        bounding_box_corner1_x = -0.20000e+00
        bounding_box_corner1_y = -0.20000e+00
        bounding_box_corner1_z = -1.00000e+00
        bounding_box_corner2_x = 1.000000001
        bounding_box_corner2_y = 1.000000001
        bounding_box_corner2_z = 1.01000e+10
        box_corner1 = Vector(3);
        box_corner1[0]=bounding_box_corner1_x; box_corner1[1]=bounding_box_corner1_y; box_corner1[2]=bounding_box_corner1_z;
        box_corner2 = Vector(3);
        box_corner2[0]=bounding_box_corner2_x; box_corner2[1]=bounding_box_corner2_y; box_corner2[2]=bounding_box_corner2_z;
        self.box_corner1 = box_corner1
        self.box_corner2 = box_corner2

        #if(domain_size == 2):
        #    self.Mesher =  TriGenPFEMModeler()
        #    self.fluid_neigh_finder = FindNodalNeighboursProcess(model_part,9,18)
        #    #this is needed if we want to also store the conditions a node belongs to
        #    self.condition_neigh_finder = FindConditionsNeighboursProcess(model_part,2, 10)
        #    self.elem_neighbor_finder = FindElementalNeighboursProcess(model_part, 2, 10)	
        #    self.Pfem2_apply_bc_process = Pfem2ApplyBCProcess(model_part);
        #    self.mark_fluid_process = MarkFluidProcess(model_part);

        #elif (domain_size == 3):
        #    self.Mesher = TetGenPfemModeler()
        #    self.fluid_neigh_finder = FindNodalNeighboursProcess(model_part,20,30)
        #    #this is needed if we want to also store the conditions a node belongs to
        #    self.condition_neigh_finder = FindConditionsNeighboursProcess(model_part,3, 20)
        #    self.elem_neighbor_finder = FindElementalNeighboursProcess(model_part, 20, 30)	
        #    self.Pfem2_apply_bc_process = Pfem2ApplyBCProcess(model_part);
        #    self.mark_fluid_process = MarkFluidProcess(model_part);

        #(self.fluid_neigh_finder).Execute();
        Hfinder  = FindNodalHProcess(model_part);
        Hfinder.Execute();
            

    #
    def Initialize(self):
        
        # creating the solution strategy
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol, self.abs_vel_tol, self.rel_pres_tol, self.abs_pres_tol)
        #self.conv_criteria = ResidualCriteria(0.0001, 0.0000001)

        #(self.conv_criteria).SetEchoLevel(self.echo_level)

        self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(self.alpha, self.move_mesh_strategy, self.domain_size)

        builder_and_solver = ResidualBasedBlockBuilderAndSolver(
            self.linear_solver)

        self.solver = ResidualBasedNewtonRaphsonStrategy(
            self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria,
            builder_and_solver, self.max_iter, self.compute_reactions, True, self.MoveMeshFlag)
        (self.solver).SetEchoLevel(self.echo_level)
        self.solver.Check()


        #self.Remesh()


# print "Initialization monolithic solver finished"
    #
    def Solve(self):
        
        #self.Streamline.MoveMesh_RKTn(self.model_part,100)

        #self.Remesh();

        #(self.fluid_neigh_finder).Execute();


        #self.Streamline.Force(self.model_part)

        (self.solver).Solve() #it dumps in this line... 20151020
        
        #self.Streamline.MoveMesh_FE(self.model_part,100)

        #self.Streamline.MoveMesh_RK(self.model_part,100)
        
        #self.Remesh();

        #(self.fluid_neigh_finder).Execute();
        
        
 

    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)


    def Clear(self):
        (self.solver).Clear()



    ##########################################
    def Remesh(self):

        for node in (self.model_part).Nodes: 
            node.SetSolutionStepValue(NODAL_H,0,0.025) #fin
            node.SetSolutionStepValue(NODAL_H,0,0.002) #fin

        
        (self.fluid_neigh_finder).Execute();
        (self.elem_neighbor_finder).Execute()
        (self.condition_neigh_finder).Execute();

        (self.Pfem2_apply_bc_process).Execute();

        self.alpha_shape = 1.4;

        #self.Pfem2Utils.MarkNodesTouchingWall(self.model_part, self.domain_size, 0.30)

        h_factor=0.15;
        #(self.fluid_neigh_finder).Execute();   

        (self.mark_outer_nodes_process).MarkOuterNodes(self.box_corner1, self.box_corner2);
        #self.Pfem2Utils.MarkLonelyNodesForErasing(self.model_part)

        #self.node_erase_process.Execute()

        #LagrangianFluidVMS2D VMS2D
#        (self.Pfem2_apply_bc_process).Execute();
        if (self.domain_size == 2):
            (self.Mesher).ReGenerateMesh("LagrangianFluidVMS2D","Condition2D", self.model_part, self.node_erase_process, True, False, self.alpha_shape, h_factor)

        elif (self.domain_size == 3):
            print("ddddddddd")
             
            (self.Mesher).ReGenerateMesh("LagrangianFluidVMS3D","Condition3D", self.model_part, self.node_erase_process, True, False, self.alpha_shape, h_factor)    

        #LagrangianFluidVMS2D


        (self.fluid_neigh_finder).Execute();
        (self.elem_neighbor_finder).Execute()
        (self.condition_neigh_finder).Execute();
        (self.Pfem2_apply_bc_process).Execute();


        ##############THIS IS FOR EMBEDDED"""""""""""""""""""""""""
        print("end of remesh function")
    ######################################################################

    def FindNeighbours(self):
        (self.neigh_finder).Execute();


def CreateSolver(model_part, config): #FOR 3D!
    fluid_solver = STMonolithicSolver(model_part, config.domain_size)
    
    if(hasattr(config, "alpha")):
        fluid_solver.alpha = config.alpha

    #if(hasattr(config, "eul_model_part")):
    #    fluid_solver.eulerian_model_part = config.eul_model_part

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

    #import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
    #if(hasattr(config, "linear_solver_config")):
    #    fluid_solver.linear_solver = linear_solver_factory.ConstructSolver(
    #        config.linear_solver_config)

    return fluid_solver
