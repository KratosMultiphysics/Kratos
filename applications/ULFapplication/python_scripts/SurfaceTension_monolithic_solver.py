from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
#from KratosMultiphysics.PFEM2Application import Pfem2Utils
from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.MeshingApplication import *
import math
import time
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


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
    model_part.AddNodalSolutionStepVariable(NODAL_LENGTH)
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
    
    #20150128 ajarauta: surface tension variables
    model_part.AddNodalSolutionStepVariable(CURVATURE)
    model_part.AddNodalSolutionStepVariable(IS_WATER)
    model_part.AddNodalSolutionStepVariable(VISCOUS_STRESSX)
    model_part.AddNodalSolutionStepVariable(VISCOUS_STRESSY)
    model_part.AddNodalSolutionStepVariable(VISCOUS_STRESSZ)
    model_part.AddNodalSolutionStepVariable(TRIPLE_POINT)
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    model_part.AddNodalSolutionStepVariable(IS_FLUID)
    model_part.AddNodalSolutionStepVariable(CONTACT_ANGLE)
    model_part.AddNodalSolutionStepVariable(FRACT_VEL)
    model_part.AddNodalSolutionStepVariable(IS_LAGRANGIAN_INLET)
    model_part.AddNodalSolutionStepVariable(MEAN_CURVATURE)
    model_part.AddNodalSolutionStepVariable(GAUSSIAN_CURVATURE)
    model_part.AddNodalSolutionStepVariable(PRINCIPAL_CURVATURE_1)
    model_part.AddNodalSolutionStepVariable(PRINCIPAL_CURVATURE_2)
    model_part.AddNodalSolutionStepVariable(PRINCIPAL_DIRECTION_1)
    model_part.AddNodalSolutionStepVariable(PRINCIPAL_DIRECTION_2)
    model_part.AddNodalSolutionStepVariable(NORMAL_GEOMETRIC)
    model_part.AddNodalSolutionStepVariable(ADHESION_FORCE)
    model_part.AddNodalSolutionStepVariable(NORMAL_CL)
    model_part.AddNodalSolutionStepVariable(NORMAL_CL_EQ)
    model_part.AddNodalSolutionStepVariable(NORMAL_EQ)
    model_part.AddNodalSolutionStepVariable(NORMAL_TP)
    model_part.AddNodalSolutionStepVariable(SOLID_FRACTION_GRADIENT)
    model_part.AddNodalSolutionStepVariable(PHASE_FRACTION_GRADIENT)

    if config is not None:
        if hasattr(config, "TurbulenceModel"):
            if config.TurbulenceModel == "Spalart-Allmaras":
                model_part.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY)
                model_part.AddNodalSolutionStepVariable(MOLECULAR_VISCOSITY)
                model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
                model_part.AddNodalSolutionStepVariable(DISTANCE)

    print("variables for the SurfaceTension monolithic solver added correctly")


def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(PRESSURE, REACTION_WATER_PRESSURE)

    if config is not None:
        if hasattr(config, "TurbulenceModel"):
            if config.TurbulenceModel == "Spalart-Allmaras":
                for node in model_part.Nodes:
                    node.AddDof(TURBULENT_VISCOSITY)

    print("dofs for the surface tension monolithic solver added correctly")


class MonolithicSolver:
    def __init__(self, model_part, domain_size, eul_model_part, gamma, contact_angle):
        self.model_part = model_part
        self.domain_size = domain_size
        # eul_model_part can be 0 (meaning that the model part is lagrangian) or 1 (eulerian)
        self.eul_model_part = eul_model_part

        self.alpha = -0.3
        if(eul_model_part == 0):
            self.move_mesh_strategy = 2
        else:
            self.move_mesh_strategy = 0

        # definition of the solvers
        try:
            from KratosMultiphysics.ExternalSolversApplication import SuperLUIterativeSolver
            self.linear_solver = SuperLUIterativeSolver()
        except:
            self.linear_solver = SkylineLUFactorizationSolver()

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
        
        self.contact_angle = contact_angle
        self.gamma = gamma        

        # default settings
        self.echo_level = 0
        self.compute_reactions = True
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = True
        self.use_slip_conditions = False

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
        
        if(self.eul_model_part == 0):
          self.UlfUtils = UlfUtils()
          #self.Pfem2Utils = Pfem2Utils()
          self.mark_outer_nodes_process = MarkOuterNodesProcess(model_part);
          self.node_erase_process = NodeEraseProcess(model_part);
          self.add_nodes=True
          self.alpha_shape = 3.5;
          self.ulf_apply_bc_process = UlfApplyBCProcess(model_part);
	  #self.mark_fluid_process = MarkFluidProcess(model_part);
	  
	  #saving the limits of the box (all the nodes external to this will be erased)
          bounding_box_corner1_x = -1.00000e+00
          bounding_box_corner1_y = -1.00000e+00
          bounding_box_corner1_z = -1.00000e+00
          bounding_box_corner2_x = 1.01000e+10
          bounding_box_corner2_y = 1.01000e+10
          bounding_box_corner2_z = 1.01000e+10
          box_corner1 = Vector(3); 
          box_corner1[0]=bounding_box_corner1_x; box_corner1[1]=bounding_box_corner1_y; box_corner1[2]=bounding_box_corner1_z;
          box_corner2 = Vector(3); 
          box_corner2[0]=bounding_box_corner2_x; box_corner2[1]=bounding_box_corner2_y; box_corner2[2]=bounding_box_corner2_z;
          self.box_corner1 = box_corner1
          self.box_corner2 = box_corner2
	  
          if(domain_size == 2):
              self.Mesher = TriGenPFEMModeler()            
              self.fluid_neigh_finder = FindNodalNeighboursProcess(model_part,9,18)
	      #this is needed if we want to also store the conditions a node belongs to
              self.condition_neigh_finder = FindConditionsNeighboursProcess(model_part,2, 10)
          #elif (domain_size == 3):
              #self.Mesher = TetGenPfemModeler()            
              #self.fluid_neigh_finder = FindNodalNeighboursProcess(model_part,20,30)
	      ##this is needed if we want to also store the conditions a node belongs to
              #self.condition_neigh_finder = FindConditionsNeighboursProcess(model_part,3, 20)
	      #self.surf_neigh_finder = FindNodalNeighboursProcess(surface_part,20,30)
	      #self.surf_cond_neigh_finder = FindConditionsNeighboursProcess(surface_part,3, 20)
	      
          (self.fluid_neigh_finder).Execute();
	  #(self.condition_neigh_finder).Execute();
          Hfinder  = FindNodalHProcess(model_part);
          Hfinder.Execute();

    #
    def Initialize(self):
        # check if slip conditions are defined
        if self.use_slip_conditions == False:
            for cond in self.model_part.Conditions:
                if cond.GetValue(IS_STRUCTURE) != 0.0:
                    self.use_slip_conditions = True
                    break

        # if we use slip conditions, calculate normals on the boundary
        if self.use_slip_conditions:
            self.normal_util = NormalCalculationUtils()
            self.normal_util.CalculateOnSimplex(
                self.model_part, self.domain_size, IS_STRUCTURE)

            for cond in self.model_part.Conditions:
                if cond.GetValue(IS_STRUCTURE) != 0.0:
                    for node in cond.GetNodes():
                        node.SetValue(IS_STRUCTURE, 1.0)

        # Turbulence model
        if self.use_spalart_allmaras:
            self.activate_spalart_allmaras()

        # creating the solution strategy
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol, self.abs_vel_tol, self.rel_pres_tol, self.abs_pres_tol)
        #self.conv_criteria = ResidualCriteria(0.0001, 0.0000001)

        if self.turbulence_model is None:
            self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent\
                (self.alpha, self.move_mesh_strategy,
                 self.domain_size)
        else:
            self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent\
                (self.alpha,
                 self.move_mesh_strategy,
                 self.domain_size,
                 self.turbulence_model)

        builder_and_solver = ResidualBasedBlockBuilderAndSolver(
            self.linear_solver)

        self.solver = ResidualBasedNewtonRaphsonStrategy(
            self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria,
            builder_and_solver, self.max_iter, self.compute_reactions, self.ReformDofSetAtEachStep, self.MoveMeshFlag)
        (self.solver).SetEchoLevel(self.echo_level)
        self.solver.Check()

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau)
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch)
        self.model_part.ProcessInfo.SetValue(M, self.regularization_coef)
        
        #self.model_part.ProcessInfo.SetValue(CONTACT_ANGLE_STATIC, self.contact_angle)
        self.model_part.ProcessInfo.SetValue(SURFTENS_COEFF, self.gamma)
        
        if(self.eul_model_part == 0):
	    #saving the structural elements
	    #(self.mark_fluid_process).Execute(); #we need this before saving the structrural elements
        
	    #marking the fluid
            (self.fluid_neigh_finder).Execute();
            (self.ulf_apply_bc_process).Execute();  
	    #(self.mark_fluid_process).Execute();
            #if (self.domain_size == 2):
                #FindTriplePoint().FindTriplePoint2D(self.model_part)
            self.Remesh()

# print "Initialization monolithic solver finished"
    #
    def Solve(self):  
        
      
        if self.divergence_clearance_steps > 0:
            print("Calculating divergence-free initial condition")
            # initialize with a Stokes solution step
            try:
                import KratosMultiphysics.ExternalSolversApplication as kes
                smoother_type = kes.AMGCLSmoother.DAMPED_JACOBI
                solver_type = kes.AMGCLIterativeSolverType.CG
                gmres_size = 50
                max_iter = 200
                tol = 1e-6
                verbosity = 0
                stokes_linear_solver = kes.AMGCLSolver(
                    smoother_type,
                    solver_type,
                    tol,
                    max_iter,
                    verbosity,
                    gmres_size)
            except:
                pPrecond = DiagonalPreconditioner()
                stokes_linear_solver = BICGSTABSolver(1e-6, 5000, pPrecond)
            stokes_process = StokesInitializationProcess(
                self.model_part,
                stokes_linear_solver,
                self.domain_size,
                PATCH_INDEX)
            # copy periodic conditions to Stokes problem
            stokes_process.SetConditions(self.model_part.Conditions)
            # execute Stokes process
            stokes_process.Execute()
            stokes_process = None
            
            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(PRESSURE, 0, 0.0)
                node.SetSolutionStepValue(ACCELERATION_X, 0, 0.0)
                node.SetSolutionStepValue(ACCELERATION_Y, 0, 0.0)
                node.SetSolutionStepValue(ACCELERATION_Z, 0, 0.0)

            self.divergence_clearance_steps = 0
            print("Finished divergence clearance")

        if self.ReformDofSetAtEachStep:
            if self.use_slip_conditions:
                self.normal_util.CalculateOnSimplex(
                    self.model_part, self.domain_size, IS_STRUCTURE)
            if self.use_spalart_allmaras:
                self.neighbour_search.Execute()

        NormalCalculationUtils().CalculateOnSimplex(self.model_part.Conditions, self.domain_size)
        for node in self.model_part.Nodes:
            if (node.GetSolutionStepValue(IS_BOUNDARY) == 0.0):# and node.GetSolutionStepValue(TRIPLE_POINT) == 0):
               node.SetSolutionStepValue(NORMAL_X,0,0.0)
               node.SetSolutionStepValue(NORMAL_Y,0,0.0)
               node.SetSolutionStepValue(NORMAL_Z,0,0.0)     
        if (self.domain_size == 2):
                #FindTriplePoint().FindTriplePoint2D(self.model_part)
            CalculateCurvature().CalculateCurvature2D(self.model_part)
            #CalculateNodalLength().CalculateNodalLength2D(self.model_part)
            #CalculateContactAngle().CalculateContactAngle2D(self.model_part)
            #self.cont_angle_cond()
        #elif (self.domain_size == 3):
            #FindTriplePoint().FindTriplePoint3D(self.model_part)
            #CalculateCurvature().CalculateCurvature3D(self.model_part)
            #for node in self.model_part.Nodes:
            #node.SetSolutionStepValue(CONTACT_ANGLE,0,0.0)
            #CalculateContactAngle().CalculateContactAngle3D(self.model_part)
            #CalculateNodalLength().CalculateNodalLength3D(self.model_part)
            #CalculateNormalEq().CalculateNormalEq3D(self.model_part)
                ##CalculateAdhesionForce().CalculateAdhesionForce3D(self.model_part)
            #self.cont_angle_cond3D()
	  
        (self.solver).Solve() #it dumps in this line... 20151020
        #AssignPointNeumannConditions().AssignPointNeumannConditions3D(self.model_part)
        
        if(self.eul_model_part == 0):
            (self.fluid_neigh_finder).Execute();
	    #(self.condition_neigh_finder).Execute(); #20160105 Commented because in 3D it dumps here
            (self.fluid_neigh_finder).Execute();
            self.Remesh();

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
    def Clear(self):
        (self.solver).Clear()

    #
    def activate_smagorinsky(self, C):
        for elem in self.model_part.Elements:
            elem.SetValue(C_SMAGORINSKY, C)

    #
    def activate_spalart_allmaras(self):

        # Spalart-Allmaras initialization
        for node in self.wall_nodes:
            node.SetValue(IS_VISITED, 1.0)
            node.SetSolutionStepValue(DISTANCE, 0, 0.0)

        if(self.domain_size == 2):
            self.redistance_utils = ParallelDistanceCalculator2D()
        else:
            self.redistance_utils = ParallelDistanceCalculator3D()

        max_levels = 100
        max_distance = 1000
        self.redistance_utils.CalculateDistancesLagrangianSurface(
            self.model_part,
            DISTANCE,
            NODAL_AREA,
            max_levels,
            max_distance)

        # Spalart-Allmaras uses the componentwise builder and solver, which
        # requires a neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(
            self.model_part, number_of_avg_elems, number_of_avg_nodes)
        self.neighbour_search.Execute()

        non_linear_tol = 0.001
        max_it = 10
        reform_dofset = self.ReformDofSetAtEachStep
        time_order = 2

        if self.spalart_allmaras_linear_solver is None:
            pPrecond = DiagonalPreconditioner()
            self.spalart_allmaras_linear_solver = BICGSTABSolver(
                1e-6, 5000, pPrecond)

        self.turbulence_model = SpalartAllmarasTurbulenceModel(
            self.model_part,
            self.spalart_allmaras_linear_solver,
            self.domain_size,
            non_linear_tol,
            max_it,
            reform_dofset,
            time_order)
        if self.use_des:
            self.turbulence_model.ActivateDES(self.Cdes)


    ##########################################
    def Remesh(self):
       
        #self.Pfem2Utils.MarkNodesTouchingWall(self.model_part, self.domain_size, 0.08)
        
        ##erase all conditions and elements prior to remeshing
        ((self.model_part).Elements).clear();
        #((self.model_part).Conditions).clear();      


        (self.mark_outer_nodes_process).MarkOuterNodes(self.box_corner1, self.box_corner2);

        h_factor=0.25;
        
        #remesh CHECK for 3D or 2D
        if (self.domain_size == 2):
             (self.Mesher).ReGenerateMesh("SurfaceTension2D","Condition2D", self.model_part, self.node_erase_process, True, True, self.alpha_shape, h_factor)
        #elif (self.domain_size == 3):
            #(self.Mesher).ReGenerateMesh("SurfaceTension3D","Condition3D", self.model_part, self.node_erase_process, True, False, self.alpha_shape, h_factor)            

        ##calculating fluid neighbours before applying boundary conditions
        (self.fluid_neigh_finder).Execute();
        (self.condition_neigh_finder).Execute();

        #print "marking fluid" and applying fluid boundary conditions
        (self.ulf_apply_bc_process).Execute();  
        #(self.mark_fluid_process).Execute();
	
        (self.UlfUtils).CalculateNodalArea(self.model_part,self.domain_size);
        self.UlfUtils.MarkLonelyNodesForErasing(self.model_part)

        self.node_erase_process.Execute()
        #(self.mark_fluid_process).Execute();

        self.UlfUtils.MarkLonelyNodesForErasing(self.model_part)#, self.domain_size)
        self.node_erase_process.Execute()

        ##############THIS IS FOR EMBEDDED"""""""""""""""""""""""""
######################################################################################################                
        #FOR PFEM
        NormalCalculationUtils().CalculateOnSimplex(self.model_part.Conditions, self.domain_size)
        if (self.domain_size == 2):
          for node in self.model_part.Nodes:
              node.SetSolutionStepValue(FLAG_VARIABLE,0,0.0)
#
          for node in self.model_part.Nodes:
              if (node.GetSolutionStepValue(IS_FREE_SURFACE) > 0.99999999):
                  node.SetSolutionStepValue(FLAG_VARIABLE, 0, 1.0)
                  node.SetSolutionStepValue(IS_INTERFACE, 0, 1.0)
              else:
                    node.SetSolutionStepValue(IS_INTERFACE, 0, 0.0)
          for node in self.model_part.Nodes:
            if (node.GetSolutionStepValue(IS_BOUNDARY) == 0.0):# and node.GetSolutionStepValue(TRIPLE_POINT) == 0):
                  node.SetSolutionStepValue(NORMAL_X,0,0.0)
                  node.SetSolutionStepValue(NORMAL_Y,0,0.0)
                  node.SetSolutionStepValue(NORMAL_Z,0,0.0)	      
           #AssignPointNeumannConditions().AssignPointNeumannConditions3D(self.model_part)
          #FindTriplePoint().FindTriplePoint2D(self.model_part)
          CalculateCurvature().CalculateCurvature2D(self.model_part)
          #CalculateNodalLength().CalculateNodalLength2D(self.model_part)
          #CalculateContactAngle().CalculateContactAngle2D(self.model_part)
          #self.cont_angle_cond()
	  #self.print_pressure()
        #elif (self.domain_size == 3):
          #FindTriplePoint().FindTriplePoint3D(self.model_part)
	  #for node in self.model_part.Nodes:
	    #if (node.GetSolutionStepValue(IS_BOUNDARY) > 0.99999999 and node.GetSolutionStepValue(IS_STRUCTURE) == 0.0):
	      ##node.SetSolutionStepValue(IS_FREE_SURFACE, 0, 1.0)
	      ##node.SetSolutionStepValue(IS_INTERFACE, 0, 1.0)
	  ##AssignPointNeumannConditions().AssignPointNeumannConditions3D(self.model_part)
          #CalculateCurvature().CalculateCurvature3D(self.model_part)
          #for node in self.model_part.Nodes:
              #node.SetSolutionStepValue(CONTACT_ANGLE,0,0.0)
          #CalculateContactAngle().CalculateContactAngle3D(self.model_part)
          #CalculateNodalLength().CalculateNodalLength3D(self.model_part)
          #CalculateNormalEq().CalculateNormalEq3D(self.model_part)
	  ##CalculateAdhesionForce().CalculateAdhesionForce3D(self.model_part)
          #self.cont_angle_cond3D()

        ##############THIS IS FOR EMBEDDED"""""""""""""""""""""""""
        #print("end of remesh function")
    #######################################################################
    #def FindNeighbours(self):
        #(self.neigh_finder).Execute();
        
    #def cont_angle_cond(self):
        #theta_adv = 105
        #theta_rec = 70
	##theta_adv = self.contact_angle + 0.5
	##theta_rec = self.contact_angle - 0.5
        #time = self.model_part.ProcessInfo.GetValue(TIME)
        #dt = self.model_part.ProcessInfo.GetValue(DELTA_TIME)
	##x_mean = 0.0
	##found_tp = 0
	################### For sessile drop examples
        #for node in self.model_part.Nodes:
            #if (node.GetSolutionStepValue(TRIPLE_POINT) != 0.0):
                #if ((node.GetSolutionStepValue(CONTACT_ANGLE) > theta_adv) or (node.GetSolutionStepValue(CONTACT_ANGLE) < theta_rec)):
                    #node.Free(VELOCITY_X)
                #else:
                    #node.SetSolutionStepValue(VELOCITY_X,0, 0.0)
                    #node.Fix(VELOCITY_X)
            #if ((node.GetSolutionStepValue(TRIPLE_POINT) == 0.0) and (node.GetSolutionStepValue(IS_STRUCTURE) != 0.0)):
                    #node.SetSolutionStepValue(VELOCITY_X,0, 0.0)
                    #node.SetSolutionStepValue(VELOCITY_Y,0, 0.0)
                    #node.Fix(VELOCITY_X)
                    #node.Fix(VELOCITY_Y)     
 
    #def cont_angle_cond3D(self):
        #theta_adv = self.contact_angle + 5.0
        #theta_rec = self.contact_angle - 5.0
        #time = self.model_part.ProcessInfo.GetValue(TIME)
        #dt = self.model_part.ProcessInfo.GetValue(DELTA_TIME)
	################### For sessile drop examples
        #if (time < 2*dt):
            #for node in self.model_part.Nodes:
                #if (node.GetSolutionStepValue(TRIPLE_POINT) != 0.0):
                    #node.SetSolutionStepValue(VELOCITY_X,0, 0.0)
                    #node.SetSolutionStepValue(VELOCITY_Y,0, 0.0)
                    #node.SetSolutionStepValue(VELOCITY_Z,0, 0.0)
                    #node.Fix(VELOCITY_X)
                    #node.Fix(VELOCITY_Y)
                    #node.Fix(VELOCITY_Z)	  
        #else:
            #for node in self.model_part.Nodes:
                #if (node.GetSolutionStepValue(TRIPLE_POINT) != 0.0):
                    #node.SetSolutionStepValue(VELOCITY_Z,0,0.0)
                    #node.Fix(VELOCITY_Z)
                    #if ((node.GetSolutionStepValue(CONTACT_ANGLE) > theta_adv) or (node.GetSolutionStepValue(CONTACT_ANGLE) < theta_rec)):		  
                        #node.Free(VELOCITY_X)
                        #node.Free(VELOCITY_Y)
                    #else:
                        #node.SetSolutionStepValue(VELOCITY_X,0, 0.0)
                        #node.SetSolutionStepValue(VELOCITY_Y,0, 0.0)
                        #node.Fix(VELOCITY_X)
                        #node.Fix(VELOCITY_Y)
                    #if (node.Z < -0.00000001):
                        #node.SetValue(TO_ERASE, true)
                #if ((node.GetSolutionStepValue(TRIPLE_POINT) == 0.0) and (node.GetSolutionStepValue(IS_STRUCTURE) != 0.0)):
                        #node.SetSolutionStepValue(VELOCITY_X,0, 0.0)
                        #node.SetSolutionStepValue(VELOCITY_Y,0, 0.0)
                        #node.SetSolutionStepValue(VELOCITY_Z,0, 0.0)
                        #node.Fix(VELOCITY_X)
                        #node.Fix(VELOCITY_Y)
                        #node.Fix(VELOCITY_Z)


def CreateSolver(model_part, config, eul_model_part, gamma, contact_angle): #FOR 3D!!!!!!!!!!
    fluid_solver = MonolithicSolver(model_part, config.domain_size, eul_model_part, gamma, contact_angle)    
      
    if(hasattr(config, "alpha")):
        fluid_solver.alpha = config.alpha

    if(hasattr(config, "eul_model_part")):
        fluid_solver.eulerian_model_part = config.eul_model_part
    
        
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
    if(hasattr(config, "oss_switch")):
        fluid_solver.oss_switch = config.oss_switch
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
    if hasattr(config, "TurbulenceModel"):
        if config.TurbulenceModel == "Spalart-Allmaras":
            fluid_solver.use_spalart_allmaras = True
        elif config.TurbulenceModel == "Smagorinsky-Lilly":
            if hasattr(config, "SmagorinskyConstant"):
                fluid_solver.activate_smagorinsky(config.SmagorinskyConstant)
            else:
                msg = """Fluid solver error: Smagorinsky model requested, but
                         the value for the Smagorinsky constant is
                         undefined."""
                raise Exception(msg)

    import linear_solver_factory
    if(hasattr(config, "linear_solver_config")):
        fluid_solver.linear_solver = linear_solver_factory.ConstructSolver(
            config.linear_solver_config)

    return fluid_solver
