from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ULFApplication import *
#from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.MeshingApplication import *
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
    model_part.AddNodalSolutionStepVariable(MEAN_CURVATURE_2D)
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
    model_part.AddNodalSolutionStepVariable(MEAN_CURVATURE_3D)
    model_part.AddNodalSolutionStepVariable(GAUSSIAN_CURVATURE)
    model_part.AddNodalSolutionStepVariable(PRINCIPAL_CURVATURE_1)
    model_part.AddNodalSolutionStepVariable(PRINCIPAL_CURVATURE_2)
    model_part.AddNodalSolutionStepVariable(PRINCIPAL_DIRECTION_1)
    model_part.AddNodalSolutionStepVariable(PRINCIPAL_DIRECTION_2)
    model_part.AddNodalSolutionStepVariable(NORMAL_GEOMETRIC)
    model_part.AddNodalSolutionStepVariable(ADHESION_FORCE)
    model_part.AddNodalSolutionStepVariable(NORMAL_CONTACT_LINE)
    model_part.AddNodalSolutionStepVariable(NORMAL_CONTACT_LINE_EQUILIBRIUM)
    model_part.AddNodalSolutionStepVariable(NORMAL_EQUILIBRIUM)
    model_part.AddNodalSolutionStepVariable(NORMAL_TRIPLE_POINT)
    model_part.AddNodalSolutionStepVariable(PHASE_FRACTION_GRADIENT)
    model_part.AddNodalSolutionStepVariable(INITIAL_MESH_SIZE)
    model_part.AddNodalSolutionStepVariable(DISSIPATIVE_FORCE_COEFF_JM_X)
    model_part.AddNodalSolutionStepVariable(DISSIPATIVE_FORCE_COEFF_BM_X)
    model_part.AddNodalSolutionStepVariable(DISSIPATIVE_FORCE_COEFF_SM_X)
    model_part.AddNodalSolutionStepVariable(DISSIPATIVE_FORCE_COEFF_JM_Y)
    model_part.AddNodalSolutionStepVariable(DISSIPATIVE_FORCE_COEFF_BM_Y)
    model_part.AddNodalSolutionStepVariable(DISSIPATIVE_FORCE_COEFF_SM_Y)
    model_part.AddNodalSolutionStepVariable(DISSIPATIVE_FORCE_COEFF_JM_Z)
    model_part.AddNodalSolutionStepVariable(DISSIPATIVE_FORCE_COEFF_BM_Z)
    model_part.AddNodalSolutionStepVariable(DISSIPATIVE_FORCE_COEFF_SM_Z)
    model_part.AddNodalSolutionStepVariable(TESTFACTA)
    model_part.AddNodalSolutionStepVariable(TESTFACTB)
    model_part.AddNodalSolutionStepVariable(TESTFACTC)
    model_part.AddNodalSolutionStepVariable(TESTFACTD)
    model_part.AddNodalSolutionStepVariable(TESTFACTE)
    model_part.AddNodalSolutionStepVariable(IS_WALL)
    



def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(PRESSURE, REACTION_WATER_PRESSURE)


    print("dofs for the SurfaceTension monolithic solver added correctly")


class STMonolithicSolver:
    #def __init__(self, model_part, domain_size, eul_model_part, gamma, contact_angle):
    def __init__(self, model_part, domain_size, eul_model_part, gamma, contact_angle, zeta_dissapative_JM_x, zeta_dissapative_BM_x, zeta_dissapative_SM_x, zeta_dissapative_JM_y, zeta_dissapative_BM_y, zeta_dissapative_SM_y, zeta_dissapative_JM_z, zeta_dissapative_BM_z, zeta_dissapative_SM_z, testfacta, testfactb, testfactc, testfactd, testfacte):
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
        ######try:
            ######from KratosMultiphysics.ExternalSolversApplication import SuperLUIterativeSolver
            ######self.linear_solver = SuperLUIterativeSolver()
        ######except:
            ######self.linear_solver = SkylineLUFactorizationSolver()
        pDiagPrecond = DiagonalPreconditioner()
        self.linear_solver = BICGSTABSolver(1e-6, 5000, pDiagPrecond)

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
        
        #dissipative force relalted variables (like true or false to either use them or no): for example, the first one is the JM model in the x direction, and so on... here we are using the value of zero or one to enable which model and in which direction we are going to use our model. For example, We will be using this when we are going to have droplet basement on the XZ plane for example so we will disable the y contribution.
        self.zeta_dissapative_JM_x = zeta_dissapative_JM_x
        self.zeta_dissapative_BM_x = zeta_dissapative_BM_x
        self.zeta_dissapative_SM_x = zeta_dissapative_SM_x
        self.zeta_dissapative_JM_y = zeta_dissapative_JM_y
        self.zeta_dissapative_BM_y = zeta_dissapative_BM_y
        self.zeta_dissapative_SM_y = zeta_dissapative_SM_y
        self.zeta_dissapative_JM_z = zeta_dissapative_JM_z
        self.zeta_dissapative_BM_z = zeta_dissapative_BM_z
        self.zeta_dissapative_SM_z = zeta_dissapative_SM_z
        #these variable are for testing purposes: For example, (testfacta) is either zero or one, this is mainly to see where I am adding the dissipative forces in the surface_tension.h file, and which one will be working fine. 
        self.testfacta = testfacta
        self.testfactb = testfactb
        self.testfactc = testfactc
        self.testfactd = testfactd
        self.testfacte = testfacte


        # default settings
        self.echo_level = 0
        self.compute_reactions = True
        self.ReformDofSetAtEachStep = True
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
        self.spalart_allmaras_linear_solver = None

        self.divergence_clearance_steps = 0

        print("Construction monolithic solver finished")

        print("after reading all the model contains:")
        print(self.model_part)

        if(self.eul_model_part == 0):
            self.UlfUtils = UlfUtils()
          #self.PfemUtils = PfemUtils()
            self.mark_outer_nodes_process = MarkOuterNodesProcess(model_part);
            self.node_erase_process = NodeEraseProcess(model_part);
            self.add_nodes=True
            self.alpha_shape = 3.5;
            #######self.ulf_apply_bc_process = UlfApplyBCProcess(model_part);
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
                self.Mesher =  TriGenDropletModeler()
                self.fluid_neigh_finder = FindNodalNeighboursProcess(model_part,9,18)
                #this is needed if we want to also store the conditions a node belongs to
                self.condition_neigh_finder = FindConditionsNeighboursProcess(model_part,2, 10)
                
            elif (domain_size == 3):
                #self.Mesher = TetGenDropletModeler()
                self.Mesher =TetGenPfemModeler()
                
                self.fluid_neigh_finder = FindNodalNeighboursProcess(model_part,20,30)
	      #this is needed if we want to also store the conditions a node belongs to
                self.condition_neigh_finder = FindConditionsNeighboursProcess(model_part,3, 20)

            (self.fluid_neigh_finder).Execute();
            Hfinder  = FindNodalHProcess(model_part);
            Hfinder.Execute();

    #
    def Initialize(self):

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

        # creating the solution strategy
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol, self.abs_vel_tol, self.rel_pres_tol, self.abs_pres_tol)
        #self.conv_criteria = ResidualCriteria(0.0001, 0.0000001)

        #(self.conv_criteria).SetEchoLevel(self.echo_level)

        self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(self.alpha, self.move_mesh_strategy, self.domain_size)

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

        self.model_part.ProcessInfo.SetValue(CONTACT_ANGLE_STATIC, self.contact_angle)
        self.model_part.ProcessInfo.SetValue(SURFACE_TENSION_COEF, self.gamma)
        
        #dissipative forces and testing related varialbes
        self.model_part.ProcessInfo.SetValue(DISSIPATIVE_FORCE_COEFF_JM_X, self.zeta_dissapative_JM_x)
        self.model_part.ProcessInfo.SetValue(DISSIPATIVE_FORCE_COEFF_BM_X, self.zeta_dissapative_BM_x)
        self.model_part.ProcessInfo.SetValue(DISSIPATIVE_FORCE_COEFF_SM_X, self.zeta_dissapative_SM_x)
        self.model_part.ProcessInfo.SetValue(DISSIPATIVE_FORCE_COEFF_JM_Y, self.zeta_dissapative_JM_y)
        self.model_part.ProcessInfo.SetValue(DISSIPATIVE_FORCE_COEFF_BM_Y, self.zeta_dissapative_BM_y)
        self.model_part.ProcessInfo.SetValue(DISSIPATIVE_FORCE_COEFF_SM_Y, self.zeta_dissapative_SM_y)
        self.model_part.ProcessInfo.SetValue(DISSIPATIVE_FORCE_COEFF_JM_Y, self.zeta_dissapative_JM_z)
        self.model_part.ProcessInfo.SetValue(DISSIPATIVE_FORCE_COEFF_BM_Y, self.zeta_dissapative_BM_z)
        self.model_part.ProcessInfo.SetValue(DISSIPATIVE_FORCE_COEFF_SM_Y, self.zeta_dissapative_SM_z)
        self.model_part.ProcessInfo.SetValue(TESTFACTA, self.testfacta)
        self.model_part.ProcessInfo.SetValue(TESTFACTB, self.testfactb)
        self.model_part.ProcessInfo.SetValue(TESTFACTC, self.testfactc)
        self.model_part.ProcessInfo.SetValue(TESTFACTD, self.testfactd)
        self.model_part.ProcessInfo.SetValue(TESTFACTE, self.testfacte)


        if(self.eul_model_part == 0):
	    #marking the fluid
            (self.fluid_neigh_finder).Execute();
            ######(self.ulf_apply_bc_process).Execute();
            if (self.domain_size == 2):
                FindTriplePoint().FindTriplePoint2D(self.model_part)
            self.Remesh()

# print "Initialization monolithic solver finished"
    #
    def Solve(self):

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
            FindTriplePoint().FindTriplePoint2D(self.model_part)
            CalculateCurvature().CalculateCurvature2D(self.model_part)
            CalculateNodalLength().CalculateNodalLength2D(self.model_part)
            CalculateContactAngle().CalculateContactAngle2D(self.model_part)
            self.cont_angle_cond()
        elif (self.domain_size == 3):
            self.UlfUtils.MarkLonelyNodesForErasing(self.model_part)
            FindTriplePoint().FindTriplePoint3D(self.model_part)
            CalculateCurvature().CalculateCurvature3D(self.model_part)
            #for node in self.model_part.Nodes:
                #node.SetSolutionStepValue(CONTACT_ANGLE,0,0.0)
            CalculateContactAngle().CalculateContactAngle3D(self.model_part)
            CalculateNodalLength().CalculateNodalLength3D(self.model_part)
            CalculateNormalEq().CalculateNormalEq3D(self.model_part)
	  #CalculateAdhesionForce().CalculateAdhesionForce3D(self.model_part)
            #self.cont_angle_cond3D()
            self.vel()
            #if (node.GetSolutionStepValue(TRIPLE_POINT) == 1.0 and node.GetSolutionStepValue(IS_STRUCTURE) == 1.0):
                #v_cl = Vector(3)
                #v_cl[0] = 0.0
                #v_cl[1] = 0.0
                #v_cl[2] = 0.0
                #v_cl[0] = node.GetSolutionStepValue(VELOCITY_X,0)
                #v_cl[1] = node.GetSolutionStepValue(VELOCITY_Y,0)
                #v_cl[2] = node.GetSolutionStepValue(VELOCITY_Z,0)
                #zeta_dissapative_JM_x = 1.0
                #zeta_dissapative_JM_y = 1.0
                #zeta_dissapative_JM_z = 1.0
                #if (v_cl[0] < 0.0):
                    #zeta_dissapative_JM_x = - 1.0 * zeta_dissapative_JM_x
                #if (v_cl[1] < 0.0):
                    #zeta_dissapative_JM_y = - 1.0 * zeta_dissapative_JM_y
                #if (v_cl[2] < 0.0):
                    #zeta_dissapative_JM_z = - 1.0 * zeta_dissapative_JM_z
                #mu  = 8.90 * 1e-4
                #gamma =  0.072
                ##capillary
                #v_clx_abs = abs (v_cl[0])
                #v_cly_abs = abs (v_cl[1])
                #v_clz_abs = abs (v_cl[2])
                #cap_x   =  mu *  v_clx_abs / gamma;
                #cap_y   =  mu *  v_cly_abs / gamma;
                #cap_z   =  mu *  v_clz_abs / gamma;
                #v_clx_2 = gamma * 2.24 * pow(cap_x, 0.54)
                #v_cly_2 = gamma * 2.24 * pow(cap_y, 0.54)
                #v_clz_2 = gamma * 2.24 * pow(cap_z, 0.54)
                #print("in the solve v_clx_2",  v_clx_2)
                #print("in the solve v_cly_2",  v_cly_2)
                #print("in the solve v_clz_2",  v_clz_2)
                ##using Jiang's Model : gamma tanh(4.96 Ca^(0.702))
                #node.SetSolutionStepValue(VELOCITY_X,0, v_clx_2)
                #node.SetSolutionStepValue(VELOCITY_Y,0, v_cly_2)
                #node.SetSolutionStepValue(VELOCITY_Z,0, v_clz_2)


        #self.solver.MoveMesh()
        (self.solver).Solve() #it dumps in this line... 20151020
        #AssignPointNeumannConditions().AssignPointNeumannConditions3D(self.model_part)

        #self.solver.MoveMesh()

        if(self.eul_model_part == 0):
            (self.fluid_neigh_finder).Execute();
            #(self.fluid_neigh_finder).Execute();
            self.Remesh();


    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)


    def Clear(self):
        (self.solver).Clear()



    ##########################################
    def Remesh(self):

        #self.UlfUtils.MarkNodesTouchingWall(self.model_part, self.domain_size, 0.08)

        ##erase all conditions and elements prior to remeshing
        ((self.model_part).Elements).clear();

        (self.mark_outer_nodes_process).MarkOuterNodes(self.box_corner1, self.box_corner2);

        #for node in self.model_part.Nodes:
            #if(node.GetSolutionStepValue(IS_BOUNDARY) > 0.9):
                #h_factor=0.5;
            #else:
                #h_factor=0.25

        h_factor=0.6

        if (self.domain_size == 2):
         (self.Mesher).ReGenerateMeshDROPLET("SurfaceTension2D","Condition2D", self.model_part, self.node_erase_process, True, True, self.alpha_shape, h_factor)
        elif (self.domain_size == 3):
         #(self.Mesher).ReGenerateMeshDROPLET3D("SurfaceTension3D","Condition3D", self.model_part, self.node_erase_process, True, False, self.alpha_shape, h_factor)      
         (self.Mesher).ReGenerateMesh("SurfaceTension3D", "Condition3D", self.model_part, self.node_erase_process, True, True, self.alpha_shape, h_factor)

         ##(self.Mesher).ReGenerateMesh("ASGSCompressible3D", "Condition3D", self.model_part,
                            ##(self.structure_model_part).Elements, self.node_erase_process, True, True, self.alpha_shape, self.h_factor)
         
        (self.fluid_neigh_finder).Execute();
        (self.condition_neigh_finder).Execute();

        #print "marking fluid" and applying fluid boundary conditions
        #######(self.ulf_apply_bc_process).Execute();
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
            FindTriplePoint().FindTriplePoint2D(self.model_part)
            CalculateCurvature().CalculateCurvature2D(self.model_part)
            CalculateNodalLength().CalculateNodalLength2D(self.model_part)
            CalculateContactAngle().CalculateContactAngle2D(self.model_part)
            self.cont_angle_cond()
            
        elif (self.domain_size == 3):
            for node in self.model_part.Nodes:
                if (node.GetSolutionStepValue(IS_FREE_SURFACE) > 0.999999999999):
                    node.SetSolutionStepValue(FLAG_VARIABLE, 0, 1.0)
                    node.SetSolutionStepValue(IS_INTERFACE, 0, 1.0)
                else:
                    node.SetSolutionStepValue(IS_INTERFACE, 0, 0.0)
            for node in self.model_part.Nodes:
                if (node.GetSolutionStepValue(IS_BOUNDARY) == 0.0):# and node.GetSolutionStepValue(TRIPLE_POINT) == 0):
                    node.SetSolutionStepValue(NORMAL_X,0,0.0)
                    node.SetSolutionStepValue(NORMAL_Y,0,0.0)
                    node.SetSolutionStepValue(NORMAL_Z,0,0.0)
            self.UlfUtils.MarkLonelyNodesForErasing(self.model_part)
            FindTriplePoint().FindTriplePoint3D(self.model_part)
            CalculateCurvature().CalculateCurvature3D(self.model_part)
            #for node in self.model_part.Nodes:
                #node.SetSolutionStepValue(CONTACT_ANGLE,0,0.0)
            CalculateContactAngle().CalculateContactAngle3D(self.model_part)
            CalculateNodalLength().CalculateNodalLength3D(self.model_part)
            CalculateNormalEq().CalculateNormalEq3D(self.model_part)
	  #CalculateAdhesionForce().CalculateAdhesionForce3D(self.model_part)
            #self.cont_angle_cond3D()
            self.vel()

        ##############THIS IS FOR EMBEDDED"""""""""""""""""""""""""
        #print("end of remesh function")
    ######################################################################

    def FindNeighbours(self):
        (self.neigh_finder).Execute();

    def cont_angle_cond(self):
        theta_adv = 105
        theta_rec = 70
	#theta_adv = self.contact_angle + 0.5
	#theta_rec = self.contact_angle - 0.5
        time = self.model_part.ProcessInfo.GetValue(TIME)
        dt = self.model_part.ProcessInfo.GetValue(DELTA_TIME)
	#x_mean = 0.0
	#found_tp = 0
	################## For sessile drop examples
        for node in self.model_part.Nodes:
            if (node.GetSolutionStepValue(TRIPLE_POINT) != 0.0):
                if ((node.GetSolutionStepValue(CONTACT_ANGLE) > theta_adv) or (node.GetSolutionStepValue(CONTACT_ANGLE) < theta_rec)):
                    node.Free(VELOCITY_X)
                else:
                    node.SetSolutionStepValue(VELOCITY_X,0, 0.0)
                    node.Fix(VELOCITY_X)
            if ((node.GetSolutionStepValue(TRIPLE_POINT) == 0.0) and (node.GetSolutionStepValue(IS_STRUCTURE) != 0.0)):
                    node.SetSolutionStepValue(VELOCITY_X,0, 0.0)
                    node.SetSolutionStepValue(VELOCITY_Y,0, 0.0)
                    node.Fix(VELOCITY_X)
                    node.Fix(VELOCITY_Y)
                    
                    
    def vel(self):
        for node in self.model_part.Nodes:
            theta_adv = self.contact_angle + 0.5
            theta_rec = self.contact_angle  - 0.5
            if (node.GetSolutionStepValue(IS_STRUCTURE) != 0.0):
                v_cl = Vector(3)
                v_cl[0] = node.GetSolutionStepValue(VELOCITY_X,0)
                v_cl[1] = node.GetSolutionStepValue(VELOCITY_Y,0)
                v_cl[2] = node.GetSolutionStepValue(VELOCITY_Z,0)
                v_cl[2] = 0.0
                node.SetSolutionStepValue(VELOCITY_X,0, v_cl[0])
                node.SetSolutionStepValue(VELOCITY_Y,0, v_cl[1])
                node.SetSolutionStepValue(VELOCITY_Z,0, v_cl[2])
                d_z = node.GetSolutionStepValue(VELOCITY_Z,0)
                d_z = 0.0
                node.SetSolutionStepValue(DISPLACEMENT_Z,0, d_z)
                node.Fix(DISPLACEMENT_Z)
                node.Free(VELOCITY_X);
                node.Free(VELOCITY_Y);
                if (node.GetSolutionStepValue(IS_BOUNDARY) != 1.0 and node.GetSolutionStepValue(TRIPLE_POINT)!= 0.0):
                    node.Set(TO_ERASE,True)
                if (node.GetSolutionStepValue(IS_BOUNDARY) != 1.0 and node.GetSolutionStepValue(IS_STRUCTURE)!= 0.0):
                    node.Set(TO_ERASE,True)
                    
                #december 18 2018
                #if (node.GetSolutionStepValue(TRIPLE_POINT) != 0.0):
                
                #######v_clx= node.GetSolutionStepValue(VELOCITY_X,0)
                #######v_cly= node.GetSolutionStepValue(VELOCITY_Y,0)
                ########v_clz= node.GetSolutionStepValue(VELOCITY_Z,0)
                #######d_clx= node.GetSolutionStepValue(DISPLACEMENT_X,0)
                #######d_cly= node.GetSolutionStepValue(DISPLACEMENT_X,0)
                ########d_clz= node.GetSolutionStepValue(VELOCITY_Z,0)
                    
                #######constx = 1.0
                #######consty = 1.0
                
                #######if (v_clx < 0):
                    #######constx = - 1.0
                #######elif (v_clx > 0):
                    #######constx = 1.0
                
                #######if (v_cly < 0):
                    #######consty = - 1.0
                #######elif (v_cly > 0):
                    #######consty = 1.0
                    
                #######v_clxabs = abs(v_clx)
                #######v_clyabs = abs(v_cly)
                    
                #######if (v_clxabs > v_clyabs):
                        #######v_clyabs = v_clxabs
                #######else:
                        #######v_clxabs = v_clyabs
                    
                #######v_clxabs = constx * v_clxabs
                #######v_clyabs = consty * v_clyabs
                
                #######node.Free(DISPLACEMENT_X)
                #######node.Free(DISPLACEMENT_Y)
                    
                    
                    #if ((node.GetSolutionStepValue(CONTACT_ANGLE) < theta_adv) and (node.GetSolutionStepValue(CONTACT_ANGLE) > theta_rec)):
                        #ca=node.GetSolutionStepValue(CONTACT_ANGLE,0)
                        #ca = ca
                        #node.SetSolutionStepValue(CONTACT_ANGLE,0,ca)
                        #node.Fix(CONTACT_ANGLE);
                        #node.Fix(VELOCITY_X);
                        #node.Fix(VELOCITY_Y);
                        #node.Fix(VELOCITY_Z);
                        #node.Fix(DISPLACEMENT_X);
                        #node.Fix(DISPLACEMENT_Y);
                        #node.Fix(DISPLACEMENT_Z);
                        #node.Fix(CONTACT_ANGLE);
                
            #if (node.GetSolutionStepValue(TRIPLE_POINT) != 0.0):
              #if ((node.GetSolutionStepValue(CONTACT_ANGLE) < theta_adv) and (node.GetSolutionStepValue(CONTACT_ANGLE) > theta_rec)):
                #v_cl = Vector(3)
                #v_cl[0] = 0.0
                #v_cl[1] = 0.0
                #v_cl[2] = 0.0
                #v_cl[0] = node.GetSolutionStepValue(VELOCITY_X,0)
                #v_cl[1] = node.GetSolutionStepValue(VELOCITY_Y,0)
                #v_cl[2] = node.GetSolutionStepValue(VELOCITY_Z,0)
                #d_cl = Vector(3)
                #d_cl[0] = 0.0
                #d_cl[1] = 0.0
                #d_cl[2] = 0.0
                #d_cl[0] = node.GetSolutionStepValue(DISPLACEMENT_X,0)
                #d_cl[1] = node.GetSolutionStepValue(DISPLACEMENT_Y,0)
                #d_cl[2] = node.GetSolutionStepValue(DISPLACEMENT_Z,0)
                #d_cl[0] = 0.0
                #d_cl[1] = 0.0
                #d_cl[2] = 0.0
                #v_cl[0] = 0.0
                #v_cl[1] = 0.0
                #v_cl[2] = 0.0
                ##node.SetSolutionStepValue(VELOCITY_X,0, v_cl[0])
                ##node.SetSolutionStepValue(VELOCITY_Y,0, v_cl[1])
                #node.SetSolutionStepValue(VELOCITY_Z,0, v_cl[2])
                ##node.SetSolutionStepValue(DISPLACEMENT_X,0, d_cl[0])
                ##node.SetSolutionStepValue(DISPLACEMENT_Y,0, d_cl[1])
                #node.SetSolutionStepValue(DISPLACEMENT_Z,0, d_cl[2])
              #if ((node.GetSolutionStepValue(CONTACT_ANGLE) < theta_adv) or (node.GetSolutionStepValue(CONTACT_ANGLE) > theta_rec)):
                #(1) case, no dissipative force is used
                #v_cl = Vector(3)
                #v_cl[0] = 0.0
                #v_cl[1] = 0.0
                #v_cl[2] = 0.0
                #v_cl[0] = node.GetSolutionStepValue(VELOCITY_X,0)
                #v_cl[1] = node.GetSolutionStepValue(VELOCITY_Y,0)
                #v_cl[2] = node.GetSolutionStepValue(VELOCITY_Z,0)
                #v_cl[2] = 0.0
                #node.SetSolutionStepValue(VELOCITY_X,0, v_cl[0])
                #node.SetSolutionStepValue(VELOCITY_Y,0, v_cl[1])
                #node.SetSolutionStepValue(VELOCITY_Z,0, v_cl[2])
                #d_z = node.GetSolutionStepValue(VELOCITY_X,0)
                #d_z = 0.0
                #node.SetSolutionStepValue(DISPLACEMENT_Z,0, d_z)
                #node.Fix(DISPLACEMENT_Z)
                #node.Free(VELOCITY_X);
                #node.Free(VELOCITY_Y);
                #node.Free(VELOCITY_Z);
                
                #case (2) first model
                #v_cl = Vector(3)
                #v_cl[0] = 0.0
                #v_cl[1] = 0.0
                #v_cl[2] = 0.0
                #v_cl[0] = node.GetSolutionStepValue(VELOCITY_X,0)
                #v_cl[1] = node.GetSolutionStepValue(VELOCITY_Y,0)
                #v_cl[2] = node.GetSolutionStepValue(VELOCITY_Z,0)
                ##using BM's Model; here we will use zeta_dissapative_JM_x and zeta_dissapative_JM_y and zeta_dissapative_JM_z to help us in getting the right direction of the dissipative foce as:
                ##zeta_dissapative_JM_x1 = 1.0
                ##zeta_dissapative_JM_y1 = 1.0
                ##zeta_dissapative_JM_z1 = 1.0
                #if (v_cl[0] < 0.0):
                    #zeta_dissapative_JM_x1 = - 1.0
                #elif (v_cl[0] >= 0.0):
                    #zeta_dissapative_JM_x1 =  1.0 
                    
                #if (v_cl[1] < 0.0):
                    #zeta_dissapative_JM_y1 = - 1.0 
                #elif (v_cl[1] >= 0.0):
                    #zeta_dissapative_JM_y1 = 1.0
                    
                ##if (v_cl[2] < 0.0):
                    ##zeta_dissapative_JM_z1 = - 1.0 
                ##elif (v_cl[2] >= 0.0):
                    ##zeta_dissapative_JM_z1 = 1.0 
                #mu  = 8.90 * 1e-4
                #gamma =  0.072
                ##capillary
                #v_clx_abs = abs (v_cl[0])
                #v_cly_abs = abs (v_cl[1])
                ##v_clz_abs = abs (v_cl[2])
                #cap_x   =  mu *  v_clx_abs / gamma;
                #cap_y   =  mu *  v_cly_abs / gamma;
                ##cap_z   =  mu *  v_clz_abs / gamma;
                ##we will multiply by (node.GetSolutionStepValue(NODAL_LENGTH,0)/2.0) below because this the integral (on the contact line)
                #v_cl[0] = gamma * 2.24 * pow(cap_x, 0.54) * zeta_dissapative_JM_x1 * node.GetSolutionStepValue(NODAL_LENGTH,0)/2.0
                #v_cl[1] = gamma * 2.24 * pow(cap_y, 0.54) * zeta_dissapative_JM_y1 * node.GetSolutionStepValue(NODAL_LENGTH,0)/2.0
                ##v_cl[2] = gamma * 2.24 * pow(cap_z, 0.54) * zeta_dissapative_JM_z1 * node.GetSolutionStepValue(NODAL_LENGTH,0)/2.0
                ##force before adding the dissipative force (f_cl is the contact line force)
                #f_cl = Vector(3)
                #f_cl[0] = 0.0
                #f_cl[1] = 0.0
                #f_cl[2] = 0.0
                #f_cl[0] = node.GetSolutionStepValue(FORCE_X,0)
                #f_cl[1] = node.GetSolutionStepValue(FORCE_Y,0)
                ##f_cl[2] = node.GetSolutionStepValue(FORCE_Z,0)
                ##sof_cl is summation of force (before adding the dissipative force + dissipative force (f_cl))
                #sof_cl = Vector(3)
                #sof_cl[0] = f_cl[0] - v_cl[0]
                #sof_cl[1] = f_cl[1] - v_cl[1]
                ##sof_cl[2] = f_cl[2] - v_cl[2]
                #node.SetSolutionStepValue(FORCE_X,0, sof_cl[0])
                #node.SetSolutionStepValue(FORCE_Y,0, sof_cl[1])
                ##node.SetSolutionStepValue(FORCE_Z,0, sof_cl[2])
                #time = self.model_part.ProcessInfo.GetValue(TIME)
                #dt   = self.model_part.ProcessInfo.GetValue(DELTA_TIME)
                #CalculateNodalLength().CalculateNodalLength3D(self.model_part)
                #h = node.GetSolutionStepValue(NODAL_LENGTH,0)
                #v_x = 0.0
                #v_y = 0.0
                ##v_z = 0.0
                ##velocity = (F / mass_line) * time ; where the mass_line = NODAL_LENGTH/2 * density (1000.0 for water)
                #v_x =  v_cl[0] * time / (h/2.0 * 1000.0) 
                #v_y =  v_cl[1] * time / (h/2.0 * 1000.0) 
                ##v_z =  v_cl[2] * time / (h/2.0 * 1000.0)
                ##now we add the dissipative velocity to the current time step velocity
                #v_x_2 = node.GetSolutionStepValue(VELOCITY_X,0) - v_x
                #v_y_2 = node.GetSolutionStepValue(VELOCITY_Y,0) - v_y
                ##v_z_2 = node.GetSolutionStepValue(VELOCITY_Z,0) - v_z
                #node.SetSolutionStepValue(VELOCITY_X,0, (v_x_2))
                #node.SetSolutionStepValue(VELOCITY_Y,0, (v_y_2))
                ##node.SetSolutionStepValue(VELOCITY_Z,0, (v_z_2))
                ##mow we add the velocity to the DISPLACEMENT
                #node.SetSolutionStepValue(DISPLACEMENT_X,0, (node.GetSolutionStepValue(DISPLACEMENT_X,0) - v_x * time))
                #node.SetSolutionStepValue(DISPLACEMENT_Y,0, (node.GetSolutionStepValue(DISPLACEMENT_Y,0) - v_y * time))
                ##here we fix the DISPLACEMENT_Z and VELOCITY_Z zero
                #disp_z=0.0
                #node.SetSolutionStepValue(DISPLACEMENT_Z,0, disp_z)
                #disp_z=0.0
                #node.Fix(DISPLACEMENT_Z)
                #node.SetSolutionStepValue(DISPLACEMENT_Z,0, disp_z)
                #d_z = node.GetSolutionStepValue(VELOCITY_Z,0)
                #d_z = 0.0
                #node.SetSolutionStepValue(VELOCITY_Z,0, d_z)
                #node.Fix(VELOCITY_Z)
                ##finally we free the X and Y velocity to move as coded with the dissipative forces
                #node.Free(DISPLACEMENT_X)
                #node.Free(VELOCITY_X)
                #node.Free(DISPLACEMENT_Y)
                #node.Free(VELOCITY_Y)
        
        
    def cont_angle_cond3D(self):
        theta_adv = self.contact_angle 
        theta_rec = self.contact_angle 
        time = self.model_part.ProcessInfo.GetValue(TIME)
        dt = self.model_part.ProcessInfo.GetValue(DELTA_TIME)
	################## For sessile drop examples
        for node in self.model_part.Nodes:
                #if (node.GetSolutionStepValue(IS_STRUCTURE) != 0.0 and node.GetSolutionStepValue(TRIPLE_POINT) != 1.0):
                    ###node.Set(TO_ERASE,False)
                    #node.SetSolutionStepValue(VELOCITY_X,0, 0.0)
                    #node.SetSolutionStepValue(VELOCITY_Y,0, 0.0)
                    #node.SetSolutionStepValue(VELOCITY_Z,0, 0.0)
                    #node.Fix(VELOCITY_X)
                    #node.Fix(VELOCITY_Y)
                    #node.Fix(VELOCITY_Z)
                #if (node.GetSolutionStepValue(TRIPLE_POINT) != 0.0):
                        #if ((node.GetSolutionStepValue(CONTACT_ANGLE) < theta_adv) and (node.GetSolutionStepValue(CONTACT_ANGLE) > theta_rec)):
                            #v_cl = Vector(3)
                            #v_cl[0] = 0.0
                            #v_cl[1] = 0.0
                            #v_cl[2] = 0.0
                            #v_cl[0] = node.GetSolutionStepValue(VELOCITY_X,0)
                            #v_cl[1] = node.GetSolutionStepValue(VELOCITY_Y,0)
                            #v_cl[2] = node.GetSolutionStepValue(VELOCITY_Z,0)
                            #d_cl = Vector(3)
                            #d_cl[0] = 0.0
                            #d_cl[1] = 0.0
                            #d_cl[2] = 0.0
                            #d_cl[0] = node.GetSolutionStepValue(DISPLACEMENT_X,0)
                            #d_cl[1] = node.GetSolutionStepValue(DISPLACEMENT_Y,0)
                            #d_cl[2] = node.GetSolutionStepValue(DISPLACEMENT_Z,0)
                            #d_cl[0] = 0.0
                            #d_cl[1] = 0.0
                            #d_cl[2] = 0.0
                            #v_cl[0] = 0.0
                            #v_cl[1] = 0.0
                            #v_cl[2] = 0.0
                            #node.SetSolutionStepValue(VELOCITY_X,0, v_cl[0])
                            #node.SetSolutionStepValue(VELOCITY_Y,0, v_cl[1])
                            #node.SetSolutionStepValue(VELOCITY_Z,0, v_cl[2])
                            #node.SetSolutionStepValue(DISPLACEMENT_X,0, d_cl[0])
                            #node.SetSolutionStepValue(DISPLACEMENT_Y,0, d_cl[1])
                            #node.SetSolutionStepValue(DISPLACEMENT_Z,0, d_cl[2])
                            
                        #else:
                        if ((node.GetSolutionStepValue(CONTACT_ANGLE) < theta_adv) or (node.GetSolutionStepValue(CONTACT_ANGLE) > theta_rec)):
                                node.Free(VELOCITY_X);
                                node.Free(VELOCITY_Y);
                                node.Free(VELOCITY_Z);
                                if (node.Z < - 1e-5):
                                    node.Set(TO_ERASE,True)
                                if (node.GetSolutionStepValue(IS_BOUNDARY) != 1.0 and node.GetSolutionStepValue(IS_STRUCTURE)!= 0.0):
                                    node.Set(TO_ERASE,True)
                        if (node.GetSolutionStepValue(IS_BOUNDARY) != 1.0 and node.GetSolutionStepValue(TRIPLE_POINT)!= 0.0):
                            node.Set(TO_ERASE,True)
                        if (node.GetSolutionStepValue(IS_BOUNDARY) != 1.0 and node.GetSolutionStepValue(IS_STRUCTURE)!= 0.0):
                            node.Set(TO_ERASE,True)


    
def CreateSolver(model_part, config, eul_model_part, gamma, contact_angle, zeta_dissapative_JM_x, zeta_dissapative_BM_x, zeta_dissapative_SM_x, zeta_dissapative_JM_y, zeta_dissapative_BM_y, zeta_dissapative_SM_y, zeta_dissapative_JM_z, zeta_dissapative_BM_z, zeta_dissapative_SM_z, testfacta, testfactb, testfactc, testfactd, testfacte): #FOR 3D!
    fluid_solver = STMonolithicSolver(model_part, config.domain_size, eul_model_part, gamma, contact_angle, zeta_dissapative_JM_x, zeta_dissapative_BM_x, zeta_dissapative_SM_x, zeta_dissapative_JM_y, zeta_dissapative_BM_y, zeta_dissapative_SM_y, zeta_dissapative_JM_z, zeta_dissapative_BM_z, zeta_dissapative_SM_z, testfacta, testfactb, testfactc, testfactd, testfacte)

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
        
    pDiagPrecond = DiagonalPreconditioner()
    #self.linear_solver = BICGSTABSolver(1e-6, 5000, pDiagPrecond)

    #import linear_solver_factory
    #if(hasattr(config, "linear_solver_config")):
        #fluid_solver.linear_solver = linear_solver_factory.ConstructSolver(
            #config.linear_solver_config)

    return fluid_solver 
