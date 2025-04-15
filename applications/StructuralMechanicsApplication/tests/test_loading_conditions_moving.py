import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import math

class TestLoadingConditionsMoving(KratosUnittest.TestCase):

    @staticmethod
    def setup_strategy(mp):
        """
        Setup default strategy for solving the problem

        Parameters
        ----------
        mp: model part

        Returns default strategy
        -------

        """

        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-1, 1e-1)
        convergence_criterion.SetEchoLevel(0)

        max_iters = 10
        compute_reactions = False
        reform_step_dofs = False
        move_mesh_flag = False
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(mp,
                                                                         scheme,
                                                                         convergence_criterion,
                                                                         builder_and_solver,
                                                                         max_iters,
                                                                         compute_reactions,
                                                                         reform_step_dofs,
                                                                         move_mesh_flag)

        strategy.SetEchoLevel(0)
        return strategy

    @staticmethod
    def calculate_rotation_of_beam(disp_left, disp_right, rot_left, rot_right, length_beam, x, rotation):
        """
        Calculate rotation at a certain position along an euler bernoulli beam

        Parameters
        ----------
        disp_left: left displacement of beam [m]
        disp_right: right displacement of beam [m]
        rot_left: left rotation of beam [rad]
        rot_right: right rotation of beam [rad]
        length_beam: length of beam [m]
        x: local position of the load within the beam [m]

        Returns expected rotation at location of load
        -------

        """
        eta_left = disp_left * math.cos(rotation)
        eta_right = disp_right * math.cos(rotation)

        n1 = 6 * (x**2 / length_beam**3) - 6 * (x / length_beam**2)
        n2 = 1 + 3*x ** 2 / length_beam ** 2 - 4 * x / length_beam
        n3 = -6 * (x**2 / length_beam**3) + 6 * (x / length_beam**2)
        n4 = 3 * x ** 2 / length_beam ** 2 - 2 * x / length_beam

        expected_rotation = n1 * eta_left + n2 * rot_left + n3 * eta_right + n4 * rot_right

        return expected_rotation

    @staticmethod
    def calculate_deflection_of_beam(disp_left, disp_right, rot_left, rot_right, length_beam, x, rotation):
        """
        Calculate vertical deflection at a certain position along an Euler Bernoulli beam

        Parameters
        ----------
        disp_left: left displacement of beam [m]
        disp_right: right displacement of beam [m]
        rot_left: left rotation of beam [rad]
        rot_right: right rotation of beam [rad]
        length_beam: length of beam [m]
        x: local position of the load within the beam [m]

        Returns expected rotation at location of load
        -------

        """
        xi_left, eta_left = disp_left * math.sin(rotation), disp_left * math.cos(rotation)
        xi_right, eta_right = disp_right * math.sin(rotation), disp_right * math.cos(rotation)

        n1 = 1 + 2*(x/length_beam)**3 - 3*(x/length_beam)**2
        n2 = x + x**3/length_beam**2 - 2*x**2/length_beam
        n3 = -2*(x/length_beam)**3 + 3*(x/length_beam)**2
        n4 = x**3/length_beam**2 - x**2/length_beam

        expected_deflection_eta = n1*eta_left + n2*rot_left + n3*eta_right + n4*rot_right
        expected_deflection_xi = (1 - x / length_beam) * xi_left + (x / length_beam) * xi_right

        expected_deflection_x = expected_deflection_xi*math.cos(rotation) - expected_deflection_eta*math.sin(rotation)
        expected_deflection_y = expected_deflection_xi*math.sin(rotation) + expected_deflection_eta*math.cos(rotation)
        return expected_deflection_x, expected_deflection_y

    def __CalculateReaction2n(self, length, angle, load, distance):
        """
        Calculate reaction forces on 2 noded 2d inclined element

        Parameters
        ----------
        length: length of the element
        angle: angle of the element
        load: load vector
        distance: local distance of the load within the element

        Returns reaction forces
        -------

        """

        cos = math.cos(angle)
        sin = math.sin(angle)

        p_shear = load[1] * sin + load[0] * -cos
        a = distance
        b = length - distance

        left_shear_reaction = p_shear * b/length
        right_shear_reaction = p_shear * a/length

        x_l = left_shear_reaction * -sin
        y_l = left_shear_reaction * cos

        x_r = right_shear_reaction * -sin
        y_r = right_shear_reaction * cos

        return (x_l, y_l), (x_r, y_r)


    def __CalculateReaction3n(self, length, angle, load, distance):
        """
        Calculate reaction forces on 3 noded 2d inclined element

        Parameters
        ----------
        length: length of the element
        angle: angle of the element
        load: load vector
        distance: local distance of the load within the element

        Returns reaction forces
        -------

        """

        cos = math.cos(angle)
        sin = math.sin(angle)

        # calculate local forces
        p_shear = load[1] * sin + load[0] * -cos
        p_normal = load[1] * cos + load[0] * sin

        # calculate local reaction forces
        left_shear_reaction = p_shear * (length**2 - 3*length*distance + 2*distance**2)/length**2
        right_shear_reaction = p_shear * (2*distance**2 - length*distance)/length**2
        mid_shear_reaction = p_shear * 4 * (length * distance - distance ** 2) / length ** 2

        left_normal_reaction = p_normal * (length**2 - 3*length*distance + 2*distance**2)/length**2
        right_normal_reaction = p_normal * (2*distance**2 - length*distance)/length**2
        mid_normal_reaction = p_normal * 4 * (length * distance - distance ** 2) / length ** 2

        # calculate global force
        x_l = left_shear_reaction * -sin + left_normal_reaction * cos
        y_l = left_shear_reaction * cos + left_normal_reaction * sin

        x_mid = mid_shear_reaction * -sin + mid_normal_reaction * cos
        y_mid = mid_shear_reaction * cos + mid_normal_reaction * sin

        x_r = right_shear_reaction * -sin + right_normal_reaction * cos
        y_r = right_shear_reaction * cos + right_normal_reaction * sin

        return (x_l, y_l), (x_mid, y_mid), (x_r, y_r)

    def __CalculateReactionWithRotation(self, length, angle, load, distance):
        """
        Calculate reaction forces including moment on a 2d inclined element

        Parameters
        ----------
        length: length of the element
        angle: angle of the element
        load: load vector
        distance: local distance of the load within the element

        Returns reaction forces
        -------

        """

        cos = math.cos(angle)
        sin = math.sin(angle)

        p_shear = load[1] * sin + load[0] * -cos
        a = distance
        b = length - distance

        expected_left_moment = p_shear * a * b ** 2 / length ** 2
        expected_right_moment = -p_shear * a ** 2 * b / length ** 2
        moment_left, moment_right = expected_left_moment, expected_right_moment

        left_shear_reaction = (p_shear * b) / length + (moment_left + moment_right) / length
        right_shear_reaction = (p_shear * a) / length - (moment_left + moment_right) / length

        x_l = left_shear_reaction * -sin
        y_l = left_shear_reaction * cos

        x_r = right_shear_reaction * -sin
        y_r = right_shear_reaction * cos

        return (x_l, y_l, moment_left), (x_r, y_r, moment_right)

    def _MovingLoadCondition3D2NRotDofZ(self):
        """"
        Check moving load on a 3d 2n line element, where the line element is inclined around z axis.
        """

        current_model = KratosMultiphysics.Model()

        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

        strategy = self.setup_strategy(mp)

        # create nodes
        second_coord = [math.sqrt(2), math.sqrt(2.0), 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0], second_coord[1], second_coord[2])

        length = 2.0
        rotation = math.atan(second_coord[1] / second_coord[0])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, mp)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z, mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition3D2N", 1, [1, 2], mp.GetProperties()[1])

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set perpendicular POINT_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 2.0
        load_on_cond[1] = -2.0
        load_on_cond[2] = 0.0

        cond.SetValue(StructuralMechanicsApplication.POINT_LOAD, load_on_cond)

        # compare results when load is located at different points
        num_tests = 5
        for i in range(num_tests):
            distance = i/(num_tests-1) * length
            (x_left, y_left, moment_left), (x_right, y_right, moment_right) = \
                self.__CalculateReactionWithRotation(length, rotation, load_on_cond, distance)

            # set load on quarter
            cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, distance)
            strategy.InitializeSolutionStep()
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

            self.assertAlmostEqual(rhs[0], x_left)
            self.assertAlmostEqual(rhs[1], y_left)
            self.assertAlmostEqual(rhs[2], 0)
            self.assertAlmostEqual(rhs[3], 0)
            self.assertAlmostEqual(rhs[4], 0)
            self.assertAlmostEqual(rhs[5], moment_left)
            self.assertAlmostEqual(rhs[6], x_right)
            self.assertAlmostEqual(rhs[7], y_right)
            self.assertAlmostEqual(rhs[8], 0)
            self.assertAlmostEqual(rhs[9], 0)
            self.assertAlmostEqual(rhs[10], 0)
            self.assertAlmostEqual(rhs[11], moment_right)

        # set load outside condition
        cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, length * 2)
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], 0)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)
        self.assertAlmostEqual(rhs[4], 0)
        self.assertAlmostEqual(rhs[5], 0)
        self.assertAlmostEqual(rhs[6], 0)
        self.assertAlmostEqual(rhs[7], 0)
        self.assertAlmostEqual(rhs[8], 0)
        self.assertAlmostEqual(rhs[9], 0)
        self.assertAlmostEqual(rhs[10], 0)
        self.assertAlmostEqual(rhs[11], 0)

    def _MovingLoadCondition3D2NRotDofY(self):
        """"
        Check moving load on a 3d 2n line element, where the line element is inclined around y axis.
        """

        current_model = KratosMultiphysics.Model()

        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

        strategy = self.setup_strategy(mp)

        # create nodes
        second_coord = [math.sqrt(2),  0.0, math.sqrt(2.0)]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, second_coord[0], second_coord[1], second_coord[2])

        length = 2.0
        rotation = math.atan(second_coord[2] / second_coord[0])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, mp)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z, mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition3D2N", 1, [1, 2], mp.GetProperties()[1])

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set perpendicular POINT_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 2.0
        load_on_cond[1] = 0.0
        load_on_cond[2] = -2.0

        cond.SetValue(StructuralMechanicsApplication.POINT_LOAD, load_on_cond)

        # compare results when load is located at different points
        num_tests = 5
        for i in range(num_tests):
            distance = i/(num_tests-1) * length
            (x_left, z_left, moment_left), (x_right, z_right, moment_right) = \
                self.__CalculateReactionWithRotation(length, rotation, [load_on_cond[0], load_on_cond[2],
                                                                      load_on_cond[1]], distance)

            # set load on quarter
            cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, distance)

            strategy.InitializeSolutionStep()
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

            self.assertAlmostEqual(rhs[0], x_left)
            self.assertAlmostEqual(rhs[1], 0)
            self.assertAlmostEqual(rhs[2], z_left)
            self.assertAlmostEqual(rhs[3], 0)
            self.assertAlmostEqual(rhs[4], moment_left)
            self.assertAlmostEqual(rhs[5], 0)
            self.assertAlmostEqual(rhs[6], x_right)
            self.assertAlmostEqual(rhs[7], 0)
            self.assertAlmostEqual(rhs[8], z_right)
            self.assertAlmostEqual(rhs[9], 0)
            self.assertAlmostEqual(rhs[10], moment_right)
            self.assertAlmostEqual(rhs[11], 0)

    def _MovingLoadCondition2D2N(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        strategy = self.setup_strategy(mp)

        # create nodes
        second_coord = [math.sqrt(2), math.sqrt(2.0), 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],second_coord[2])
        length = 2.0
        rotation = math.atan(second_coord[1] / second_coord[0])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1,2], mp.GetProperties()[1])

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # set perpendicular POINT_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 2.0
        load_on_cond[1] = -2.0
        load_on_cond[2] = 0.0 # note that this is a 2D condition

        num_tests = 5
        for i in range(num_tests):
            distance = i / (num_tests - 1) * length

            (x_left, y_left), (x_right, y_right) = self.__CalculateReaction2n(length, rotation, load_on_cond, distance)

            cond.SetValue(StructuralMechanicsApplication.POINT_LOAD,load_on_cond)

            # set load on node
            cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, distance)
            strategy.InitializeSolutionStep()
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

            self.assertAlmostEqual(rhs[0], x_left)
            self.assertAlmostEqual(rhs[1], y_left)
            self.assertAlmostEqual(rhs[2], x_right)
            self.assertAlmostEqual(rhs[3], y_right)

        # set load outside condition
        cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, length * 2)
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], 0)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)

    def _MovingLoadCondition2D2NRot(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

        strategy = self.setup_strategy(mp)
        # create nodes
        second_coord = [math.sqrt(2), math.sqrt(2.0), 0.0]
        mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        mp.CreateNewNode(2, math.sqrt(2), math.sqrt(2.0), 0.0)
        length = 2.0
        rotation = math.atan(second_coord[1] / second_coord[0])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,
                                                  mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1, 2], mp.GetProperties()[1])

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0, 0)
        rhs = KratosMultiphysics.Vector(0)

        # set perpendicular POINT_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 2.0
        load_on_cond[1] = -2.0
        load_on_cond[2] = 0.0  # note that this is a 2D condition

        cond.SetValue(StructuralMechanicsApplication.POINT_LOAD, load_on_cond)

        num_tests = 5
        for i in range(num_tests):
            distance = i / (num_tests - 1) * length
            (x_left, y_left, moment_left), (x_right, y_right, moment_right) = \
                self.__CalculateReactionWithRotation(length, rotation, [load_on_cond[0], load_on_cond[1],
                                                                        load_on_cond[2]], distance)

            # set load on node
            cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, distance)
            strategy.InitializeSolutionStep()
            strategy.Solve()
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

            self.assertAlmostEqual(rhs[0], x_left)
            self.assertAlmostEqual(rhs[1], y_left)
            self.assertAlmostEqual(rhs[2], moment_left)
            self.assertAlmostEqual(rhs[3], x_right)
            self.assertAlmostEqual(rhs[4], y_right)
            self.assertAlmostEqual(rhs[5], moment_right)

    def _MovingLoadCondition2D2NRot_disp_rot_at_load(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

        strategy = self.setup_strategy(mp)

        # create nodes
        coords_node_2 = [math.sqrt(2.0), math.sqrt(2.0), 0.0]

        node_1 = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        node_2 = mp.CreateNewNode(2, coords_node_2[0], coords_node_2[1], coords_node_2[2])

        rotation = math.atan(coords_node_2[1] / coords_node_2[0])
        length = math.sqrt(coords_node_2[0]**2 + coords_node_2[1]**2)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,
                                                  mp)

        disp_left, disp_right = 2, 0.0
        node_1.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, disp_left)
        node_2.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, disp_right)
        node_1.Fix(KratosMultiphysics.DISPLACEMENT_Y)
        node_2.Fix(KratosMultiphysics.DISPLACEMENT_Y)

        location_load = 0.5

        expected_deflection_x, expected_deflection_y = self.calculate_deflection_of_beam(disp_left, disp_right,
                                                                                         0, 0, length,
                                                                                         location_load, rotation)

        expected_rotation = self.calculate_rotation_of_beam(disp_left, disp_right, 0, 0, length,
                                                            location_load, rotation)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1, 2], mp.GetProperties()[1])

        cond.SetValue(StructuralMechanicsApplication.POINT_LOAD, [0, 1, 0])
        cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, location_load)

        # solve system and calculate displacement and rotation at the location of the load
        strategy.InitializeSolutionStep()
        strategy.SolveSolutionStep()

        # get calculated displacement at the load
        disp = cond.GetValue(KratosMultiphysics.DISPLACEMENT)

        self.assertAlmostEqual(disp[0], expected_deflection_x)
        self.assertAlmostEqual(disp[1], expected_deflection_y)

        # get calculated rotation at the load
        rot = cond.GetValue(KratosMultiphysics.ROTATION)

        self.assertAlmostEqual(rot[2], expected_rotation)

    def _MovingLoadCondition2D2NImposedRot_disp_rot_at_load(self):
        """
        Tests the MovingLoadCondition2D2N with imposed rotation and checks the displacement and rotation at the position
        of the load
        """

        vertical_point_load = 1.0
        EI = 300

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

        strategy = self.setup_strategy(mp)

        # create nodes
        coords_node_2 = [5, 0.0, 0.0]

        node_1 = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        node_2 = mp.CreateNewNode(2, coords_node_2[0], coords_node_2[1], coords_node_2[2])

        rotation = math.atan(coords_node_2[1] / coords_node_2[0])
        length = math.sqrt(coords_node_2[0]**2 + coords_node_2[1]**2)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,
                                                  mp)

        disp_left, disp_right = 0.0, 0.0
        node_1.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, disp_left)
        node_2.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, disp_right)
        node_1.Fix(KratosMultiphysics.DISPLACEMENT_Y)
        node_2.Fix(KratosMultiphysics.DISPLACEMENT_Y)

        # set imposed rotation based on euler beam theory with a point load in the middle of a simply supported beam
        rot_left = vertical_point_load * length**2 / (16*EI)
        rot_right = -vertical_point_load * length**2 / (16*EI)

        node_1.SetSolutionStepValue(KratosMultiphysics.ROTATION_Z, 0, rot_left)
        node_2.SetSolutionStepValue(KratosMultiphysics.ROTATION_Z, 0, rot_right)
        node_1.Fix(KratosMultiphysics.ROTATION_Z)
        node_2.Fix(KratosMultiphysics.ROTATION_Z)

        location_load = length/2

        # calculate analytical deflection and rotation at the location of the load
        expected_deflection_x, expected_deflection_y = self.calculate_deflection_of_beam(disp_left, disp_right,
                                                                                         rot_left, rot_right, length,
                                                                                         location_load, rotation)

        expected_rotation = self.calculate_rotation_of_beam(disp_left, disp_right, rot_left, rot_right, length,
                                                            location_load, rotation)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D2N", 1, [1, 2], mp.GetProperties()[1])

        cond.SetValue(StructuralMechanicsApplication.POINT_LOAD, [0, vertical_point_load, 0])
        cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, location_load)

        # solve system and calculate displacement and rotation at the location of the load
        strategy.InitializeSolutionStep()
        strategy.SolveSolutionStep()

        # get calculated displacement at the load
        disp = cond.GetValue(KratosMultiphysics.DISPLACEMENT)

        self.assertAlmostEqual(disp[0], expected_deflection_x)
        self.assertAlmostEqual(disp[1], expected_deflection_y)

        # get calculated rotation at the load
        rot = cond.GetValue(KratosMultiphysics.ROTATION)

        self.assertAlmostEqual(rot[2], expected_rotation)

    def _MovingLoadCondition3D2NRot_disp_rot_at_load(self, second_coordinates, axis_1, axis_2):

        out_of_plane_axis = list({0, 1, 2} - {axis_1, axis_2})[0]

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

        strategy = self.setup_strategy(mp)

        # create nodes
        node_1 = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        node_2 = mp.CreateNewNode(2, second_coordinates[0], second_coordinates[1], second_coordinates[2])

        length = math.sqrt(second_coordinates[0] ** 2 + second_coordinates[1] ** 2 + second_coordinates[2] ** 2)

        rotation = math.atan(second_coordinates[axis_2] / second_coordinates[axis_1])

        # add dofs
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, mp)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,
                                                  mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,
                                                  mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,
                                                  mp)

        # set displacement
        disp_left, disp_right = 2, 0.0

        if axis_2 == 0:
            var_disp = KratosMultiphysics.DISPLACEMENT_X
        elif axis_2 == 1:
            var_disp = KratosMultiphysics.DISPLACEMENT_Y
        elif axis_2 == 2:
            var_disp = KratosMultiphysics.DISPLACEMENT_Z
        else:
            raise Exception("Wrong axis_2, must be 0, 1 or 2")

        node_1.SetSolutionStepValue(var_disp, 0, disp_left)
        node_2.SetSolutionStepValue(var_disp, 0, disp_right)
        node_1.Fix(var_disp)
        node_2.Fix(var_disp)

        # set rotation
        rotation_left, rotation_right = 1.0, 0.0

        if out_of_plane_axis == 0:
            var_rot = KratosMultiphysics.ROTATION_X
        elif out_of_plane_axis == 1:
            var_rot = KratosMultiphysics.ROTATION_Y
        elif out_of_plane_axis == 2:
            var_rot = KratosMultiphysics.ROTATION_Z
        else:
            raise Exception("Wrong out of plane axis, must be 0, 1 or 2")

        node_1.SetSolutionStepValue(var_rot, 0, rotation_left)
        node_2.SetSolutionStepValue(var_rot, 0, rotation_right)
        node_1.Fix(var_disp)
        node_2.Fix(var_disp)

        location_load = 0.5

        # calculate expected deflection and rotation at the location of the load
        expected_deflection_axis_1, expected_deflection_axis_2 = (
            self.calculate_deflection_of_beam(disp_left, disp_right, rotation_left, rotation_right, length,
                                              location_load, rotation))

        expected_deflection = [0, 0, 0]
        expected_deflection[axis_1] = expected_deflection_axis_1
        expected_deflection[axis_2] = expected_deflection_axis_2

        expected_rotation_plane = self.calculate_rotation_of_beam(disp_left, disp_right, rotation_left, rotation_right,
                                                                  length,location_load, rotation)
        expected_rotation = [0, 0, 0]
        expected_rotation[out_of_plane_axis] = expected_rotation_plane
        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition3D2N", 1, [1, 2], mp.GetProperties()[1])

        cond.SetValue(StructuralMechanicsApplication.POINT_LOAD, [0, 1, 0])
        cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, location_load)

        # solve system and calculate displacement and rotation at the location of the load
        strategy.InitializeSolutionStep()
        strategy.SolveSolutionStep()

        # get calculated displacement at the load position
        disp = cond.GetValue(KratosMultiphysics.DISPLACEMENT)

        self.assertVectorAlmostEqual(disp, expected_deflection)

        # get calculated rotation at the load position
        rot = cond.GetValue(KratosMultiphysics.ROTATION)
        self.assertVectorAlmostEqual(rot, expected_rotation)

    def _MovingLoadCondition3D3N_disp_rot_at_load(self, second_coordinates, axis_1, axis_2):

        # get out of plane axis
        out_of_plane_axis = list({0, 1, 2} - {axis_1, axis_2})[0]

        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)

        strategy = self.setup_strategy(mp)

        # create nodes

        node_1 = mp.CreateNewNode(1, 0.0, 0.0, 0.0)
        node_2 = mp.CreateNewNode(2, second_coordinates[0], second_coordinates[1], second_coordinates[2])
        node_3 = mp.CreateNewNode(3, second_coordinates[0]/2, second_coordinates[1]/2, second_coordinates[2]/2)

        length = math.sqrt(second_coordinates[0] ** 2 + second_coordinates[1] ** 2 + second_coordinates[2] ** 2)

        # add dofs
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, mp)

        # set displacement
        disp_left, disp_right = 2, 0.0

        if axis_2 == 0:
            var_disp = KratosMultiphysics.DISPLACEMENT_X
        elif axis_2 == 1:
            var_disp = KratosMultiphysics.DISPLACEMENT_Y
        elif axis_2 == 2:
            var_disp = KratosMultiphysics.DISPLACEMENT_Z
        else:
            raise Exception("Wrong axis_2, must be 0, 1 or 2")

        node_1.SetSolutionStepValue(var_disp, 0, disp_left)
        node_2.SetSolutionStepValue(var_disp, 0, disp_right)
        node_3.SetSolutionStepValue(var_disp, 0, (disp_left + disp_right)/2)
        node_1.Fix(var_disp)
        node_2.Fix(var_disp)
        node_3.Fix(var_disp)

        location_load = 0.5

        # calculate expected deflection and rotation at the location of the load
        expected_deflection = [0, 0, 0]
        expected_deflection[axis_2] = disp_left * (1 - location_load / length) + disp_right * location_load / length

        expected_rotation = (disp_right - disp_left)/length

        cond = mp.CreateNewCondition("MovingLoadCondition3D3N", 1, [1, 2, 3], mp.GetProperties()[1])

        cond.SetValue(StructuralMechanicsApplication.POINT_LOAD, [0, 1, 0])
        cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, location_load)

        # solve system and calculate displacement and rotation at the location of the load
        strategy.InitializeSolutionStep()
        strategy.SolveSolutionStep()

        # get calculated displacement at the load position
        disp = cond.GetValue(KratosMultiphysics.DISPLACEMENT)
        rot = cond.GetValue(KratosMultiphysics.ROTATION)
        self.assertVectorAlmostEqual(disp, expected_deflection)
        self.assertAlmostEqual(rot[out_of_plane_axis], expected_rotation)

    def _MovingLoadCondition2D3N(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        strategy = self.setup_strategy(mp)

        # create nodes
        second_coord = [math.sqrt(2), math.sqrt(2.0), 0.0]
        third_coord = [math.sqrt(2)/2, math.sqrt(2.0)/2, 0.0]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])
        length = 2.0
        rotation = math.atan(second_coord[1] / second_coord[0])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition2D3N", 1, [1,2,3], mp.GetProperties()[1])

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # set perpendicular POINT_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 2.0
        load_on_cond[1] = -3.0
        load_on_cond[2] = 0.0 # note that this is a 2D condition

        num_tests = 5
        for i in range(num_tests):
            distance = i / (num_tests - 1) * length

            (x_left, y_left),(x_mid, y_mid), (x_right, y_right) = \
                self.__CalculateReaction3n(length, rotation, load_on_cond, distance)

            cond.SetValue(StructuralMechanicsApplication.POINT_LOAD,load_on_cond)

            # set load on node
            cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, distance)
            strategy.InitializeSolutionStep()
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

            self.assertAlmostEqual(rhs[0], x_left)
            self.assertAlmostEqual(rhs[1], y_left)
            self.assertAlmostEqual(rhs[2], x_right)
            self.assertAlmostEqual(rhs[3], y_right)
            self.assertAlmostEqual(rhs[4], x_mid)
            self.assertAlmostEqual(rhs[5], y_mid)

        # set load outside condition
        cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, length * 2)
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], 0)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)
        self.assertAlmostEqual(rhs[4], 0)
        self.assertAlmostEqual(rhs[5], 0)

    def _MovingLoadCondition3D3N(self):
        current_model = KratosMultiphysics.Model()
        mp = current_model.CreateModelPart("solid_part")
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        strategy = self.setup_strategy(mp)
        # create nodes
        second_coord = [math.sqrt(2), 0.0, math.sqrt(2.0)]
        third_coord = [math.sqrt(2)/2, 0.0, math.sqrt(2.0)/2]
        mp.CreateNewNode(1,0.0,0.0,0.0)
        mp.CreateNewNode(2,second_coord[0],second_coord[1],second_coord[2])
        mp.CreateNewNode(3, third_coord[0], third_coord[1], third_coord[2])
        length = 2.0
        rotation = math.atan(second_coord[2] / second_coord[0])

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)

        # create condition
        cond = mp.CreateNewCondition("MovingLoadCondition3D3N", 1, [1,2,3], mp.GetProperties()[1])

        # initialise matrices
        lhs = KratosMultiphysics.Matrix(0,0)
        rhs = KratosMultiphysics.Vector(0)

        # set perpendicular POINT_LOAD to the condition
        load_on_cond = KratosMultiphysics.Vector(3)
        load_on_cond[0] = 2.0
        load_on_cond[1] = 0.0
        load_on_cond[2] = -3.0

        num_tests = 5
        for i in range(num_tests):
            distance = i / (num_tests - 1) * length

            (x_left, z_left), (x_mid, z_mid), (x_right, z_right) = \
                self.__CalculateReaction3n(length, rotation, [load_on_cond[0], load_on_cond[2], load_on_cond[1]],
                                           distance)

            cond.SetValue(StructuralMechanicsApplication.POINT_LOAD,load_on_cond)

            # set load on node
            cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, distance)
            strategy.InitializeSolutionStep()
            cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)

            self.assertAlmostEqual(rhs[0], x_left)
            self.assertAlmostEqual(rhs[1], 0)
            self.assertAlmostEqual(rhs[2], z_left)
            self.assertAlmostEqual(rhs[3], x_right)
            self.assertAlmostEqual(rhs[4], 0)
            self.assertAlmostEqual(rhs[5], z_right)
            self.assertAlmostEqual(rhs[6], x_mid)
            self.assertAlmostEqual(rhs[7], 0)
            self.assertAlmostEqual(rhs[8], z_mid)

        # set load outside condition
        cond.SetValue(StructuralMechanicsApplication.MOVING_LOAD_LOCAL_DISTANCE, length * 2)
        strategy.InitializeSolutionStep()
        cond.CalculateLocalSystem(lhs, rhs, mp.ProcessInfo)
        self.assertAlmostEqual(rhs[0], 0)
        self.assertAlmostEqual(rhs[1], 0)
        self.assertAlmostEqual(rhs[2], 0)
        self.assertAlmostEqual(rhs[3], 0)
        self.assertAlmostEqual(rhs[4], 0)
        self.assertAlmostEqual(rhs[5], 0)
        self.assertAlmostEqual(rhs[6], 0)
        self.assertAlmostEqual(rhs[7], 0)
        self.assertAlmostEqual(rhs[8], 0)

    def test_MovingLoadCondition3D2NRotDofZ(self):
        self._MovingLoadCondition3D2NRotDofZ()

    def test_MovingLoadCondition3D2NRotDofY(self):
        self._MovingLoadCondition3D2NRotDofY()

    def test_MovingLoadCondition2D2N(self):
        self._MovingLoadCondition2D2N()
    def test_MovingLoadCondition2D2NRot(self):
        self._MovingLoadCondition2D2NRot()

    def test_MovingLoadCondition2D3N(self):
        self._MovingLoadCondition2D3N()

    def test_MovingLoadCondition3D3N(self):
        self._MovingLoadCondition3D3N()

    def test_MovingLoadCondition2D2NRot_disp_rot_at_load(self):
        self._MovingLoadCondition2D2NRot_disp_rot_at_load()

    def test_MovingLoadCondition2D2NImposedRot_disp_rot_at_load(self):
        """
        Tests the MovingLoadCondition2D2N with imposed rotation and checks the displacement and rotation at the position
        of the load
        """
        self._MovingLoadCondition2D2NImposedRot_disp_rot_at_load()

    def test_MovingLoadCondition3D2NRot_disp_rot_at_load_x_y_plane(self):
        coords_node_2 = [math.sqrt(2.0), math.sqrt(2.0), 0.0]
        self._MovingLoadCondition3D2NRot_disp_rot_at_load(coords_node_2, 0,1)

    def test_MovingLoadCondition3D2NRot_disp_rot_at_load_x_z_plane(self):
        coords_node_2 = [math.sqrt(2.0), 0.0, math.sqrt(2.0)]
        self._MovingLoadCondition3D2NRot_disp_rot_at_load(coords_node_2, 0,2)

    def test_MovingLoadCondition3D2NRot_disp_rot_at_load_z_y_plane(self):
        coords_node_2 = [0.0, math.sqrt(2.0), math.sqrt(2.0)]
        self._MovingLoadCondition3D2NRot_disp_rot_at_load(coords_node_2, 2,1)

    def test_MovingLoadCondition3D2N_disp_rot_at_load_x_y_plane(self):
        coords_node_2 = [2, 3, 0]
        self._MovingLoadCondition3D3N_disp_rot_at_load(coords_node_2, 0, 1)

    def test_MovingLoadCondition3D2N_disp_rot_at_load_z_y_plane(self):
        coords_node_2 = [0, 3, 2]
        self._MovingLoadCondition3D3N_disp_rot_at_load(coords_node_2, 2, 1)

if __name__ == '__main__':
    KratosUnittest.main()
